#include "DroneTrajectory.h"

#include <cmath>
#include <Eigen/Dense>

double PIDctrl(PIDParameters params, Eigen::Vector<double, NUM_PID_STATES> state);
Eigen::Vector<double, NUM_PID_STATES> updatePIDstate(Eigen::Vector<double, NUM_PID_STATES> currVal, double currSig, double prevSig, double ref, double timestep, double time);
std::string printPIDstate(Eigen::Vector<double, NUM_PID_STATES> state);
double capAngle(double angle);
double rad2Deg(double angle) { return angle / M_PI * 180.0; }
double deg2Rad(double angle) { return angle / 180.0 * M_PI; }
// These are different because during the simulation we might not know if its converging or not converging
// aka we don't know what it's doing yet, in which case both of these functions will return false

DroneTrajectory::DroneTrajectory(
    Logger & logger, 
    std::array<double(*)(double), NUM_DIST_STATES> const& dist,
    std::array<double(*)(double), NUM_REF_STATES> const& ref,
    std::array<PIDParameters, NUM_PIDS> ctrlParams,
    DroneParameters droneParameters, 
    double simTimestep, double finalTime, double sampleRate, double cutoffFreq, bool fixedNumIterations) : 
    m_logger(logger), 
    m_ctrlParams(ctrlParams),
    m_droneParams(droneParameters), 
    m_dist(dist), 
    m_ref(ref), 
    m_simTimestep(simTimestep),
    m_finalTime(finalTime),
    m_fixedNumIterations(fixedNumIterations) {

        double fr = sampleRate/cutoffFreq;
        double ohm = std::tan(M_PI/fr);
        double c = 1.0+2.0*std::cos(M_PI/4.0)*ohm+ohm*ohm;
        lpf_b0 = ohm*ohm/c;
        lpf_b1 = 2.0*lpf_b0;
        lpf_b2 = lpf_b0;
        lpf_a1 = 2.0*(ohm*ohm-1.0f)/c;
        lpf_a2 = (1.0-2.0*std::cos(M_PI/4.0)*ohm+ohm*ohm)/c;
       
    }

SimResults DroneTrajectory::Trajectory(SystemState initialState)
{
    SimResults simResults;
    simResults.stateProgression.push_back(initialState);
    double time = 0;
    simResults.time.push_back(time);
    while(time <= m_finalTime && simResults.stable){
        SystemState prev = simResults.stateProgression.back();
        Timestep state1 = simulateTimestep(prev, time, m_simTimestep);
        simResults.stateProgression.push_back(state1.state);
        simResults.stable &= state1.stable;

        time += m_simTimestep;        
        simResults.time.push_back(time);

        if (!m_fixedNumIterations)
        {
            Timestep state2 = simulateTimestep(state1.state, time, m_simTimestep);
            Timestep stateDouble = simulateTimestep(prev, time, 2*m_simTimestep);
            simResults.stateProgression.push_back(state2.state);
            simResults.stable &= state2.stable;

            time += m_simTimestep;        
            simResults.time.push_back(time);

            if ((stateDouble.state.plant - state2.state.plant).norm() < 1e-7 && stateDouble.stable) {
                m_simTimestep *= 2;
            } else if ((stateDouble.state.plant - state2.state.plant).norm() > 1e-4) {
                if (m_simTimestep > 1e-3){
                    m_simTimestep /= 2;
                }
            }
            
            // TOOD: check not converging and end earlier
            if (isConverging(state2.state, m_ref, time)){
                // m_logger << "converging :D" << std::endl;
                simResults.converged = true;
                return simResults;
            } else if (isNotConverging(state2.state, m_ref, time)) {
                // m_logger << "not converging D:" << std::endl;
                return simResults;
            }
        }
    }
    // m_logger << "convergence status: " << simResults.converged << std::endl;
    return simResults;
}

Timestep DroneTrajectory::simulateTimestep(SystemState prev, double time, double timestep)
{
    double tol = 1e-12;
    SystemState guess = prev;
    guess.plant += timestep*f(prev, time-timestep);
    int count = 0;
    int max_iterations = 100;
    bool stable = true;
    for(; count < max_iterations; count++ ){
        Eigen::Vector<double, NUM_PLANT_STATES> h = H(prev, guess.plant, time, timestep);
        if (h.norm() < tol) {
            break;
        }
        Eigen::MatrixX<double> dh = DH(guess, timestep);
        if (dh.determinant() != 0) {
            guess.plant = guess.plant - dh.inverse()*h;
        } else {
            stable = false;
            m_logger.warn(std::string("determinant is zero")); 
        }
    }

    Eigen::Vector<double, NUM_ALGE_STATES> algeStates = CascadedPIDController(guess.plant, prev.plant, guess.alge, time, timestep);

    guess.alge = algeStates;
    stable = stable && count < max_iterations;   
    if (count >= max_iterations) {
        m_logger << "WARN: took more than " << max_iterations << " iterations for sim to converge" << std::endl;
    }
    return {guess, stable};
}

Eigen::Vector<double, NUM_PLANT_STATES> DroneTrajectory::H(SystemState prev, Eigen::Vector<double, NUM_PLANT_STATES> guess, double time, double timestep)
{
    Eigen::Vector<double, NUM_PLANT_STATES> fprev = f(prev, time-timestep);
    Eigen::Vector<double, NUM_PLANT_STATES> fguess = f({guess, prev.alge}, time);
    return prev.plant - guess + timestep/2*(fprev + fguess);
}

Eigen::MatrixX<double> DroneTrajectory::DH(SystemState state, double timestep)
{
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dh = -1*Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>::Identity();
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus = dfdx(state);
    return dh + timestep/2*dfdx_plus;
}

// dynamics and updated control state
Eigen::Vector<double, NUM_PLANT_STATES> DroneTrajectory::f(SystemState state, double time)
{
    Eigen::Vector<double, NUM_PLANT_STATES> dot = Eigen::Vector<double, NUM_PLANT_STATES>::Zero();
    dot(x) = state.plant(xdot);
    dot(y) = state.plant(ydot);
    dot(z) = state.plant(zdot);

    dot(phi) = state.plant(p) 
                + state.plant(r) * std::cos(state.plant(phi)) * std::tan(state.plant(theta))
                + state.plant(q) * std::sin(state.plant(phi)) * std::tan(state.plant(theta));
    dot(theta) = state.plant(q) * std::cos(state.plant(phi)) 
                - state.plant(r) * std::sin(state.plant(phi));
    dot(psi) = state.plant(r) * std::cos(state.plant(phi)) / std::cos(state.plant(theta)) 
                + state.plant(q) * std::sin(state.plant(phi)) / std::cos(state.plant(theta));

    dot(xdot) = (state.alge(ft) + m_dist.at(Fwx)(time)) / m_droneParams.mass * 
                (std::sin(state.plant(phi)) * std::sin(state.plant(psi)) 
               + std::cos(state.plant(phi))*std::cos(state.plant(psi))*std::sin(state.plant(theta)));

    dot(ydot) = (state.alge(ft) + m_dist.at(Fwy)(time)) / m_droneParams.mass * 
                (std::cos(state.plant(phi))*std::sin(state.plant(psi))*std::sin(state.plant(theta)) 
               - std::cos(state.plant(psi))*std::sin(state.plant(phi)));

    dot(zdot) = - m_droneParams.g + (state.alge(ft) + m_dist.at(Fwz)(time))/m_droneParams.mass * 
                (std::cos(state.plant(phi)) * std::cos(state.plant(theta)));

    dot(p) = (m_droneParams.Iy-m_droneParams.Iz)/m_droneParams.Ix*state.plant(r)*state.plant(q) 
            + (state.alge(tx) + m_dist.at(Twx)(time))/m_droneParams.Ix;
    dot(q) = (m_droneParams.Iz-m_droneParams.Ix)/m_droneParams.Iy*state.plant(p)*state.plant(r)
            + (state.alge(ty) + m_dist.at(Twy)(time))/m_droneParams.Iy;
    dot(r) = (m_droneParams.Ix-m_droneParams.Iy)/m_droneParams.Iz*state.plant(p)*state.plant(q)
            + (state.alge(tz) + m_dist.at(Twz)(time))/m_droneParams.Iz;
    return dot;
}


// https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L193
Eigen::Vector<double, NUM_ALGE_STATES> DroneTrajectory::CascadedPIDController( 
    Eigen::Vector<double, NUM_PLANT_STATES> plantState,
    Eigen::Vector<double, NUM_PLANT_STATES> prevPlantState,
    Eigen::Vector<double, NUM_ALGE_STATES> currAlgeStates,
    double time, double timestep)
{
    Eigen::Vector<double, NUM_ALGE_STATES> algeStates;

    //https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L202
    double cosyaw = std::cos(plantState(psi));
    double sinyaw = std::sin(plantState(psi));
    double prevcosyaw = std::cos(prevPlantState(psi));
    double prevsinyaw = std::sin(prevPlantState(psi));

    // Convert global x y z into body x y z assuming CCW yaw = +ve
    algeStates(eix) = currAlgeStates(eix) + timestep*((m_ref.at(refx)(time) - plantState(x))*cosyaw + (m_ref.at(refy)(time) - plantState(y))*sinyaw);
    algeStates(eiy) = currAlgeStates(eiy) + timestep*((m_ref.at(refx)(time) - plantState(x))*(-sinyaw) + (m_ref.at(refy)(time) - plantState(y))*cosyaw);
    algeStates(eiz) = currAlgeStates(eiz) + timestep*(m_ref.at(refz)(time) - plantState(z));

    algeStates(edx) = -1/timestep*( plantState(x)*cosyaw - prevPlantState(x)*prevcosyaw + plantState(y)*sinyaw - prevPlantState(y)*prevsinyaw);
    algeStates(edy) = -1/timestep*(-plantState(x)*sinyaw + prevPlantState(x)*prevsinyaw + plantState(y)*cosyaw - prevPlantState(y)*prevcosyaw);
    algeStates(edz) = -1/timestep*(plantState(z) - prevPlantState(z));
    
    algeStates(desVelX) = PIDctrl(m_ctrlParams.at(posX), {(m_ref.at(refx)(time) - plantState(x))*cosyaw + (m_ref.at(refy)(time) - plantState(y))*sinyaw, algeStates(eix), algeStates(edx)});
    algeStates(desVelY) = PIDctrl(m_ctrlParams.at(posY), {(m_ref.at(refx)(time) - plantState(x))*(-sinyaw) + (m_ref.at(refy)(time) - plantState(y))*cosyaw, algeStates(eiy), algeStates(edy)});
    algeStates(desVelZ) = PIDctrl(m_ctrlParams.at(posZ), {m_ref.at(refz)(time) - plantState(z), algeStates(eiz), algeStates(edz)});

    algeStates(eixdot) = currAlgeStates(eixdot) + timestep*(algeStates(desVelX) - plantState(xdot)*cosyaw - plantState(ydot)*sinyaw);
    algeStates(eiydot) = currAlgeStates(eiydot) + timestep*(algeStates(desVelY) + plantState(xdot)*sinyaw - plantState(ydot)*cosyaw);
    algeStates(eizdot) = currAlgeStates(eizdot) + timestep*(algeStates(desVelZ) - plantState(zdot));
    
    algeStates(edxdot) = -1/timestep*(plantState(xdot)*cosyaw + plantState(ydot)*sinyaw - prevPlantState(xdot)*prevcosyaw - prevPlantState(ydot)*prevsinyaw);
    algeStates(edydot) = -1/timestep*(-plantState(xdot)*sinyaw + plantState(ydot)*cosyaw + prevPlantState(xdot)*prevsinyaw - prevPlantState(ydot)*prevcosyaw);
    algeStates(edzdot) = -1/timestep*(plantState(zdot) - prevPlantState(zdot));
    
    algeStates(desThrust) = PIDctrl(m_ctrlParams.at(velZ), {algeStates(desVelZ) - plantState(zdot), algeStates(eizdot), algeStates(edzdot)})*m_thrustScale+m_thrustBase;
    
    // Need a positive pitch to move forward in the x direction
    // Need a negative roll to move forward in the y direction
    algeStates(desRoll) = -std::clamp(PIDctrl(m_ctrlParams.at(velY), {algeStates(desVelY) + plantState(xdot)*sinyaw - plantState(ydot)*cosyaw, algeStates(eiydot), algeStates(edydot)}), -m_droneParams.pid_vel_pitch_max,  m_droneParams.pid_vel_pitch_max);
    algeStates(desPitch) = std::clamp(PIDctrl(m_ctrlParams.at(velX), {algeStates(desVelX) - plantState(xdot)*cosyaw - plantState(ydot)*sinyaw, algeStates(eixdot), algeStates(edxdot)}), -m_droneParams.pid_vel_roll_max,  m_droneParams.pid_vel_roll_max);

    algeStates(eiphi) = currAlgeStates(eiphi) + timestep*(algeStates(desRoll) - 180/M_PI*plantState(phi));
    algeStates(eitheta) = currAlgeStates(eitheta) + timestep*(algeStates(desPitch) - 180/M_PI*plantState(theta));
    algeStates(eipsi) = currAlgeStates(eipsi) + timestep*(m_ref.at(refyaw)(time) - 180/M_PI*plantState(psi));

    algeStates(edphi) = -1/timestep*180/M_PI*(plantState(phi) - prevPlantState(phi));
    algeStates(edtheta) = -1/timestep*180/M_PI*(plantState(theta) - prevPlantState(theta));
    algeStates(edpsi) = -1/timestep*180/M_PI*(plantState(psi) - prevPlantState(psi));

    algeStates(desRollRate) = PIDctrl(m_ctrlParams.at(roll), {algeStates(desRoll) - 180/M_PI*plantState(phi), algeStates(eiphi), algeStates(edphi)});
    algeStates(desPitchRate) = PIDctrl(m_ctrlParams.at(pitch), {algeStates(desPitch) - 180/M_PI*plantState(theta), algeStates(eitheta), algeStates(edtheta)});
    algeStates(desYawRate) = PIDctrl(m_ctrlParams.at(yaw), {m_ref.at(refyaw)(time) - 180/M_PI*plantState(psi), algeStates(eipsi), algeStates(edpsi)});
    
    // LPF for rollRate and pitchRate
    algeStates(delay_1_rollRate) = -rad2Deg(plantState(p)-prevPlantState(p))/timestep - currAlgeStates(delay_1_rollRate)*lpf_a1 - currAlgeStates(delay_2_rollRate)*lpf_a2;
    algeStates(delay_1_pitchRate) = -rad2Deg(plantState(q)-prevPlantState(q))/timestep - currAlgeStates(delay_1_pitchRate)*lpf_a1 - currAlgeStates(delay_2_pitchRate)*lpf_a2;
    algeStates(delay_2_rollRate) = currAlgeStates(delay_1_rollRate);
    algeStates(delay_2_pitchRate) = currAlgeStates(delay_1_pitchRate);

    algeStates(eip) = currAlgeStates(eip) + timestep*(algeStates(desRollRate) - 180/M_PI*plantState(p));
    algeStates(eiq) = currAlgeStates(eiq) + timestep*(algeStates(desPitchRate) - 180/M_PI*plantState(q));
    algeStates(eir) = currAlgeStates(eir) + timestep*(algeStates(desYawRate) - 180/M_PI*plantState(r));

    algeStates(edp) = lpf_b0 * algeStates(delay_1_rollRate) + currAlgeStates(delay_1_rollRate) * lpf_b1 + currAlgeStates(delay_2_rollRate) * lpf_b2;
    algeStates(edq) = lpf_b0 * algeStates(delay_1_pitchRate) + currAlgeStates(delay_1_pitchRate) * lpf_b1 + currAlgeStates(delay_2_pitchRate) * lpf_b2;
    algeStates(edr) = -1/timestep*180/M_PI*(plantState(r) - prevPlantState(r));

    algeStates(desRollOutput)  = PIDctrl(m_ctrlParams.at(rollRate), {algeStates(desRollRate) - 180/M_PI*plantState(p), algeStates(eip), algeStates(edp)});
    algeStates(desPitchOutput) = PIDctrl(m_ctrlParams.at(pitchRate), {algeStates(desPitchRate) - 180/M_PI*plantState(q), algeStates(eiq), algeStates(edq)});
    algeStates(desYawOutput) = PIDctrl(m_ctrlParams.at(yawRate), {algeStates(desYawRate) - 180/M_PI*plantState(r), algeStates(eir), algeStates(edr)});

    // TODO: implement saturation

    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/power_distribution_quadrotor.c
    // line 86 they divide roll and pitch by 2 
    
    algeStates(w1) = m_alpha * (algeStates(desThrust) - 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput));
    algeStates(w2) = m_alpha * (algeStates(desThrust) - 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput));
    algeStates(w3) = m_alpha * (algeStates(desThrust) + 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput));
    algeStates(w4) = m_alpha * (algeStates(desThrust) + 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput));

    // https://giuseppesilano.net/publications/rosChapter19.pdf
    // the crazyflie firmware using x wing configuration 
    algeStates(ft) = m_droneParams.kf * (algeStates(w1)*algeStates(w1) + algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tx) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (- algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(ty) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (+ algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) - algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tz) = m_droneParams.km * (algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) - algeStates(w4)*algeStates(w4));

    return algeStates;
}

Eigen::Vector<double, NUM_ALGE_STATES> DroneTrajectory::h( 
    Eigen::Vector<double, NUM_PLANT_STATES> plantState,
    Eigen::Vector<double, NUM_PLANT_STATES> prevPlantState,
    Eigen::Vector<double, NUM_ALGE_STATES> currAlgeStates,
    double time, double timestep, int index, double delta)
{
    Eigen::Vector<double, NUM_ALGE_STATES> algeStates;

    //https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L202
    double cosyaw = std::cos(plantState(psi));
    double sinyaw = std::sin(plantState(psi));
    double prevcosyaw = std::cos(prevPlantState(psi));
    double prevsinyaw = std::sin(prevPlantState(psi));

    // Convert global x y z into body x y z assuming CCW yaw = +ve
    algeStates(eix) = currAlgeStates(eix) + timestep*((m_ref.at(refx)(time) - plantState(x))*cosyaw + (m_ref.at(refy)(time) - plantState(y))*sinyaw);
    algeStates(eiy) = currAlgeStates(eiy) + timestep*((m_ref.at(refx)(time) - plantState(x))*(-sinyaw) + (m_ref.at(refy)(time) - plantState(y))*cosyaw);
    algeStates(eiz) = currAlgeStates(eiz) + timestep*(m_ref.at(refz)(time) - plantState(z));

    if (index == eix) { algeStates(eix) += delta;} 
    if (index == eiy) { algeStates(eiy) += delta;} 
    if (index == eiz) { algeStates(eiz) += delta;} 

    algeStates(edx) = -1/timestep*( plantState(x)*cosyaw - prevPlantState(x)*prevcosyaw + plantState(y)*sinyaw - prevPlantState(y)*prevsinyaw);
    algeStates(edy) = -1/timestep*(-plantState(x)*sinyaw + prevPlantState(x)*prevsinyaw + plantState(y)*cosyaw - prevPlantState(y)*prevcosyaw);
    algeStates(edz) = -1/timestep*(plantState(z) - prevPlantState(z));

    if (index == edx) { algeStates(edx) += delta;} 
    if (index == edy) { algeStates(edy) += delta;} 
    if (index == edz) { algeStates(edz) += delta;}
    
    algeStates(desVelX) = PIDctrl(m_ctrlParams.at(posX), {(m_ref.at(refx)(time) - plantState(x))*cosyaw + (m_ref.at(refy)(time) - plantState(y))*sinyaw, algeStates(eix), algeStates(edx)});
    algeStates(desVelY) = PIDctrl(m_ctrlParams.at(posY), {(m_ref.at(refx)(time) - plantState(x))*(-sinyaw) + (m_ref.at(refy)(time) - plantState(y))*cosyaw, algeStates(eiy), algeStates(edy)});
    algeStates(desVelZ) = PIDctrl(m_ctrlParams.at(posZ), {m_ref.at(refz)(time) - plantState(z), algeStates(eiz), algeStates(edz)});

    if (index == desVelX) { algeStates(desVelX) += delta;} 
    if (index == desVelY) { algeStates(desVelY) += delta;} 
    if (index == desVelZ) { algeStates(desVelZ) += delta;}

    algeStates(eixdot) = currAlgeStates(eixdot) + timestep*(algeStates(desVelX) - plantState(xdot)*cosyaw - plantState(ydot)*sinyaw);
    algeStates(eiydot) = currAlgeStates(eiydot) + timestep*(algeStates(desVelY) + plantState(xdot)*sinyaw - plantState(ydot)*cosyaw);
    algeStates(eizdot) = currAlgeStates(eizdot) + timestep*(algeStates(desVelZ) - plantState(zdot));

    if (index == eixdot) { algeStates(eixdot) += delta;} 
    if (index == eiydot) { algeStates(eiydot) += delta;} 
    if (index == eizdot) { algeStates(eizdot) += delta;}
    
    algeStates(edxdot) = -1/timestep*(plantState(xdot)*cosyaw + plantState(ydot)*sinyaw - prevPlantState(xdot)*prevcosyaw - prevPlantState(ydot)*prevsinyaw);
    algeStates(edydot) = -1/timestep*(-plantState(xdot)*sinyaw + plantState(ydot)*cosyaw + prevPlantState(xdot)*prevsinyaw - prevPlantState(ydot)*prevcosyaw);
    algeStates(edzdot) = -1/timestep*(plantState(zdot) - prevPlantState(zdot));

    if (index == edxdot) { algeStates(edxdot) += delta;} 
    if (index == edydot) { algeStates(edydot) += delta;} 
    if (index == edzdot) { algeStates(edzdot) += delta;}
    
    algeStates(desThrust) = PIDctrl(m_ctrlParams.at(velZ), {algeStates(desVelZ) - plantState(zdot), algeStates(eizdot), algeStates(edzdot)})*m_thrustScale+m_thrustBase;
    algeStates(desRoll) = -std::clamp(PIDctrl(m_ctrlParams.at(velY), {algeStates(desVelY) + plantState(xdot)*sinyaw - plantState(ydot)*cosyaw, algeStates(eiydot), algeStates(edydot)}), -m_droneParams.pid_vel_pitch_max,  m_droneParams.pid_vel_pitch_max);
    algeStates(desPitch) = std::clamp(PIDctrl(m_ctrlParams.at(velX), {algeStates(desVelX) - plantState(xdot)*cosyaw - plantState(ydot)*sinyaw, algeStates(eixdot), algeStates(edxdot)}), -m_droneParams.pid_vel_roll_max,  m_droneParams.pid_vel_roll_max);

    if (index == desThrust) { algeStates(desThrust) += delta;} 
    if (index == desRoll && algeStates(desRoll) > -20 && algeStates(desRoll) < 20 ) { algeStates(desRoll) += delta;} 
    if (index == desPitch  && algeStates(desPitch) > -20 && algeStates(desPitch) < 20 ) { algeStates(desPitch) += delta;}

    algeStates(eiphi) = currAlgeStates(eiphi) + timestep*(algeStates(desRoll) - 180/M_PI*plantState(phi));
    algeStates(eitheta) = currAlgeStates(eitheta) + timestep*(algeStates(desPitch) - 180/M_PI*plantState(theta));
    algeStates(eipsi) = currAlgeStates(eipsi) + timestep*(m_ref.at(refyaw)(time) - 180/M_PI*plantState(psi));

    if (index == eiphi) { algeStates(eiphi) += delta;} 
    if (index == eitheta) { algeStates(eitheta) += delta;} 
    if (index == eipsi) { algeStates(eipsi) += delta;}

    algeStates(edphi) = -1/timestep*180/M_PI*(plantState(phi) - prevPlantState(phi));
    algeStates(edtheta) = -1/timestep*180/M_PI*(plantState(theta) - prevPlantState(theta));
    algeStates(edpsi) = -1/timestep*180/M_PI*(plantState(psi) - prevPlantState(psi));

    if (index == edphi) { algeStates(edphi) += delta;} 
    if (index == edtheta) { algeStates(edtheta) += delta;} 
    if (index == edpsi) { algeStates(edpsi) += delta;}

    algeStates(desRollRate) = PIDctrl(m_ctrlParams.at(roll), {algeStates(desRoll) - 180/M_PI*plantState(phi), algeStates(eiphi), algeStates(edphi)});
    algeStates(desPitchRate) = PIDctrl(m_ctrlParams.at(pitch), {algeStates(desPitch) - 180/M_PI*plantState(theta), algeStates(eitheta), algeStates(edtheta)});
    algeStates(desYawRate) = PIDctrl(m_ctrlParams.at(yaw), {m_ref.at(refyaw)(time) - 180/M_PI*plantState(psi), algeStates(eipsi), algeStates(edpsi)});
    
    if (index == desRollRate) { algeStates(desRollRate) += delta;} 
    if (index == desPitchRate) { algeStates(desPitchRate) += delta;} 
    if (index == desYawRate) { algeStates(desYawRate) += delta;}

    // LPF for rollRate and pitchRate
    algeStates(delay_1_rollRate) = -180/M_PI*(plantState(p)-prevPlantState(p))/timestep - currAlgeStates(delay_1_rollRate)*lpf_a1 - currAlgeStates(delay_2_rollRate)*lpf_a2;
    algeStates(delay_1_pitchRate) = -180/M_PI*(plantState(q)-prevPlantState(q))/timestep - currAlgeStates(delay_1_pitchRate)*lpf_a1 - currAlgeStates(delay_2_pitchRate)*lpf_a2;
    algeStates(delay_2_rollRate) = currAlgeStates(delay_1_rollRate);
    algeStates(delay_2_pitchRate) = currAlgeStates(delay_1_pitchRate);

    if (index == delay_1_rollRate) { algeStates(delay_1_rollRate) += delta;} 
    if (index == delay_1_pitchRate) { algeStates(delay_1_pitchRate) += delta;} 
    if (index == delay_2_rollRate) { algeStates(delay_2_rollRate) += delta;}
    if (index == delay_2_pitchRate) { algeStates(delay_2_pitchRate) += delta;}

    algeStates(eip) = currAlgeStates(eip) + timestep*(algeStates(desRollRate) - 180/M_PI*plantState(p));
    algeStates(eiq) = currAlgeStates(eiq) + timestep*(algeStates(desPitchRate) - 180/M_PI*plantState(q));
    algeStates(eir) = currAlgeStates(eir) + timestep*(algeStates(desYawRate) - 180/M_PI*plantState(r));

    if (index == eip) { algeStates(eip) += delta;} 
    if (index == eiq) { algeStates(eiq) += delta;} 
    if (index == eir) { algeStates(eir) += delta;}

    algeStates(edp) = lpf_b0 * algeStates(delay_1_rollRate) + currAlgeStates(delay_1_rollRate) * lpf_b1 + currAlgeStates(delay_2_rollRate) * lpf_b2;
    algeStates(edq) = lpf_b0 * algeStates(delay_1_pitchRate) + currAlgeStates(delay_1_pitchRate) * lpf_b1 + currAlgeStates(delay_2_pitchRate) * lpf_b2;
    algeStates(edr) = -1/timestep*180/M_PI*(plantState(r) - prevPlantState(r));

    if (index == edp) { algeStates(edp) += delta;} 
    if (index == edq) { algeStates(edq) += delta;} 
    if (index == edr) { algeStates(edr) += delta;}

    algeStates(desRollOutput)  = PIDctrl(m_ctrlParams.at(rollRate), {algeStates(desRollRate) - 180/M_PI*plantState(p), algeStates(eip), algeStates(edp)});
    algeStates(desPitchOutput) = PIDctrl(m_ctrlParams.at(pitchRate), {algeStates(desPitchRate) - 180/M_PI*plantState(q), algeStates(eiq), algeStates(edq)});
    algeStates(desYawOutput) = PIDctrl(m_ctrlParams.at(yawRate), {algeStates(desYawRate) - 180/M_PI*plantState(r), algeStates(eir), algeStates(edr)});

    if (index == desRollOutput) { algeStates(desRollOutput) += delta;} 
    if (index == desPitchOutput) { algeStates(desPitchOutput) += delta;} 
    if (index == desYawOutput) { algeStates(desYawOutput) += delta;}

    algeStates(w1) = m_alpha * (algeStates(desThrust) - 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput));
    algeStates(w2) = m_alpha * (algeStates(desThrust) - 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput));
    algeStates(w3) = m_alpha * (algeStates(desThrust) + 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput));
    algeStates(w4) = m_alpha * (algeStates(desThrust) + 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput));

    if (index == w1) { algeStates(w1) += delta;} 
    if (index == w2) { algeStates(w2) += delta;} 
    if (index == w3) { algeStates(w3) += delta;}
    if (index == w4) { algeStates(w4) += delta;}

    algeStates(ft) = m_droneParams.kf * (algeStates(w1)*algeStates(w1) + algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tx) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (- algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(ty) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (+ algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) - algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tz) = m_droneParams.km * (algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) - algeStates(w4)*algeStates(w4));

    if (index == ft) { algeStates(ft) += delta;} 
    if (index == tx) { algeStates(tx) += delta;} 
    if (index == ty) { algeStates(ty) += delta;}
    if (index == tz) { algeStates(tz) += delta;}

    return algeStates;
}


double PIDctrl(PIDParameters params, Eigen::Vector<double, NUM_PID_STATES> state)
{ return params.ki * state(ki_error) + params.kp * state(kp_error) + params.kd * state(kd_error); }

Eigen::Vector<double, NUM_PID_STATES> updatePIDstate(Eigen::Vector<double, NUM_PID_STATES> currVal, double currSig, double prevSig, double ref, double timestep, double time)
{
    Eigen::Vector<double, NUM_PID_STATES> newPIDstate = Eigen::Vector<double, NUM_PID_STATES>::Zero();;
    newPIDstate(kp_error) = ref - currSig;
    newPIDstate(ki_error) = currVal(ki_error) + newPIDstate(kp_error)*timestep; 
   
    if (time != 0) {
        // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/utils/src/pid.c
        // prevent derivative kick
        newPIDstate(kd_error) = -1.0/timestep*(currSig - prevSig);     
    }
   
    return newPIDstate;
}

std::string printPIDstate(Eigen::Vector<double, NUM_PID_STATES> state)
{
    return  std::string("kp: ") + std::to_string(state(kp_error)) + 
            std::string(" ki: ") + std::to_string(state(ki_error)) +
            std::string(" kd: ") + std::to_string(state(kd_error));
}

double capAngle(double angle)
{
  double result = angle;

  while (result > 180.0) {
    result -= 2*180.0;
  }

  while (result < -180.0) {
    result += 2*180.0;
  }

  return result;
}

bool DroneTrajectory::isConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time)
{
    double posTol = 5e-1;
    double angleTol = 1e-3;
    double velTol = 1e-2;

    double ftTol = 1e-2;
    double hover = m_droneParams.mass * m_droneParams.g;
    double torqueTol = 1e-6;
    return  std::abs(state.plant(x) - ref.at(refx)(time)) < posTol && std::abs(state.plant(y) - ref.at(refy)(time)) < posTol &&  std::abs(state.plant(z) - ref.at(refz)(time)) < posTol && 
            // Assuming the reference is a fixed set point - if it becomes a trajectory, then we need to calculate how much the angle and velocities should be
            std::abs(state.plant(phi)) < angleTol && std::abs(state.plant(theta)) < angleTol &&  std::abs(state.plant(psi) - ref.at(refyaw)(time)) < angleTol && 
            std::abs(state.plant(xdot)) < velTol && std::abs(state.plant(ydot)) < velTol &&  std::abs(state.plant(zdot)) < velTol && 
            std::abs(state.plant(p)) < velTol && std::abs(state.plant(q)) < velTol &&  std::abs(state.plant(r)) < velTol &&
            std::abs(state.alge(ft) - hover) < ftTol && std::abs(state.alge(tx)) < torqueTol && std::abs(state.alge(ty)) < torqueTol && std::abs(state.alge(tz)) < torqueTol;
}   

bool DroneTrajectory::isNotConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time)
{
    double posDist = 100;
    double angleDist = 90;
    double pidStateLimit = 1e7;

    bool notConverging =std::abs(state.plant(x) - ref.at(refx)(time)) > posDist ||std::abs(state.plant(y) - ref.at(refy)(time)) > posDist || std::abs(state.plant(z) - ref.at(refz)(time)) > posDist ||
            // Assuming the reference is a fixed set point - if it becomes a trajectory, then we need to calculate how much the angle and velocities should be
           std::abs(state.plant(phi)) > deg2Rad(angleDist) ||std::abs(state.plant(theta)) > deg2Rad(angleDist);

    notConverging = notConverging ||std::abs(m_ctrlParams.at(posX).ki*state.alge(eix)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posX).kd*state.alge(edx)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(posY).ki*state.alge(eiy)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posY).kd*state.alge(edy)) > pidStateLimit;
    
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posZ).ki*state.alge(eiz)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posZ).kd*state.alge(edz)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velX).ki*state.alge(eixdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velX).kd*state.alge(edxdot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velY).ki*state.alge(eiydot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velY).kd*state.alge(edydot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velZ).ki*state.alge(eizdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velZ).kd*state.alge(edzdot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(roll).ki*state.alge(eiphi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(roll).kd*state.alge(edphi)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitch).ki*state.alge(eitheta)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitch).kd*state.alge(edtheta)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(yaw).ki*state.alge(eipsi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yaw).kd*state.alge(edpsi)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(rollRate).ki*state.alge(eip)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(rollRate).kd*state.alge(edp)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitchRate).ki*state.alge(eiq)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitchRate).kd*state.alge(edq)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(yawRate).ki*state.alge(eir)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yawRate).kd*state.alge(edr)) > pidStateLimit;
    
    return notConverging;
}

