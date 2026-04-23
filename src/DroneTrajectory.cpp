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

    // Convert global x y z into body x y z assuming CCW yaw = +ve
    algeStates(setp_body_x) = m_ref.at(refx)(time) * cosyaw + m_ref.at(refy)(time) * sinyaw;
    algeStates(setp_body_y) = -m_ref.at(refx)(time) * sinyaw + m_ref.at(refy)(time) * cosyaw;
    algeStates(state_body_x) = plantState(x) * cosyaw + plantState(y) * sinyaw;
    algeStates(state_body_y) = -plantState(x) * sinyaw + plantState(y) * cosyaw;    

    algeStates.segment(epx, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epx, NUM_PID_STATES), algeStates(state_body_x), currAlgeStates(state_body_x),  algeStates(setp_body_x), timestep, time);
    algeStates.segment(epy, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epy, NUM_PID_STATES), algeStates(state_body_y), currAlgeStates(state_body_y), algeStates(setp_body_y), timestep, time);
    algeStates.segment(epz, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epz, NUM_PID_STATES), plantState(z), prevPlantState(z), m_ref.at(refz)(time), timestep, time);
    
    algeStates(desVelX) = PIDctrl(m_ctrlParams.at(posX), algeStates.segment(epx, NUM_PID_STATES));
    algeStates(desVelY) = PIDctrl(m_ctrlParams.at(posY), algeStates.segment(epy, NUM_PID_STATES));
    algeStates(desVelZ) = PIDctrl(m_ctrlParams.at(posZ), algeStates.segment(epz, NUM_PID_STATES));
    
    algeStates(state_body_vx) = plantState(xdot) * cosyaw + plantState(ydot) * sinyaw;
    algeStates(state_body_vy) = -plantState(xdot) * sinyaw + plantState(ydot) * cosyaw;

    algeStates.segment(epxdot, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epxdot, NUM_PID_STATES), algeStates(state_body_vx), currAlgeStates(state_body_vx), algeStates(desVelX), timestep, time);
    algeStates.segment(epydot, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epydot, NUM_PID_STATES), algeStates(state_body_vy), currAlgeStates(state_body_vy), algeStates(desVelY), timestep, time);
    algeStates.segment(epzdot, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epzdot, NUM_PID_STATES), plantState(zdot), prevPlantState(zdot), algeStates(desVelZ), timestep, time);
 
    // Need a positive pitch to move forward in the x direction
    // Need a negative roll to move forward in the y direction
    algeStates(desRoll) = -std::clamp(PIDctrl(m_ctrlParams.at(velY), algeStates.segment(epydot, NUM_PID_STATES)), -m_droneParams.pid_vel_pitch_max,  m_droneParams.pid_vel_pitch_max);
    algeStates(desPitch) = std::clamp(PIDctrl(m_ctrlParams.at(velX), algeStates.segment(epxdot, NUM_PID_STATES)), -m_droneParams.pid_vel_roll_max,  m_droneParams.pid_vel_roll_max);
    // double desPitch = PIDctrl(m_ctrlParams.at(velX), newCtrlState.at(velX));
    // double desRoll = -PIDctrl(m_ctrlParams.at(velY), newCtrlState.at(velY));

    // TODO: Can add constraints later

    algeStates(desThrust) = PIDctrl(m_ctrlParams.at(velZ), algeStates.segment(epzdot, NUM_PID_STATES))*m_thrustScale+m_thrustBase;
    // desThrust = std::clamp(desThrust, 20000.0, (double)UINT16_MAX);

    algeStates.segment(epphi, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epphi, NUM_PID_STATES), rad2Deg(plantState(phi)), rad2Deg(prevPlantState(phi)), algeStates(desRoll), timestep, time);
    algeStates.segment(eptheta, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(eptheta, NUM_PID_STATES), rad2Deg(plantState(theta)), rad2Deg(prevPlantState(theta)), algeStates(desPitch), timestep, time);
    algeStates.segment(eppsi, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(eppsi, NUM_PID_STATES), rad2Deg(plantState(psi)), rad2Deg(prevPlantState(psi)), rad2Deg(m_ref.at(refyaw)(time)), timestep, time);
        
    algeStates(desRollRate) = PIDctrl(m_ctrlParams.at(roll), algeStates.segment(epphi, NUM_PID_STATES));
    algeStates(desPitchRate) = PIDctrl(m_ctrlParams.at(pitch), algeStates.segment(eptheta, NUM_PID_STATES));
    algeStates(desYawRate) = PIDctrl(m_ctrlParams.at(yaw), algeStates.segment(eppsi, NUM_PID_STATES));
    
    // Assume that [phi dot theta dot psi dot] = [p q r]
    // This assumption holds true for small angles of movement
    algeStates.segment(epp, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epp, NUM_PID_STATES), rad2Deg(plantState(p)), rad2Deg(prevPlantState(p)), algeStates(desRollRate), timestep, time);
    algeStates.segment(epq, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epq, NUM_PID_STATES), rad2Deg(plantState(q)), rad2Deg(prevPlantState(q)), algeStates(desPitchRate), timestep, time);
    algeStates.segment(epr, NUM_PID_STATES) = updatePIDstate(currAlgeStates.segment(epr, NUM_PID_STATES), rad2Deg(plantState(r)), rad2Deg(prevPlantState(r)), algeStates(desYawRate), timestep, time);

    // LPF for rollRate and pitchRate
    algeStates(delay_1_rollRate) = -rad2Deg(plantState(p)-prevPlantState(p))/timestep - currAlgeStates(delay_1_rollRate)*lpf_a1 - currAlgeStates(delay_2_rollRate)*lpf_a2;
    algeStates(delay_1_pitchRate) = -rad2Deg(plantState(q)-prevPlantState(q))/timestep - currAlgeStates(delay_1_pitchRate)*lpf_a1 - currAlgeStates(delay_2_pitchRate)*lpf_a2;
    algeStates(delay_2_rollRate) = currAlgeStates(delay_1_rollRate);
    algeStates(delay_2_pitchRate) = currAlgeStates(delay_1_pitchRate);

    algeStates(edp) = lpf_b0 * algeStates(delay_1_rollRate) + currAlgeStates(delay_1_rollRate) * lpf_b1 + currAlgeStates(delay_2_rollRate) * lpf_b2;
    algeStates(edq) = lpf_b0 * algeStates(delay_1_pitchRate) + currAlgeStates(delay_1_pitchRate) * lpf_b1 + currAlgeStates(delay_2_pitchRate) * lpf_b2;

    algeStates(desRollOutput)  = PIDctrl(m_ctrlParams.at(rollRate), algeStates.segment(epp, NUM_PID_STATES));
    algeStates(desPitchOutput) = PIDctrl(m_ctrlParams.at(pitchRate), algeStates.segment(epq, NUM_PID_STATES));
    algeStates(desYawOutput) = PIDctrl(m_ctrlParams.at(yawRate), algeStates.segment(epr, NUM_PID_STATES));

    // TODO: implement saturation

    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/power_distribution_quadrotor.c
    // line 86 they divide roll and pitch by 2 

    double m1 = algeStates(desThrust) - 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput);
    double m2 = algeStates(desThrust) - 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput);
    double m3 = algeStates(desThrust) + 0.5*algeStates(desRollOutput) - 0.5*algeStates(desPitchOutput) + algeStates(desYawOutput);
    double m4 = algeStates(desThrust) + 0.5*algeStates(desRollOutput) + 0.5*algeStates(desPitchOutput) - algeStates(desYawOutput);

    
    algeStates(w1) = m_alpha * m1;
    algeStates(w2) = m_alpha * m2;
    algeStates(w3) = m_alpha * m3;
    algeStates(w4) = m_alpha * m4;

    // https://giuseppesilano.net/publications/rosChapter19.pdf
    // the crazyflie firmware using x wing configuration 
    algeStates(ft) = m_droneParams.kf * (algeStates(w1)*algeStates(w1) + algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tx) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (- algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(ty) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (+ algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) - algeStates(w3)*algeStates(w3) + algeStates(w4)*algeStates(w4));
    algeStates(tz) = m_droneParams.km * (algeStates(w1)*algeStates(w1) - algeStates(w2)*algeStates(w2) + algeStates(w3)*algeStates(w3) - algeStates(w4)*algeStates(w4));

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

    notConverging = notConverging ||std::abs(m_ctrlParams.at(posX).kp*state.alge(epx)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posX).ki*state.alge(eix)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posX).kd*state.alge(edx)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(posY).kp*state.alge(epy)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posY).ki*state.alge(eiy)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posY).kd*state.alge(edy)) > pidStateLimit;
    
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posZ).kp*state.alge(epz)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posZ).ki*state.alge(eiz)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(posZ).kd*state.alge(edz)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velX).kp*state.alge(epxdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velX).ki*state.alge(eixdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velX).kd*state.alge(edxdot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velY).kp*state.alge(epydot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velY).ki*state.alge(eiydot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velY).kd*state.alge(edydot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(velZ).kp*state.alge(epzdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velZ).ki*state.alge(eizdot)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(velZ).kd*state.alge(edzdot)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(roll).kp*state.alge(epphi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(roll).ki*state.alge(eiphi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(roll).kd*state.alge(edphi)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitch).kp*state.alge(eptheta)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitch).ki*state.alge(eitheta)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitch).kd*state.alge(edtheta)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(yaw).kp*state.alge(eppsi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yaw).ki*state.alge(eipsi)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yaw).kd*state.alge(edpsi)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(rollRate).kp*state.alge(epp)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(rollRate).ki*state.alge(eip)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(rollRate).kd*state.alge(edp)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitchRate).kp*state.alge(epq)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitchRate).ki*state.alge(eiq)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(pitchRate).kd*state.alge(edq)) > pidStateLimit;

    notConverging = notConverging ||std::abs(m_ctrlParams.at(yawRate).kp*state.alge(epr)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yawRate).ki*state.alge(eir)) > pidStateLimit;
    notConverging = notConverging ||std::abs(m_ctrlParams.at(yawRate).kd*state.alge(edr)) > pidStateLimit;
    
    return notConverging;
}