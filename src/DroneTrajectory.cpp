#include "DroneTrajectory.h"

#include <cmath>
#include <Eigen/Dense>

double PIDctrl(PIDParameters params, PIDstate state);
PIDstate updatePIDstate(PIDstate currVal, double currSig, double ref, double timestep);

DroneTrajectory::DroneTrajectory(
    Logger & logger, 
    std::array<double(*)(double), NUM_DIST_STATES> const& dist,
    std::array<double(*)(double), NUM_REF_STATES> const& ref,
    PIDCtrllers ctrlParams,
    DroneParameters droneParameters, 
    double simTimestep, double finalTime) : 
    m_logger(logger), 
    m_ctrlParams(ctrlParams),
    m_droneParams(droneParameters), 
    m_dist(dist), 
    m_ref(ref), 
    m_simTimestep(simTimestep),
    m_finalTime(finalTime) {}

SimResults DroneTrajectory::Trajectory(SystemState initialState)
{
    SimResults simResults;
    CtrlOut ctrlOut = CascadedPIDController(initialState.plant, initialState.ctrl, 0);
    initialState.ctrl = ctrlOut.ctrlStates;
    initialState.alge = ctrlOut.algeStates;
    simResults.stateProgression.push_back(initialState);
    double time = 0;
    simResults.time.push_back(time);
    while(time <= m_finalTime && simResults.stateProgression.back().stable){
        SystemState prev = simResults.stateProgression.back();
        SystemState state1 = simulateTimestep(prev, time);
        simResults.stateProgression.push_back(state1);
        simResults.stable = simResults.stable & state1.stable;
        time += m_simTimestep;
        simResults.time.push_back(time);

        m_logger << "INFO - timestep: " << time << std::endl;
        m_logger << "INFO - CTRL EFFORT: " << state1.alge(ft) << std::endl;
        m_logger << "INFO - PLANT STATE:" << state1.plant(x) << " " << state1.plant(y) << " " << state1.plant(z) << std::endl;
        m_logger << "INFO - STABLE:" << state1.stable << std::endl;
        
        // SystemState state2 = simulateTimestep(state1, time + m_simTimestep);

        // m_logger << "INFO - timestep: " << time << std::endl;
        // m_logger << "INFO - CTRL EFFORT: " << state1.alge(ft) << std::endl;
        // m_logger << "INFO - PLANT STATE:" << state1.plant(x) << " " << state1.plant(y) << " " << state1.plant(z) << std::endl;
        // m_logger << "INFO - STABLE:" << state1.stable << std::endl;


        // // SystemState stateDouble = simulateTimestep(prev, ref, m_dist, m_ctrlParams, time, m_simTimestep, m_droneParameters);

        // simResults.stateProgression.push_back(state2);
        // simResults.stable = simResults.stable & state2.stable;
        // time += m_simTimestep;
        // simResults.time.push_back(time);

        // double diff = (stateDouble.plant - state2.plant).norm();
        // if( diff < 1e-7 && stateDouble.stable){
        //     m_simTimestep = m_simTimestep * 2;
        // } else if (diff > 1e-4 || !state1.stable || !state2.stable){
        //     if (m_simTimestep > 1e-3) {
        //         m_simTimestep = m_simTimestep / 2;
        //     }
        // }
        // check convergence for stable and end early or unstable and not converging
    }
    return simResults;
}

SystemState DroneTrajectory::simulateTimestep(SystemState prev, double time)
{
    double tol = 1e-12;
    SystemState guess = prev;
    guess.plant += m_simTimestep*f(prev, time-m_simTimestep);
    int count = 0;
    int max_iterations = 100;

    for(; count < max_iterations; count++ ){
        Eigen::Vector<double, NUM_PLANT_STATES> h = H(prev, guess.plant, time);
        if (h.norm() < tol) {
            break;
        }
        Eigen::MatrixX<double> dh = DH(guess);
        if (dh.determinant() != 0) {
            guess.plant = guess.plant - dh.inverse()*h;
        } else {
            guess.stable = false;
            m_logger.warn(std::string("determinant is zero")); 
        }
    }

    CtrlOut ctrlOut = CascadedPIDController(guess.plant, guess.ctrl, time);
    
    guess.ctrl = ctrlOut.ctrlStates;
    guess.alge = ctrlOut.algeStates;
    guess.stable = guess.stable && count < max_iterations;   
    if (count >= max_iterations) {
        m_logger << "WARN: took more than " << max_iterations << " iterations for sim to converge" << std::endl;
    }
    return guess;
}

Eigen::Vector<double, NUM_PLANT_STATES> DroneTrajectory::H(SystemState prev, Eigen::Vector<double, NUM_PLANT_STATES> guess, double time)
{
    Eigen::Vector<double, NUM_PLANT_STATES> fprev = f(prev, time-m_simTimestep);
    Eigen::Vector<double, NUM_PLANT_STATES> fguess = f({guess, prev.ctrl, prev.alge, prev.stable}, time);
    return prev.plant - guess + m_simTimestep/2*(fprev + fguess);
}

Eigen::MatrixX<double> DroneTrajectory::DH(SystemState state)
{
    int n = m_droneParams.numPlantStates;
    Eigen::MatrixX<double> dh = -1*Eigen::MatrixX<double>::Identity(n, n);
    Eigen::MatrixX<double> dfdx_plus = dfdx(state);
    return dh + m_simTimestep/2*dfdx_plus;
}

Eigen::MatrixX<double> DroneTrajectory::dfdx(SystemState state)
{
    int n = m_droneParams.numPlantStates;
    Eigen::MatrixX<double> dfdx = Eigen::MatrixX<double>::Zero(n, n);
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    dfdx(0,6) = 1;
    dfdx(1,7) = 1;
    dfdx(2,8) = 1;

    dfdx(3,3) = plant(q)*cos(plant(phi))*tan(plant(theta)) 
              - plant(r)*sin(plant(phi))*tan(plant(theta));
    dfdx(3,4) = plant(r)*cos(plant(phi))*(pow(tan(plant(theta)), 2) + 1) 
              + plant(q)*sin(plant(phi))*(pow(tan(plant(theta)), 2) + 1);
    dfdx(3,9) = 1;
    dfdx(3,10) = sin(plant(phi))*tan(plant(theta));
    dfdx(3,11) = cos(plant(phi))*tan(plant(theta));
    dfdx(4,3) = - plant(r)*cos(plant(phi)) - plant(q)*sin(plant(phi));
    dfdx(4,10) = cos(plant(phi));
    dfdx(4,11) = -sin(plant(phi));
    dfdx(5,3) = (plant(q)*cos(plant(phi)))/cos(plant(theta)) 
              - (plant(r)*sin(plant(phi)))/cos(plant(theta));
    dfdx(5,4) = (plant(r)*cos(plant(phi))*sin(plant(theta)))/pow(cos(plant(theta)), 2) 
              + (plant(q)*sin(plant(phi))*sin(plant(theta)))/pow(cos(plant(theta)), 2);
    dfdx(5,10) = sin(plant(phi))/cos(plant(theta));
    dfdx(5,11) = cos(plant(phi))/cos(plant(theta));

    dfdx(6,3) = (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass;
    dfdx(6,4) = (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    dfdx(6,5) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    dfdx(7,3) = -(alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    dfdx(7,4) = (alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    dfdx(7,5) = (alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    dfdx(8,3) = -(alge(ft)*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    dfdx(8,4) = -(alge(ft)*cos(plant(phi))*sin(plant(theta)))/m_droneParams.mass;

    dfdx(9,10) = (plant(r)*(m_droneParams.Iy - m_droneParams.Iz))/m_droneParams.Ix;
    dfdx(9,11) = (plant(q)*(m_droneParams.Iy - m_droneParams.Iz))/m_droneParams.Ix;
    dfdx(10,9) = -(plant(r)*(m_droneParams.Ix - m_droneParams.Iz))/m_droneParams.Iy;
    dfdx(10,11) = -(plant(p)*(m_droneParams.Ix - m_droneParams.Iz))/m_droneParams.Iy;
    dfdx(11,9) = (plant(q)*(m_droneParams.Ix - m_droneParams.Iy))/m_droneParams.Iz;
    dfdx(11,10) = (plant(p)*(m_droneParams.Ix - m_droneParams.Iy))/m_droneParams.Iz;

    return dfdx;
}

// dynamics and updated control state
Eigen::Vector<double, NUM_PLANT_STATES> DroneTrajectory::f(SystemState state, double time)
{
    Eigen::Vector<double, NUM_PLANT_STATES> dot = Eigen::Vector<double, NUM_PLANT_STATES>::Zero();
    dot(x) = state.plant(xdot);
    dot(y) = state.plant(ydot);
    dot(z) = state.plant(zdot);

    dot(phi) = state.plant(p) 
                + state.plant(r) * cos(state.plant(phi)) * tan(state.plant(theta))
                + state.plant(q) * sin(state.plant(phi)) * tan(state.plant(theta));
    dot(theta) = state.plant(q) * cos(state.plant(phi)) 
                - state.plant(r) * sin(state.plant(phi));
    dot(psi) = state.plant(r) * cos(state.plant(phi)) / cos(state.plant(theta)) 
                + state.plant(q) * sin(state.plant(phi)) / cos(state.plant(theta));

    dot(xdot) = (state.alge(ft) + m_dist.at(Fwx)(time)) / m_droneParams.mass * 
                (sin(state.plant(phi)) * sin(state.plant(psi)) 
               + cos(state.plant(phi))*cos(state.plant(psi))*sin(state.plant(theta)));
    dot(ydot) = (state.alge(ft)  + m_dist.at(Fwy)(time)) / m_droneParams.mass * 
                (cos(state.plant(phi))*sin(state.plant(psi))*sin(state.plant(theta)) 
               - cos(state.plant(psi))*sin(state.plant(phi)));
    dot(zdot) = - m_droneParams.g + (state.alge(ft)  + m_dist.at(Fwz)(time))/m_droneParams.mass * 
                (cos(state.plant(phi)) * cos(state.plant(theta)));

    dot(p) = (m_droneParams.Iy-m_droneParams.Iz)/m_droneParams.Ix*state.plant(r)*state.plant(q) 
            + (state.alge(tx)  + m_dist.at(Twx)(time))/m_droneParams.Ix;
    dot(q) = (m_droneParams.Iz-m_droneParams.Ix)/m_droneParams.Iy*state.plant(p)*state.plant(r)
            + (state.alge(ty) + m_dist.at(Twy)(time))/m_droneParams.Iy;
    dot(r) = (m_droneParams.Ix-m_droneParams.Iy)/m_droneParams.Iz*state.plant(p)*state.plant(q)
            + (state.alge(tz) + m_dist.at(Twz)(time))/m_droneParams.Iz;
    return dot;
}


// https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L193
CtrlOut DroneTrajectory::CascadedPIDController( 
    Eigen::Vector<double, NUM_PLANT_STATES> plantState,
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> ctrlState, 
    double time)
{
    CtrlOut ctrlOut;
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> newCtrlState;
    newCtrlState(posX) = updatePIDstate(ctrlState(posX), plantState(x), m_ref.at(refx)(time), m_simTimestep);
    newCtrlState(posY) = updatePIDstate(ctrlState(posY), plantState(y), m_ref.at(refy)(time), m_simTimestep);
    newCtrlState(posZ) = updatePIDstate(ctrlState(posZ), plantState(z), m_ref.at(refz)(time), m_simTimestep);
    double desVelX = PIDctrl(m_ctrlParams.posX, newCtrlState(posX));
    double desVelY = PIDctrl(m_ctrlParams.posY, newCtrlState(posY));
    double desVelZ = PIDctrl(m_ctrlParams.posZ, newCtrlState(posZ));

    newCtrlState(velX) = updatePIDstate(ctrlState(velX), plantState(xdot), desVelX, m_simTimestep);
    newCtrlState(velY) = updatePIDstate(ctrlState(velY), plantState(ydot), desVelY, m_simTimestep);
    newCtrlState(velZ) = updatePIDstate(ctrlState(velZ), plantState(zdot), desVelZ, m_simTimestep);

    // output of vel goes into attitude but has negative signs for pitch and roll
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c
    // line 236
    double desAttPhi   = -PIDctrl(m_ctrlParams.velX, newCtrlState(velX));
    double desAttTheta = -PIDctrl(m_ctrlParams.velY, newCtrlState(velY));

    // TODO: Can add constraints later

    // For some reason they scale their thrust by 1000 and add thrust base
    double thrustScale = 1000;
    double thrustBase = std::sqrt(m_droneParams.g*m_droneParams.mass/m_droneParams.kf/4);
    // double thrustBase = 36000;
    double desThrust = PIDctrl(m_ctrlParams.velZ, newCtrlState(velZ))*thrustScale+thrustBase;
    // double desThrust = PIDctrl(m_ctrlParams.velZ, newCtrlState(velZ));

    // Also they have minimum thrust as 20000 ?? 
    // if (desThrust < 20000){
    //     desThrust = 20000;
    // }
    // else if (desThrust > UINT16_MAX) {
    //     desThrust = UINT16_MAX;
    // }
    newCtrlState(attX) = updatePIDstate(ctrlState(attX), plantState(phi), desAttPhi, m_simTimestep);
    newCtrlState(attY) = updatePIDstate(ctrlState(attY), plantState(theta), desAttTheta, m_simTimestep);
    // TODO: Can implement ref yaw, for now fixed
    newCtrlState(attZ) = updatePIDstate(ctrlState(attZ), plantState(psi), 0, m_simTimestep);

    double desAttRateX = PIDctrl(m_ctrlParams.attX, newCtrlState(attX));
    double desAttRateY = PIDctrl(m_ctrlParams.attY, newCtrlState(attY));
    double desAttRateZ = PIDctrl(m_ctrlParams.attZ, newCtrlState(attZ));

    // Assume that [phi dot theta dot psi dot] = [p q r]
    // This assumption holds true for small angles of movement
    newCtrlState(attRateX) = updatePIDstate(ctrlState(attRateX), plantState(p), desAttRateX, m_simTimestep);
    newCtrlState(attRateY) = updatePIDstate(ctrlState(attRateY), plantState(q), desAttRateY, m_simTimestep);
    newCtrlState(attRateZ) = updatePIDstate(ctrlState(attRateZ), plantState(r), desAttRateZ, m_simTimestep);
    ctrlOut.ctrlStates = newCtrlState;
    
    double pitch = PIDctrl(m_ctrlParams.attRateX, newCtrlState(attRateX));
    double roll  = PIDctrl(m_ctrlParams.attRateY, newCtrlState(attRateY));
    double yaw   = PIDctrl(m_ctrlParams.attRateZ, newCtrlState(attRateZ));

    // TODO: implement saturation
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/power_distribution_quadrotor.c
    // line 86 they divide roll and pitch by 2 
    double r = roll / 2;
    double p = pitch / 2;
    double m1 = desThrust - r + p + yaw;
    double m2 = desThrust - r - p - yaw;
    double m3 = desThrust + r - p + yaw;
    double m4 = desThrust + r + p - yaw;

    // https://www.diva-portal.org/smash/get/diva2:860649/FULLTEXT01.pdf
    // b = kf, d = km from https://arxiv.org/pdf/2512.14450
    ctrlOut.algeStates(ft) = m_droneParams.kf * (m1*m1 + m2*m2 + m3*m3 + m4*m4);
    ctrlOut.algeStates(tx) = m_droneParams.kf * m_droneParams.length * (m3*m3 - m1*m1);
    ctrlOut.algeStates(ty) = m_droneParams.kf * m_droneParams.length * (m4*m4 - m2*m2);
    ctrlOut.algeStates(tz) = m_droneParams.km * (m2*m2 + m4*m4 - m1*m1 - m3*m3);

    m_logger << "ctrl ft " <<  ctrlOut.algeStates(ft) << std::endl;

    return ctrlOut;
}

double PIDctrl(PIDParameters params, PIDstate state)
{ return params.ki * state.ki_error + params.kp * state.kp_error + params.kd * state.kd_error;}

PIDstate updatePIDstate(PIDstate currVal, double currSig, double ref, double timestep)
{
    PIDstate newPIDstate;
    newPIDstate.kp_error = ref - currSig;
    newPIDstate.ki_error = currVal.ki_error + newPIDstate.kp_error*timestep;
    newPIDstate.kd_error = 1.0/timestep*(newPIDstate.kp_error - currVal.kp_error);
    return newPIDstate;
}


