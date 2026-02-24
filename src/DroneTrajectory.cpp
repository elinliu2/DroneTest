#include "DroneTrajectory.h"
#include <cmath>
#include <Eigen/Dense>

DroneTrajectory::DroneTrajectory(
    PIDCtrllers ctrlParams,
    DroneParameters droneParameters, 
    std::array<double(*)(double), NUM_DIST_STATES> const& dist,
    double simTimestep, double finalTime) : 
    m_ctrlParams(ctrlParams),
    m_droneParameters(droneParameters), 
    m_dist(dist), 
    m_simTimestep(simTimestep),
    m_finalTime(finalTime) {}

SimResults DroneTrajectory::Trajectory(std::array<double(*)(double), NUM_REF_STATES> const& ref, SystemState const & initialState)
{
    SimResults simResults;
    simResults.stateProgression.push_back(initialState);
    double time = 0;
    simResults.time.push_back(time);
    while(time <= m_finalTime){
        SystemState prev = simResults.stateProgression.back();
        // Simulate two timesteps
        // Simulate double timestep
        // Possibly adjust timestep
        // Update simResults
        time += m_simTimestep;
        simResults.time.push_back(time);
        // check convergence
    }
}

SystemState simulateTimestep(Eigen::Vector<double, NUM_PLANT_STATES> prev, double timestep, 
    std::array<double(*)(double), NUM_REF_STATES> const& ref, 
    std::array<double(*)(double), NUM_DIST_STATES> const& dist, 
    PIDCtrllers ctrlParams, 
    DroneParameters droneParams)
{
    double tol = 1e-12;
    SystemState guess;
    
}

SystemState H(SystemState prev, 
        Eigen::Vector<double, NUM_PLANT_STATES> guess, 
        std::array<double(*)(double), NUM_REF_STATES> const& ref, 
        std::array<double(*)(double), NUM_DIST_STATES> const& dist, 
        PIDCtrllers ctrlParams,
        double time, 
        double timestep,
        DroneParameters droneParams)
{
    SystemState fprev = f(prev, ref, dist, ctrlParams, time, timestep, droneParams);
    SystemState fguess = f({guess, prev.ctrl, true}, ref, dist, ctrlParams, time, timestep, droneParams);

    return {prev.plant - guess + timestep/2*(fprev.plant - fguess.plant), fguess.ctrl, true};
}

Eigen::MatrixX<double> DH(SystemState state, DroneParameters const& droneParams, double timestep)
{
    int n = droneParams.numPlantStates;
    Eigen::MatrixX<double> dh(n*n);
    for (int i = 0; i < n; i++){
        dh(i,i) = -1;
    }
    Eigen::MatrixX<double> dfdx = dfdx(state, droneParams);
    dh = dh + timestep/2*dfdx;
    return dh;
}

Eigen::MatrixX<double> dfdx(SystemState state, DroneParameters const& droneParams)
{
    int n = droneParams.numPlantStates;
    Eigen::MatrixX<double> dfdx(n*n);
    
    return dfdx;
}

// dynamics and updated control state
SystemState f(SystemState state, 
        std::array<double(*)(double), NUM_REF_STATES> const& ref, 
        std::array<double(*)(double), NUM_DIST_STATES> const& dist, 
        PIDCtrllers ctrlParams,
        double time, 
        double timestep,
        DroneParameters const& droneParams
        )
{
    SystemState newState;
    Eigen::Vector<double, NUM_PLANT_STATES> dot;
    dot(plantIndex::x) = state.plant(plantIndex::xdot);
    dot(plantIndex::y) = state.plant(plantIndex::ydot);
    dot(plantIndex::z) = state.plant(plantIndex::zdot);

    dot(plantIndex::phi) = state.plant(plantIndex::p) 
                + state.plant(plantIndex::r) * std::cos(state.plant(plantIndex::phi)) * std::tan(state.plant(plantIndex::theta))
                + state.plant(plantIndex::q) * std::sin(state.plant(plantIndex::phi)) * std::tan(state.plant(plantIndex::theta));
    dot(plantIndex::theta) = state.plant(plantIndex::q) * std::cos(state.plant(plantIndex::phi)) 
                - state.plant(plantIndex::r) * std::sin(state.plant(plantIndex::phi));
    dot(plantIndex::psi) = state.plant(plantIndex::r) * std::cos(state.plant(plantIndex::phi)) / std::cos(state.plant(plantIndex::theta)) 
                + state.plant(plantIndex::q) * std::sin(state.plant(plantIndex::phi)) / std::sin(state.plant(plantIndex::theta));

    CtrlOut ctrlOut = CascadedPIDController(ref, dist, ctrlParams, time, timestep, state.plant, state.ctrl, droneParams);

    dot(plantIndex::xdot) = -(ctrlOut.ft + dist.at(distIndex::Fwx)(time)) / droneParams.mass * 
                (std::sin(state.plant(plantIndex::phi)) * std::sin(state.plant(plantIndex::psi)) 
               + std::cos(state.plant(plantIndex::phi))*std::cos(state.plant(plantIndex::psi))*std::sin(state.plant(plantIndex::theta)));
    dot(plantIndex::ydot) = -(ctrlOut.ft + dist.at(distIndex::Fwx)(time)) / droneParams.mass * 
                (std::cos(state.plant(plantIndex::phi))*std::sin(state.plant(plantIndex::psi))*std::sin(state.plant(plantIndex::theta)) 
               - std::cos(state.plant(plantIndex::psi))*std::sin(state.plant(plantIndex::phi)));
    dot(plantIndex::zdot) = droneParams.g -(ctrlOut.ft + dist.at(distIndex::Fwx)(time))/droneParams.mass * 
                (std::cos(state.plant(plantIndex::phi)) * std::cos(state.plant(plantIndex::theta)));
    dot(plantIndex::p) = (droneParams.Iy-droneParams.Iz)/droneParams.Ix*state.plant(plantIndex::r)*state.plant(plantIndex::q) 
            + (ctrlOut.tx + dist.at(distIndex::Twx)(time))/droneParams.Ix;
    dot(plantIndex::q) = (droneParams.Iz-droneParams.Ix)/droneParams.Iy*state.plant(plantIndex::p)*state.plant(plantIndex::r)
            + (ctrlOut.ty + dist.at(distIndex::Twy)(time))/droneParams.Iy;
    dot(plantIndex::r) = (droneParams.Ix-droneParams.Iy)/droneParams.Iz*state.plant(plantIndex::p)*state.plant(plantIndex::q)
            + (ctrlOut.tz + dist.at(distIndex::Twz)(time))/droneParams.Iz;
    newState.plant = dot;
    newState.ctrl = ctrlOut.ctrlState;
    return newState;
}

CtrlOut CascadedPIDController( 
    std::array<double(*)(double), NUM_REF_STATES> const& ref, 
    std::array<double(*)(double), NUM_DIST_STATES> const& dist,
    PIDCtrllers ctrlParams,
    double time,
    double timestep,
    Eigen::Vector<double, NUM_PLANT_STATES> plantState,
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> ctrlState, 
    DroneParameters const& droneParams)
{
    CtrlOut ctrlOut;
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> newCtrlState;
    newCtrlState(ctrlIndex::posX) = updatePIDstate(ctrlState(ctrlIndex::posX), plantState(plantIndex::x), ref.at(refIndex::x)(time), timestep);
    newCtrlState(ctrlIndex::posY) = updatePIDstate(ctrlState(ctrlIndex::posY), plantState(plantIndex::y), ref.at(refIndex::y)(time), timestep);
    newCtrlState(ctrlIndex::posZ) = updatePIDstate(ctrlState(ctrlIndex::posZ), plantState(plantIndex::z), ref.at(refIndex::z)(time), timestep);

    double desVelX = PIDctrl(ctrlParams.posX, newCtrlState(ctrlIndex::posX));
    double desVelY = PIDctrl(ctrlParams.posY, newCtrlState(ctrlIndex::posY));
    double desVelZ = PIDctrl(ctrlParams.posZ, newCtrlState(ctrlIndex::posZ));

    newCtrlState(ctrlIndex::velX) = updatePIDstate(ctrlState(ctrlIndex::velX), plantState(plantIndex::xdot), desVelX, timestep);
    newCtrlState(ctrlIndex::velY) = updatePIDstate(ctrlState(ctrlIndex::velY), plantState(plantIndex::ydot), desVelY, timestep);
    newCtrlState(ctrlIndex::velZ) = updatePIDstate(ctrlState(ctrlIndex::velZ), plantState(plantIndex::zdot), desVelZ, timestep);

    // output of vel goes into attitude but has negative signs for pitch and roll
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c
    // line 236
    double desAttPhi   = -PIDctrl(ctrlParams.velX, newCtrlState(ctrlIndex::velX));
    double desAttTheta = -PIDctrl(ctrlParams.velY, newCtrlState(ctrlIndex::velY));

    // TODO: Can add constraints later

    // For some reason they scale their thrust by 1000 and add thrust base
    double thrustScale = 1000;
    double thrustBase = droneParams.g*droneParams.mass/4;
    double desThrust = PIDctrl(ctrlParams.velZ, newCtrlState(ctrlIndex::velZ))*thrustScale+thrustBase;

    // Also they have minimum thrust as 20000 ?? 
    // TODO: Can add constraints later 
    if (desThrust < 2000)
    {
        desThrust = 2000;
    }
    newCtrlState(ctrlIndex::attX) = updatePIDstate(ctrlState(ctrlIndex::attX), plantState(plantIndex::phi), desAttPhi, timestep);
    newCtrlState(ctrlIndex::attY) = updatePIDstate(ctrlState(ctrlIndex::attY), plantState(plantIndex::theta), desAttTheta, timestep);
    // TODO: Can implement ref yaw, for now fixed
    newCtrlState(ctrlIndex::attZ) = updatePIDstate(ctrlState(ctrlIndex::attZ), plantState(plantIndex::psi), 0, timestep);

    double desAttRateX = PIDctrl(ctrlParams.attX, newCtrlState(ctrlIndex::attX));
    double desAttRateY = PIDctrl(ctrlParams.attY, newCtrlState(ctrlIndex::attY));
    double desAttRateZ = PIDctrl(ctrlParams.attZ, newCtrlState(ctrlIndex::attZ));

    // Assume that [phi dot theta dot psi dot] = [p q r]
    // This assumption holds true for small angles of movement
    newCtrlState(ctrlIndex::attRateX) = updatePIDstate(ctrlState(ctrlIndex::attRateX), plantState(plantIndex::p), desAttRateX, timestep);
    newCtrlState(ctrlIndex::attRateY) = updatePIDstate(ctrlState(ctrlIndex::attRateY), plantState(plantIndex::q), desAttRateY, timestep);
    newCtrlState(ctrlIndex::attRateZ) = updatePIDstate(ctrlState(ctrlIndex::attRateZ), plantState(plantIndex::r), desAttRateZ, timestep);
    ctrlOut.ctrlState = newCtrlState;
    
    double roll =   PIDctrl(ctrlParams.attRateX, newCtrlState(ctrlIndex::attRateX));
    double pitch =  PIDctrl(ctrlParams.attRateY, newCtrlState(ctrlIndex::attRateY));
    double yaw =    PIDctrl(ctrlParams.attRateZ, newCtrlState(ctrlIndex::attRateZ));

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
    ctrlOut.ft = droneParams.kf * (m1*m1 + m2*m2 + m3*m3 + m4*m4);
    ctrlOut.tx = droneParams.kf * droneParams.length * (m3*m3 - m1*m1);
    ctrlOut.ty = droneParams.kf * droneParams.length * (m4*m4 - m2*m2);
    ctrlOut.tz = droneParams.km * (m2*m2 + m4*m4 - m1*m1 - m3*m3);

    return ctrlOut;
}

double PIDctrl(PIDParameters params, PIDstate state)
{ return params.ki * state.ki_error + params.kp * state.kp_error + params.kd * state.kd_error;}

PIDstate updatePIDstate(PIDstate currVal, double currSig, double ref, double timestep)
{
    PIDstate newPIDstate;
    newPIDstate.kp_error = ref - currSig;
    newPIDstate.ki_error = currVal.ki_error + newPIDstate.kp_error;
    newPIDstate.kd_error = 1.0/timestep*(currVal.kp_error - newPIDstate.kp_error);
    return newPIDstate;
}


