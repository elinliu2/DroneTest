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

Eigen::Vector<double, NUM_PLANT_STATES> H(SystemState prev, 
        Eigen::Vector<double, NUM_PLANT_STATES> guess, 
        std::array<double(*)(double), NUM_REF_STATES> const& ref, 
        std::array<double(*)(double), NUM_DIST_STATES> const& dist, 
        PIDCtrllers ctrlParams,
        double time, 
        double timestep,
        DroneParameters droneParams)
{
    Eigen::Vector<double, NUM_PLANT_STATES> fprev = f(prev, ref, dist, ctrlParams, time, timestep, droneParams);
    Eigen::Vector<double, NUM_PLANT_STATES> fguess = f({guess, prev.ctrl, prev.alge, true}, ref, dist, ctrlParams, time, timestep, droneParams);

    return prev.plant - guess + timestep/2*(fprev - fguess);
}

Eigen::MatrixX<double> DH(SystemState state, DroneParameters const& droneParams, double timestep)
{
    int n = droneParams.numPlantStates;
    Eigen::MatrixX<double> dh = -1*Eigen::MatrixX<double>::Identity(n, n);
    Eigen::MatrixX<double> dfdx_plus = dfdx(state, droneParams);
    return dh + timestep/2*dfdx_plus;
}

Eigen::MatrixX<double> dfdx(SystemState state, DroneParameters const& droneParams)
{
    int n = droneParams.numPlantStates;
    Eigen::MatrixX<double> dfdx(n, n);
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_PLANT_STATES> alge = state.alge;
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

    dfdx(6,3) = -(alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/droneParams.mass;
    dfdx(6,4) = -(alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/droneParams.mass;
    dfdx(6,5) = -(alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/droneParams.mass;
    dfdx(7,3) = (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/droneParams.mass;
    dfdx(7,4) = -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/droneParams.mass;
    dfdx(7,5) = -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/droneParams.mass;
    dfdx(8,3) = (alge(ft)*cos(plant(theta))*sin(plant(phi)))/droneParams.mass;
    dfdx(8,4) = (alge(ft)*cos(plant(phi))*sin(plant(theta)))/droneParams.mass;

    dfdx(9,10) = (plant(r)*(droneParams.Iy - droneParams.Iz))/droneParams.Ix;
    dfdx(9,11) = (plant(q)*(droneParams.Iy - droneParams.Iz))/droneParams.Ix;
    dfdx(10,9) = -(plant(r)*(droneParams.Ix - droneParams.Iz))/droneParams.Iy;
    dfdx(10,11) = -(plant(p)*(droneParams.Ix - droneParams.Iz))/droneParams.Iy;
    dfdx(11,9) = (plant(q)*(droneParams.Ix - droneParams.Iy))/droneParams.Iz;
    dfdx(11,10) = (plant(p)*(droneParams.Ix - droneParams.Iy))/droneParams.Iz;

    return dfdx;
}

// dynamics and updated control state
 Eigen::Vector<double, NUM_PLANT_STATES> f(SystemState state, 
        std::array<double(*)(double), NUM_REF_STATES> const& ref, 
        std::array<double(*)(double), NUM_DIST_STATES> const& dist, 
        PIDCtrllers ctrlParams,
        double time, 
        double timestep,
        DroneParameters const& droneParams
        )
{
    Eigen::Vector<double, NUM_PLANT_STATES> dot;
    CtrlOut ctrlOut = CascadedPIDController(ref, dist, ctrlParams, time, timestep, state.plant, state.ctrl, droneParams);

    dot(x) = state.plant(xdot);
    dot(y) = state.plant(ydot);
    dot(z) = state.plant(zdot);

    dot(phi) = state.plant(p) 
                + state.plant(r) * cos(state.plant(phi)) * tan(state.plant(theta))
                + state.plant(q) * sin(state.plant(phi)) * tan(state.plant(theta));
    dot(theta) = state.plant(q) * cos(state.plant(phi)) 
                - state.plant(r) * sin(state.plant(phi));
    dot(psi) = state.plant(r) * cos(state.plant(phi)) / cos(state.plant(theta)) 
                + state.plant(q) * sin(state.plant(phi)) / sin(state.plant(theta));

    dot(xdot) = -(ctrlOut.algeStates(ft) + dist.at(Fwx)(time)) / droneParams.mass * 
                (sin(state.plant(phi)) * sin(state.plant(psi)) 
               + cos(state.plant(phi))*cos(state.plant(psi))*sin(state.plant(theta)));
    dot(ydot) = -(ctrlOut.algeStates(ft) + dist.at(Fwx)(time)) / droneParams.mass * 
                (cos(state.plant(phi))*sin(state.plant(psi))*sin(state.plant(theta)) 
               - cos(state.plant(psi))*sin(state.plant(phi)));
    dot(zdot) = droneParams.g -(ctrlOut.algeStates(ft) + dist.at(Fwx)(time))/droneParams.mass * 
                (cos(state.plant(phi)) * cos(state.plant(theta)));
    
                dot(p) = (droneParams.Iy-droneParams.Iz)/droneParams.Ix*state.plant(r)*state.plant(q) 
            + (ctrlOut.algeStates(tx) + dist.at(Twx)(time))/droneParams.Ix;
    dot(q) = (droneParams.Iz-droneParams.Ix)/droneParams.Iy*state.plant(p)*state.plant(r)
            + (ctrlOut.algeStates(ty) + dist.at(Twy)(time))/droneParams.Iy;
    dot(r) = (droneParams.Ix-droneParams.Iy)/droneParams.Iz*state.plant(p)*state.plant(q)
            + (ctrlOut.algeStates(tz) + dist.at(Twz)(time))/droneParams.Iz;
    return dot;
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
    newCtrlState(::posX) = updatePIDstate(ctrlState(::posX), plantState(x), ref.at(refx)(time), timestep);
    newCtrlState(::posY) = updatePIDstate(ctrlState(::posY), plantState(y), ref.at(refy)(time), timestep);
    newCtrlState(::posZ) = updatePIDstate(ctrlState(::posZ), plantState(z), ref.at(refz)(time), timestep);

    double desVelX = PIDctrl(ctrlParams.posX, newCtrlState(::posX));
    double desVelY = PIDctrl(ctrlParams.posY, newCtrlState(::posY));
    double desVelZ = PIDctrl(ctrlParams.posZ, newCtrlState(::posZ));

    newCtrlState(::velX) = updatePIDstate(ctrlState(::velX), plantState(xdot), desVelX, timestep);
    newCtrlState(::velY) = updatePIDstate(ctrlState(::velY), plantState(ydot), desVelY, timestep);
    newCtrlState(::velZ) = updatePIDstate(ctrlState(::velZ), plantState(zdot), desVelZ, timestep);

    // output of vel goes into attitude but has negative signs for pitch and roll
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c
    // line 236
    double desAttPhi   = -PIDctrl(ctrlParams.velX, newCtrlState(::velX));
    double desAttTheta = -PIDctrl(ctrlParams.velY, newCtrlState(::velY));

    // TODO: Can add constraints later

    // For some reason they scale their thrust by 1000 and add thrust base
    double thrustScale = 1000;
    double thrustBase = droneParams.g*droneParams.mass/4;
    double desThrust = PIDctrl(ctrlParams.velZ, newCtrlState(::velZ))*thrustScale+thrustBase;

    // Also they have minimum thrust as 20000 ?? 
    // TODO: Can add constraints later 
    if (desThrust < 2000)
    {
        desThrust = 2000;
    }
    newCtrlState(::attX) = updatePIDstate(ctrlState(::attX), plantState(phi), desAttPhi, timestep);
    newCtrlState(::attY) = updatePIDstate(ctrlState(::attY), plantState(theta), desAttTheta, timestep);
    // TODO: Can implement ref yaw, for now fixed
    newCtrlState(::attZ) = updatePIDstate(ctrlState(::attZ), plantState(psi), 0, timestep);

    double desAttRateX = PIDctrl(ctrlParams.attX, newCtrlState(::attX));
    double desAttRateY = PIDctrl(ctrlParams.attY, newCtrlState(::attY));
    double desAttRateZ = PIDctrl(ctrlParams.attZ, newCtrlState(::attZ));

    // Assume that [phi dot theta dot psi dot] = [p q r]
    // This assumption holds true for small angles of movement
    newCtrlState(::attRateX) = updatePIDstate(ctrlState(::attRateX), plantState(p), desAttRateX, timestep);
    newCtrlState(::attRateY) = updatePIDstate(ctrlState(::attRateY), plantState(q), desAttRateY, timestep);
    newCtrlState(::attRateZ) = updatePIDstate(ctrlState(::attRateZ), plantState(r), desAttRateZ, timestep);
    ctrlOut.ctrlStates = newCtrlState;
    
    double pitch = PIDctrl(ctrlParams.attRateX, newCtrlState(::attRateX));
    double roll  = PIDctrl(ctrlParams.attRateY, newCtrlState(::attRateY));
    double yaw   = PIDctrl(ctrlParams.attRateZ, newCtrlState(::attRateZ));

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
    ctrlOut.algeStates(ft) = droneParams.kf * (m1*m1 + m2*m2 + m3*m3 + m4*m4);
    ctrlOut.algeStates(tx) = droneParams.kf * droneParams.length * (m3*m3 - m1*m1);
    ctrlOut.algeStates(ty) = droneParams.kf * droneParams.length * (m4*m4 - m2*m2);
    ctrlOut.algeStates(tz) = droneParams.km * (m2*m2 + m4*m4 - m1*m1 - m3*m3);

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

double calcFt()
{

}


