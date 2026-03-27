#include "DroneTrajectory.h"

#include <cmath>
#include <Eigen/Dense>

double PIDctrl(PIDParameters params, PIDstate state);
PIDstate updateYawPIDstate(PIDParameters params, PIDstate currVal, double currSig, double ref, double timestep, double time);
PIDstate updatePIDstate(PIDParameters params, PIDstate currVal, double currSig, double ref, double timestep, double time);
std::string printPIDstate(PIDstate state);
Eigen::Vector<PIDstate, NUM_CTRL_STATES> resetPIDstate();
double capAngle(double angle);
double rad2Deg(double angle) { return angle / M_PI * 180.0; }
double deg2Rad(double angle) { return angle / 180.0 * M_PI; }

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
        m_logger << "INFO - PLANT STATE:"   << "x: " << state1.plant(x) << " y: " << state1.plant(y) << " z: " << state1.plant(z) << " phi: "
                                            << state1.plant(phi) << " theta: " << state1.plant(theta) << " psi: " << state1.plant(psi) << " xdot: "
                                            << state1.plant(xdot) << " ydot: " << state1.plant(ydot) << " zdot: " << state1.plant(zdot) << " p: "
                                            << state1.plant(p) << " q: " << state1.plant(q) << " r: " << state1.plant(r) 
                                            << std::endl;
        // TODO: adaptive timestep
        // TOOD: check convergence for stable and end early or unstable and not converging
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
    dfdx(6,4) = (alge(ft)* cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    dfdx(6,5) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;

    dfdx(7,3) = -(alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    dfdx(7,4) = (alge(ft)* cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
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

    dot(ydot) = (state.alge(ft) + m_dist.at(Fwy)(time)) / m_droneParams.mass * 
                (cos(state.plant(phi))*sin(state.plant(psi))*sin(state.plant(theta)) 
               - cos(state.plant(psi))*sin(state.plant(phi)));

    dot(zdot) = - m_droneParams.g + (state.alge(ft) + m_dist.at(Fwz)(time))/m_droneParams.mass * 
                (cos(state.plant(phi)) * cos(state.plant(theta)));

    dot(p) = (m_droneParams.Iy-m_droneParams.Iz)/m_droneParams.Ix*state.plant(r)*state.plant(q) 
            + (state.alge(tx) + m_dist.at(Twx)(time))/m_droneParams.Ix;
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

    //https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L202
    double cosyaw = cos(plantState(psi));
    double sinyaw = sin(plantState(psi));

    // Convert global x y z into body x y z assuming CCW yaw = +ve
    double setp_body_x = m_ref.at(refx)(time) * cosyaw + m_ref.at(refy)(time) * sinyaw;
    double setp_body_y = -m_ref.at(refx)(time) * sinyaw + m_ref.at(refy)(time) * cosyaw;

    double state_body_x = plantState(x) * cosyaw + plantState(y) * sinyaw;
    double state_body_y = -plantState(x) * sinyaw + plantState(y) * cosyaw;    

    newCtrlState(posX) = updatePIDstate(m_ctrlParams.posX, ctrlState(posX), state_body_x, setp_body_x, m_simTimestep, time);
    newCtrlState(posY) = updatePIDstate(m_ctrlParams.posY, ctrlState(posY), state_body_y, setp_body_y, m_simTimestep, time);
    newCtrlState(posZ) = updatePIDstate(m_ctrlParams.posZ, ctrlState(posZ), plantState(z), m_ref.at(refz)(time), m_simTimestep, time);

    double desVelX = PIDctrl(m_ctrlParams.posX, newCtrlState(posX));
    double desVelY = PIDctrl(m_ctrlParams.posY, newCtrlState(posY));
    double desVelZ = PIDctrl(m_ctrlParams.posZ, newCtrlState(posZ));

    double state_body_vx = plantState(xdot) * cosyaw + plantState(ydot) * sinyaw;
    double state_body_vy = -plantState(xdot) * sinyaw + plantState(ydot) * cosyaw;

    newCtrlState(velX) = updatePIDstate(m_ctrlParams.velX, ctrlState(velX), state_body_vx, desVelX, m_simTimestep, time);
    newCtrlState(velY) = updatePIDstate(m_ctrlParams.velY, ctrlState(velY), state_body_vy, desVelY, m_simTimestep, time);
    newCtrlState(velZ) = updatePIDstate(m_ctrlParams.velZ, ctrlState(velZ), plantState(zdot), desVelZ, m_simTimestep, time);
 
    // Need a positive pitch to move forward in the x direction
    // Need a negative roll to move forward in the y direction
    double desPitch   = std::clamp(PIDctrl(m_ctrlParams.velX, newCtrlState(velX)), -m_droneParams.pid_vel_roll_max,  m_droneParams.pid_vel_roll_max);
    double desRoll = -std::clamp(PIDctrl(m_ctrlParams.velY, newCtrlState(velY)), -m_droneParams.pid_vel_pitch_max,  m_droneParams.pid_vel_pitch_max);

    // TODO: Can add constraints later

    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L260
    double thrustScale = 1000;
    double thrustBase = 36000;
    double desThrust = PIDctrl(m_ctrlParams.velZ, newCtrlState(velZ))*thrustScale+thrustBase;
    desThrust = std::clamp(desThrust, 20000.0, (double)UINT16_MAX);

    newCtrlState(roll) = updatePIDstate(m_ctrlParams.roll, ctrlState(roll), rad2Deg(plantState(phi)), desRoll, m_simTimestep, time);
    newCtrlState(pitch) = updatePIDstate(m_ctrlParams.pitch, ctrlState(pitch), rad2Deg(plantState(theta)), desPitch, m_simTimestep, time);
    newCtrlState(yaw) = updateYawPIDstate(m_ctrlParams.yaw, ctrlState(yaw), rad2Deg(plantState(psi)), rad2Deg(m_ref.at(refyaw)(time)), m_simTimestep, time);
        
    double desRollRate = PIDctrl(m_ctrlParams.roll, newCtrlState(roll));
    double desPitchRate = PIDctrl(m_ctrlParams.pitch, newCtrlState(pitch));
    double desYawRate = PIDctrl(m_ctrlParams.yaw, newCtrlState(yaw));
    
    // Assume that [phi dot theta dot psi dot] = [p q r]
    // This assumption holds true for small angles of movement
    newCtrlState(rollRate) = updatePIDstate(m_ctrlParams.rollRate, ctrlState(rollRate), rad2Deg(plantState(p)), desRollRate, m_simTimestep, time);
    newCtrlState(pitchRate) = updatePIDstate(m_ctrlParams.pitchRate, ctrlState(pitchRate), rad2Deg(plantState(q)), desPitchRate, m_simTimestep, time);
    newCtrlState(yawRate) = updatePIDstate(m_ctrlParams.yawRate, ctrlState(yawRate), rad2Deg(plantState(r)), desYawRate, m_simTimestep, time);
    
    ctrlOut.ctrlStates = newCtrlState;

    double rollOutput  = PIDctrl(m_ctrlParams.rollRate, newCtrlState(rollRate));
    double pitchOutput = PIDctrl(m_ctrlParams.pitchRate, newCtrlState(pitchRate));
    double yawOutput = PIDctrl(m_ctrlParams.yawRate, newCtrlState(yawRate));

    // TODO: implement saturation

    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/power_distribution_quadrotor.c
    // line 86 they divide roll and pitch by 2 
    double r = rollOutput / 2;
    double p = pitchOutput / 2;

    double m1 = desThrust - r + p + yawOutput;
    double m2 = desThrust - r - p - yawOutput;
    double m3 = desThrust + r - p + yawOutput;
    double m4 = desThrust + r + p - yawOutput;

    // need to convert motor pwm signal to angular velocity
    // want thrustBase = 36000 PWM == hover --> ft = mg
    // des_wi = sqrt(mg/(kf*4)) 
    double alpha = sqrt(m_droneParams.mass*m_droneParams.g/(m_droneParams.kf*4))/thrustBase;
    
    double w1 = alpha * m1;
    double w2 = alpha * m2;
    double w3 = alpha * m3;
    double w4 = alpha * m4;

    // https://giuseppesilano.net/publications/rosChapter19.pdf
    // the crazyflie firmware using x wing configuration 
    ctrlOut.algeStates(ft) = m_droneParams.kf * (w1*w1 + w2*w2 + w3*w3 + w4*w4);
    ctrlOut.algeStates(tx) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (- w1*w1 - w2*w2 + w3*w3 + w4*w4);
    ctrlOut.algeStates(ty) = m_droneParams.kf * m_droneParams.length * 1/sqrt(2) * (+ w1*w1 - w2*w2 - w3*w3 + w4*w4);
    ctrlOut.algeStates(tz) = m_droneParams.km * (w1*w1 - w2*w2 + w3*w3 - w4*w4);

    return ctrlOut;
}

double PIDctrl(PIDParameters params, PIDstate state)
{ return params.ki * state.ki_error + params.kp * state.kp_error + params.kd * state.kd_error;}

PIDstate updatePIDstate(PIDParameters params, PIDstate currVal, double currSig, double ref, double timestep, double time)
{
    PIDstate newPIDstate;
    newPIDstate.kp_error = ref - currSig;

    double nextIntegral = currVal.ki_error + newPIDstate.kp_error*timestep;
    if (params.integration_limit !=0 && abs(nextIntegral) > params.integration_limit){
        newPIDstate.ki_error = params.integration_limit * nextIntegral/abs(nextIntegral);
    } else {
        newPIDstate.ki_error = currVal.ki_error + newPIDstate.kp_error*timestep; 
    } 
   
    if (time != 0) {
        // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/utils/src/pid.c
        // prevent derivative kick
        double delta = -1.0/timestep*(currSig - currVal.prev_sig);
        if (params.lpf_en){
            newPIDstate.kd_error = currVal.lpf.lpf2pApply(delta);
        } else {
            newPIDstate.kd_error = delta;
        }
        
        
    }
    newPIDstate.prev_sig = currSig;
    newPIDstate.lpf = currVal.lpf;
    return newPIDstate;
}

PIDstate updateYawPIDstate(PIDParameters params, PIDstate currVal, double currSig, double ref, double timestep, double time)
{
    PIDstate newPIDstate;
    newPIDstate.kp_error = capAngle(ref - currSig);
    
    double nextIntegral = currVal.ki_error + newPIDstate.kp_error*timestep;
    if (params.integration_limit !=0 && abs(nextIntegral) > params.integration_limit){
        newPIDstate.ki_error = params.integration_limit * nextIntegral/abs(nextIntegral);
    } else {
        newPIDstate.ki_error = currVal.ki_error + newPIDstate.kp_error*timestep; 
    } 

    if (time != 0) {
        // newPIDstate.kd_error = 1.0/timestep*(newPIDstate.kp_error - currVal.kp_error);
        // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/utils/src/pid.c
        // prevent derivative kick
        newPIDstate.kd_error = -1.0/timestep*(currSig - currVal.prev_sig);
        // newPIDstate.kd_error = 1.0/timestep*(newPIDstate.kp_error - currVal.kp_error);
    }

    newPIDstate.prev_sig = currSig;
    return newPIDstate;
}

std::string printPIDstate(PIDstate state)
{
    return  std::string("kp: ") + std::to_string(state.kp_error) + 
            std::string(" ki: ") + std::to_string(state.ki_error) +
            std::string(" kd: ") + std::to_string(state.kd_error);
}

Eigen::Vector<PIDstate, NUM_CTRL_STATES> resetPIDstate()
{
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> resetState;
    for (int i = 0; i < NUM_CTRL_STATES; i++){
        PIDstate empty;
        resetState(i) = empty;
    }
    return resetState;
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