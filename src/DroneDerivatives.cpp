#include "DroneTrajectory.h"

Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::dfdx(SystemState state)
{
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx = Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>::Zero();
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

Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> DroneTrajectory::dfdz(SystemState state)
{
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz = Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES>::Zero();
    dfdz(7, 1) = 1.0/m_droneParams.mass * 
                (sin(state.plant(phi)) * sin(state.plant(psi)) 
               + cos(state.plant(phi))*cos(state.plant(psi))*sin(state.plant(theta)));

    dfdz(8, 1) = 1.0/m_droneParams.mass * 
                (cos(state.plant(phi))*sin(state.plant(psi))*sin(state.plant(theta)) 
               - cos(state.plant(psi))*sin(state.plant(phi)));

    dfdz(9, 1) = 1.0/m_droneParams.mass * 
                (cos(state.plant(phi)) * cos(state.plant(theta)));

    dfdz(10, 2) = m_droneParams.Ix;
    dfdz(10, 3) = m_droneParams.Iy;
    dfdz(10, 4) = m_droneParams.Iz;

    return dfdz;
}

Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> DroneTrajectory::dgdz(SystemState state)
{
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> dgdz = Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES>::Zero();
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdz(1, epydot) = m_ctrlParams.at(velY).kp;
        dgdz(1, eiydot) = m_ctrlParams.at(velY).ki;
        dgdz(1, edydot) = m_ctrlParams.at(velY).kd;
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdz(1, epxdot) = m_ctrlParams.at(velX).kp;
        dgdz(1, eixdot) = m_ctrlParams.at(velX).ki;
        dgdz(1, edxdot) = m_ctrlParams.at(velX).kd;
    }
    return dgdz;
}

Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_PLANT_STATES> DroneTrajectory::dhdx(SystemState state, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_PLANT_STATES> dhdx = Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_PLANT_STATES>::Zero();
    // TODO: what if ref is not constant?

    dhdx(state_body_x, x) = cos(state.plant(psi));
    dhdx(state_body_x, y) = sin(state.plant(psi));
    dhdx(state_body_x, phi) = -state.plant(x)*sin(state.plant(psi)) + state.plant(y)*cos(state.plant(psi));

    dhdx(state_body_y, x) = -sin(state.plant(psi));
    dhdx(state_body_y, y) = cos(state.plant(psi));
    dhdx(state_body_y, phi) = -state.plant(x)*cos(state.plant(psi)) - state.plant(y)*sin(state.plant(psi));

    dhdx(epz, z) = -1;
    dhdx(edz, z) = 1/timestep;
    dhdx(edz, NUM_PLANT_STATES + z) = -1/timestep;

    dhdx(state_body_vx, xdot) = cos(state.plant(psi));
    dhdx(state_body_vx, ydot) = sin(state.plant(psi));
    dhdx(state_body_vx, phi) = -state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi));

    dhdx(state_body_vy, xdot) = -sin(state.plant(psi));
    dhdx(state_body_vy, ydot) = cos(state.plant(psi));
    dhdx(state_body_vy, phi) = -state.plant(xdot)*cos(state.plant(psi)) - state.plant(ydot)*sin(state.plant(psi));

    dhdx(epzdot, zdot) = -1;
    dhdx(edzdot, zdot) = 1/timestep;
    dhdx(edzdot, NUM_PLANT_STATES + zdot) = -1/timestep;

    dhdx(epphi, phi) = -180/M_PI;
    dhdx(edphi, phi) = 180/(timestep*M_PI);
    dhdx(edphi, NUM_PLANT_STATES + phi) = -1/(timestep*M_PI);

    dhdx(eptheta, theta) = -180/M_PI;
    dhdx(edtheta, theta) = 180/(timestep*M_PI);
    dhdx(edtheta, NUM_PLANT_STATES + theta) = -1/(timestep*M_PI);

    dhdx(eppsi, psi) = -180/M_PI;
    dhdx(edpsi, psi) = 180/(timestep*M_PI);
    dhdx(edpsi, NUM_PLANT_STATES + psi) = -1/(timestep*M_PI);

    dhdx(delay_1_rollRate, p) = 180/(timestep*M_PI);
    dhdx(delay_1_rollRate, NUM_PLANT_STATES+p) = -180/(timestep*M_PI);

    dhdx(delay_1_pitchRate, q) = 180/(timestep*M_PI);
    dhdx(delay_1_pitchRate, NUM_PLANT_STATES+q) = -180/(timestep*M_PI);

    dhdx(epp, p) = -180/M_PI;
    dhdx(edp, p) = lpf_b0*180/(timestep*M_PI);
    dhdx(edp, NUM_PLANT_STATES + p) = -lpf_b0*1/(timestep*M_PI);

    dhdx(epq, q) = -180/M_PI;
    dhdx(epq, q) = lpf_b0*180/(timestep*M_PI);
    dhdx(epq, NUM_PLANT_STATES + q) = -lpf_b0*1/(timestep*M_PI);

    dhdx(epr, r) = -180/M_PI;
    dhdx(edr, r) = 180/(timestep*M_PI);
    dhdx(edr, NUM_PLANT_STATES + r) = -1/(timestep*M_PI);

    return dhdx;
}

Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_Z_STATES> DroneTrajectory::dhdz(SystemState state, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_Z_STATES> dhdz = Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_Z_STATES>::Zero();
    dhdz(epx, setp_body_x) = 1;
    dhdz(epx, state_body_x) = -1;
    dhdz(eix, epx) = timestep;
    dhdz(eix, NUM_Z_STATES+eix) = 1;
    dhdz(desVelX, epx) = m_ctrlParams.at(posX).kp;
    dhdz(desVelX, eix) = m_ctrlParams.at(posX).ki;
    dhdz(desVelX, edx) = m_ctrlParams.at(posX).kd;

    dhdz(epy, setp_body_y) = 1;
    dhdz(epy, state_body_y) = -1;
    dhdz(eiy, epy) = timestep;
    dhdz(eiy, NUM_Z_STATES+eiy) = 1;
    dhdz(desVelY, epy) = m_ctrlParams.at(posY).kp;
    dhdz(desVelY, eiy) = m_ctrlParams.at(posY).ki;
    dhdz(desVelY, edy) = m_ctrlParams.at(posY).kd;

    dhdz(eiz, epz) = timestep;
    dhdz(eiz, NUM_Z_STATES+eiz) = 1;
    dhdz(desVelZ, epz) = m_ctrlParams.at(posZ).kp;
    dhdz(desVelZ, eiz) = m_ctrlParams.at(posZ).ki;
    dhdz(desVelZ, edz) = m_ctrlParams.at(posZ).kd;

    dhdz(epxdot, desVelX) = 1;
    dhdz(epxdot, state_body_vx) = -1;
    dhdz(eixdot, epxdot) = timestep;
    dhdz(eixdot, NUM_Z_STATES+eixdot) = 1;
    dhdz(edxdot, state_body_vx) = 1/timestep;
    dhdz(edxdot, NUM_Z_STATES+state_body_vx) = -1/timestep;

    dhdz(epydot, desVelY) = 1;
    dhdz(epydot, state_body_vy) = -1;
    dhdz(eiydot, epydot) = timestep;
    dhdz(eiydot, NUM_Z_STATES+eiydot) = 1;
    dhdz(edydot, state_body_vy) = 1/timestep;
    dhdz(edydot, NUM_Z_STATES+state_body_vy) = -1/timestep;

    dhdz(epzdot, desVelZ) = 1;
    dhdz(eizdot, epzdot) = timestep;
    dhdz(eizdot, NUM_Z_STATES+eizdot) = 1;
    dhdz(desThrust, epzdot) = m_ctrlParams.at(velZ).kp*m_thrustScale;
    dhdz(desThrust, eizdot) = m_ctrlParams.at(velZ).ki*m_thrustScale;
    dhdz(desThrust, edzdot) = m_ctrlParams.at(velZ).kd*m_thrustScale;

    dhdz(eiphi, epphi) = timestep;
    dhdz(eiphi, NUM_Z_STATES+eiphi) = 1;
    dhdz(desRollRate, epphi) = m_ctrlParams.at(roll).kp;
    dhdz(desRollRate, eiphi) = m_ctrlParams.at(roll).ki;
    dhdz(desRollRate, edphi) = m_ctrlParams.at(roll).kd;

    dhdz(eitheta, eptheta) = timestep;
    dhdz(eitheta, NUM_Z_STATES+eitheta) = 1;
    dhdz(desPitchRate, eptheta) = m_ctrlParams.at(pitch).kp;
    dhdz(desPitchRate, eitheta) = m_ctrlParams.at(pitch).ki;
    dhdz(desPitchRate, edtheta) = m_ctrlParams.at(pitch).kd;

    dhdz(eipsi, eppsi) = timestep;
    dhdz(eipsi, NUM_Z_STATES+eipsi) = 1;
    dhdz(desYawRate, eppsi) = m_ctrlParams.at(yaw).kp;
    dhdz(desYawRate, eipsi) = m_ctrlParams.at(yaw).ki;
    dhdz(desYawRate, edpsi) = m_ctrlParams.at(yaw).kd;

    dhdz(delay_1_rollRate, NUM_Z_STATES+delay_1_rollRate) = -lpf_a1;
    dhdz(delay_1_rollRate, NUM_Z_STATES+delay_2_rollRate) = -lpf_a2;
    dhdz(delay_2_rollRate, NUM_Z_STATES+delay_1_rollRate) = -1;
    
    dhdz(delay_1_pitchRate, NUM_Z_STATES+delay_1_pitchRate) = -lpf_a1;
    dhdz(delay_1_pitchRate, NUM_Z_STATES+delay_2_pitchRate) = -lpf_a2;
    dhdz(delay_2_pitchRate, NUM_Z_STATES+delay_1_pitchRate) = -1;

    dhdz(epp, desRollRate) = 1;
    dhdz(eip, epp) = timestep;
    dhdz(eip, NUM_Z_STATES+eip) = 1;
    dhdz(edp, NUM_Z_STATES+delay_1_rollRate) = -lpf_b0*lpf_a1+lpf_b1;
    dhdz(edp, NUM_Z_STATES+delay_2_rollRate) = -lpf_b0*lpf_a2+lpf_b2;
    dhdz(desRollOutput, epp) = m_ctrlParams.at(rollRate).kp;
    dhdz(desRollOutput, eip) = m_ctrlParams.at(rollRate).ki;
    dhdz(desRollOutput, edp) = m_ctrlParams.at(rollRate).kd;

    dhdz(epq, desPitchRate) = 1;
    dhdz(eiq, epq) = timestep;
    dhdz(eiq, NUM_Z_STATES+eiq) = 1;
    dhdz(edq, NUM_Z_STATES+delay_1_pitchRate) = -lpf_b0*lpf_a1+lpf_b1;
    dhdz(edq, NUM_Z_STATES+delay_2_pitchRate) = -lpf_b0*lpf_a2+lpf_b2;
    dhdz(desPitchOutput, epq) = m_ctrlParams.at(pitchRate).kp;
    dhdz(desPitchOutput, eiq) = m_ctrlParams.at(pitchRate).ki;
    dhdz(desPitchOutput, edq) = m_ctrlParams.at(pitchRate).kd;

    dhdz(epr, desYawRate) = 1;
    dhdz(eir, epr) = timestep;
    dhdz(eir, NUM_Z_STATES+eir) = 1;
    dhdz(desYawOutput, epr) = m_ctrlParams.at(yawRate).kp;
    dhdz(desYawOutput, eir) = m_ctrlParams.at(yawRate).ki;
    dhdz(desYawOutput, edr) = m_ctrlParams.at(yawRate).kd;

    dhdz(w1, desThrust)      = m_alpha;
    dhdz(w1, desRollOutput)  = -m_alpha/2;
    dhdz(w1, desPitchOutput) = m_alpha/2;
    dhdz(w1, desYawOutput)   = m_alpha;

    dhdz(w2, desThrust)      = m_alpha;
    dhdz(w2, desRollOutput)  = -m_alpha/2;
    dhdz(w2, desPitchOutput) = -m_alpha/2;
    dhdz(w2, desYawOutput)   = -m_alpha;

    dhdz(w3, desThrust)      = m_alpha;
    dhdz(w3, desRollOutput)  = m_alpha/2;
    dhdz(w3, desPitchOutput) = -m_alpha/2;
    dhdz(w3, desYawOutput)   = m_alpha;

    dhdz(w4, desThrust)      = m_alpha;
    dhdz(w4, desRollOutput)  = m_alpha/2;
    dhdz(w4, desPitchOutput) = m_alpha/2;
    dhdz(w4, desYawOutput)   = -m_alpha;
    
    dhdz(ft, w1) = m_droneParams.kf * 2 *state.alge(w1);
    dhdz(ft, w2) = m_droneParams.kf * 2 *state.alge(w2);
    dhdz(ft, w3) = m_droneParams.kf * 2 *state.alge(w3);
    dhdz(ft, w4) = m_droneParams.kf * 2 *state.alge(w4);

    dhdz(tx, w1) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w1);
    dhdz(tx, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w2);
    dhdz(tx, w3) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w3);
    dhdz(tx, w4) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w4);

    dhdz(ty, w1) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w1);
    dhdz(ty, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w2);
    dhdz(ty, w3) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w3);
    dhdz(ty, w4) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w4);

    dhdz(tz, w1) = m_droneParams.km * state.alge(w1);
    dhdz(tz, w2) = - m_droneParams.km * state.alge(w2);
    dhdz(tz, w3) = m_droneParams.km * state.alge(w3);
    dhdz(tz, w4) = - m_droneParams.km * state.alge(w4);

    return dhdz;
}

Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> DroneTrajectory::dhdy()
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> dhdy = Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES>::Zero();
    dhdy(epphi, 1) = 1;
    dhdy(eptheta, 2) = 1;
    return dhdy;
}

