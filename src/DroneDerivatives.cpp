#include "DroneTrajectory.h"
typedef Eigen::Triplet<double> T;

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

Eigen::SparseMatrix<double> DroneTrajectory::dfdz(SystemState state)
{
    std::vector<T> dfdz;
    dfdz.reserve(6);

    dfdz.push_back(T(7, 1, 1.0/m_droneParams.mass * 
                (std::sin(state.plant(phi)) * std::sin(state.plant(psi)) 
               + std::cos(state.plant(phi))*std::cos(state.plant(psi))*std::sin(state.plant(theta)))));

    dfdz.push_back(T(8, 1, 1.0/m_droneParams.mass * 
                (std::cos(state.plant(phi))*std::sin(state.plant(psi))*std::sin(state.plant(theta)) 
               - std::cos(state.plant(psi))*std::sin(state.plant(phi)))));

    dfdz.push_back(T(9, 1, 1.0/m_droneParams.mass * 
                (std::cos(state.plant(phi)) * std::cos(state.plant(theta)))));

    dfdz.push_back(T(10, 2, m_droneParams.Ix));
    dfdz.push_back(T(10, 3, m_droneParams.Iy));
    dfdz.push_back(T(10, 4, m_droneParams.Iz));

    Eigen::SparseMatrix<double> dfdz_mat(NUM_PLANT_STATES, NUM_Z_STATES);
    dfdz_mat.setFromTriplets(dfdz.begin(), dfdz.end());
    return dfdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dgdz(SystemState state)
{
    std::vector<T> dgdz;
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdz.reserve(3);
        dgdz.push_back(T(1, epydot, m_ctrlParams.at(velY).kp));
        dgdz.push_back(T(1, eiydot, m_ctrlParams.at(velY).ki));
        dgdz.push_back(T(1, edydot, m_ctrlParams.at(velY).kd));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdz.reserve(3);
        dgdz.push_back(T(1, epxdot, m_ctrlParams.at(velX).kp));
        dgdz.push_back(T(1, eixdot, m_ctrlParams.at(velX).ki));
        dgdz.push_back(T(1, edxdot, m_ctrlParams.at(velX).kd));
    }
    Eigen::SparseMatrix<double> dgdz_mat(NUM_Y_STATES, NUM_Z_STATES);
    dgdz_mat.setFromTriplets(dgdz.begin(), dgdz.end());
    return dgdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxPlus(SystemState state, double timestep)
{
    std::vector<T> dhdx;
    dhdx.reserve(30);

    // TODO: what if ref is not constant?
    dhdx.push_back(T(state_body_x, x, std::cos(state.plant(psi))));
    dhdx.push_back(T(state_body_x, y, std::sin(state.plant(psi))));
    dhdx.push_back(T(state_body_x, phi, -state.plant(x)*std::sin(state.plant(psi)) + state.plant(y)*std::cos(state.plant(psi))));

    dhdx.push_back(T(state_body_y, x, -std::sin(state.plant(psi))));
    dhdx.push_back(T(state_body_y, y, std::cos(state.plant(psi))));
    dhdx.push_back(T(state_body_y, phi, -state.plant(x)*std::cos(state.plant(psi)) - state.plant(y)*std::sin(state.plant(psi))));

    dhdx.push_back(T(epz, z, -1));
    dhdx.push_back(T(edz, z, 1/timestep));

    dhdx.push_back(T(state_body_vx, xdot, std::cos(state.plant(psi))));
    dhdx.push_back(T(state_body_vx, ydot, std::sin(state.plant(psi))));
    dhdx.push_back(T(state_body_vx, phi, -state.plant(xdot)*std::sin(state.plant(psi)) + state.plant(ydot)*std::cos(state.plant(psi))));

    dhdx.push_back(T(state_body_vy, xdot, -std::sin(state.plant(psi))));
    dhdx.push_back(T(state_body_vy, ydot, std::cos(state.plant(psi))));
    dhdx.push_back(T(state_body_vy, phi, -state.plant(xdot)*std::cos(state.plant(psi)) - state.plant(ydot)*std::sin(state.plant(psi))));

    dhdx.push_back(T(epzdot, zdot, -1));
    dhdx.push_back(T(edzdot, zdot, 1/timestep));

    dhdx.push_back(T(epphi, phi, -180/M_PI));
    dhdx.push_back(T(edphi, phi, 180/(timestep*M_PI)));

    dhdx.push_back(T(eptheta, theta, -180/M_PI));
    dhdx.push_back(T(edtheta, theta, 180/(timestep*M_PI)));

    dhdx.push_back(T(eppsi, psi, -180/M_PI));
    dhdx.push_back(T(edpsi, psi, 180/(timestep*M_PI)));

    dhdx.push_back(T(delay_1_rollRate, p, 180/(timestep*M_PI)));
    dhdx.push_back(T(delay_1_pitchRate, q, 180/(timestep*M_PI)));

    dhdx.push_back(T(epp, p, -180/M_PI));
    dhdx.push_back(T(edp, p, lpf_b0*180/(timestep*M_PI)));

    dhdx.push_back(T(epq, q, -180/M_PI));
    dhdx.push_back(T(edq, q, lpf_b0*180/(timestep*M_PI)));

    dhdx.push_back(T(epr, r, -180/M_PI));
    dhdx.push_back(T(edr, r, 180/(timestep*M_PI)));

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxCurr(double timestep)
{
    // TODO: what if ref is not constant?
    std::vector<T> dhdx;
    dhdx.reserve(10);
    
    dhdx.push_back(T(edz, z, -1/timestep));
    dhdx.push_back(T(edzdot, zdot, -1/timestep));
    dhdx.push_back(T(edphi, phi, -1/(timestep*M_PI)));
    dhdx.push_back(T(edtheta, theta, -1/(timestep*M_PI)));
    dhdx.push_back(T(edpsi, psi, -1/(timestep*M_PI)));
    dhdx.push_back(T(delay_1_rollRate, p, -180/(timestep*M_PI)));
    dhdx.push_back(T(delay_1_pitchRate, q, -180/(timestep*M_PI)));
    dhdx.push_back(T(edp, p, -lpf_b0*1/(timestep*M_PI)));
    dhdx.push_back(T(epq, q, -lpf_b0*1/(timestep*M_PI)));
    dhdx.push_back(T(edr, r, -1/(timestep*M_PI)));

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzPlus(SystemState state, double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(89);
    
    dhdz.push_back(T(epx, setp_body_x, 1));
    dhdz.push_back(T(epx, state_body_x, -1));
    dhdz.push_back(T(eix, epx, timestep));
    dhdz.push_back(T(desVelX, epx, m_ctrlParams.at(posX).kp));
    dhdz.push_back(T(desVelX, eix, m_ctrlParams.at(posX).ki));
    dhdz.push_back(T(desVelX, edx, m_ctrlParams.at(posX).kd));

    dhdz.push_back(T(epy, setp_body_y, 1));
    dhdz.push_back(T(epy, state_body_y, -1));
    dhdz.push_back(T(eiy, epy, timestep));
    dhdz.push_back(T(desVelY, epy, m_ctrlParams.at(posY).kp));
    dhdz.push_back(T(desVelY, eiy, m_ctrlParams.at(posY).ki));
    dhdz.push_back(T(desVelY, edy, m_ctrlParams.at(posY).kd));

    dhdz.push_back(T(eiz, epz, timestep));
    dhdz.push_back(T(desVelZ, epz, m_ctrlParams.at(posZ).kp));
    dhdz.push_back(T(desVelZ, eiz, m_ctrlParams.at(posZ).ki));
    dhdz.push_back(T(desVelZ, edz, m_ctrlParams.at(posZ).kd));

    dhdz.push_back(T(epxdot, desVelX, 1));
    dhdz.push_back(T(epxdot, state_body_vx, -1));
    dhdz.push_back(T(eixdot, epxdot, timestep));
    dhdz.push_back(T(edxdot, state_body_vx, 1/timestep));

    dhdz.push_back(T(epydot, desVelY, 1));
    dhdz.push_back(T(epydot, state_body_vy, -1));
    dhdz.push_back(T(eiydot, epydot, timestep));
    dhdz.push_back(T(edydot, state_body_vy, 1/timestep));

    dhdz.push_back(T(epzdot, desVelZ, 1));
    dhdz.push_back(T(eizdot, epzdot, timestep));
    dhdz.push_back(T(desThrust, epzdot, m_ctrlParams.at(velZ).kp*m_thrustScale));
    dhdz.push_back(T(desThrust, eizdot, m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(desThrust, edzdot, m_ctrlParams.at(velZ).kd*m_thrustScale));

    dhdz.push_back(T(eiphi, epphi, timestep));
    dhdz.push_back(T(desRollRate, epphi, m_ctrlParams.at(roll).kp));
    dhdz.push_back(T(desRollRate, eiphi, m_ctrlParams.at(roll).ki));
    dhdz.push_back(T(desRollRate, edphi, m_ctrlParams.at(roll).kd));

    dhdz.push_back(T(eitheta, eptheta, timestep));
    dhdz.push_back(T(desPitchRate, eptheta, m_ctrlParams.at(pitch).kp));
    dhdz.push_back(T(desPitchRate, eitheta, m_ctrlParams.at(pitch).ki));
    dhdz.push_back(T(desPitchRate, edtheta, m_ctrlParams.at(pitch).kd));

    dhdz.push_back(T(eipsi, eppsi, timestep));
    dhdz.push_back(T(desYawRate, eppsi, m_ctrlParams.at(yaw).kp));
    dhdz.push_back(T(desYawRate, eipsi, m_ctrlParams.at(yaw).ki));
    dhdz.push_back(T(desYawRate, edpsi, m_ctrlParams.at(yaw).kd));

    dhdz.push_back(T(epp, desRollRate, 1));
    dhdz.push_back(T(eip, epp, timestep));
    dhdz.push_back(T(desRollOutput, epp, m_ctrlParams.at(rollRate).kp));
    dhdz.push_back(T(desRollOutput, eip, m_ctrlParams.at(rollRate).ki));
    dhdz.push_back(T(desRollOutput, edp, m_ctrlParams.at(rollRate).kd));

    dhdz.push_back(T(epq, desPitchRate, 1));
    dhdz.push_back(T(eiq, epq, timestep));
    dhdz.push_back(T(desPitchOutput, epq, m_ctrlParams.at(pitchRate).kp));
    dhdz.push_back(T(desPitchOutput, eiq, m_ctrlParams.at(pitchRate).ki));
    dhdz.push_back(T(desPitchOutput, edq, m_ctrlParams.at(pitchRate).kd));

    dhdz.push_back(T(epr, desYawRate, 1));
    dhdz.push_back(T(eir, epr, timestep));
    dhdz.push_back(T(desYawOutput, epr, m_ctrlParams.at(yawRate).kp));
    dhdz.push_back(T(desYawOutput, eir, m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(desYawOutput, edr, m_ctrlParams.at(yawRate).kd));

    dhdz.push_back(T(w1, desThrust, m_alpha));
    dhdz.push_back(T(w1, desRollOutput, -m_alpha/2));
    dhdz.push_back(T(w1, desPitchOutput, m_alpha/2));
    dhdz.push_back(T(w1, desYawOutput, m_alpha));

    dhdz.push_back(T(w2, desThrust, m_alpha));
    dhdz.push_back(T(w2, desRollOutput, -m_alpha/2));
    dhdz.push_back(T(w2, desPitchOutput, -m_alpha/2));
    dhdz.push_back(T(w2, desYawOutput, -m_alpha));

    dhdz.push_back(T(w3, desThrust, m_alpha));
    dhdz.push_back(T(w3, desRollOutput, m_alpha/2));
    dhdz.push_back(T(w3, desPitchOutput, -m_alpha/2));
    dhdz.push_back(T(w3, desYawOutput, m_alpha));

    dhdz.push_back(T(w4, desThrust, m_alpha));
    dhdz.push_back(T(w4, desRollOutput, m_alpha/2));
    dhdz.push_back(T(w4, desPitchOutput, m_alpha/2));
    dhdz.push_back(T(w4, desYawOutput, -m_alpha));
    
    dhdz.push_back(T(ft, w1, m_droneParams.kf * 2 *state.alge(w1)));
    dhdz.push_back(T(ft, w2, m_droneParams.kf * 2 *state.alge(w2)));
    dhdz.push_back(T(ft, w3, m_droneParams.kf * 2 *state.alge(w3)));
    dhdz.push_back(T(ft, w4, m_droneParams.kf * 2 *state.alge(w4)));

    dhdz.push_back(T(tx, w1, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w1)));
    dhdz.push_back(T(tx, w2, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w2)));
    dhdz.push_back(T(tx, w3, m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w3)));
    dhdz.push_back(T(tx, w4, m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w4)));

    dhdz.push_back(T(ty, w1, m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w1)));
    dhdz.push_back(T(ty, w2, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w2)));
    dhdz.push_back(T(ty, w3, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w3)));
    dhdz.push_back(T(ty, w4, m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * state.alge(w4)));

    dhdz.push_back(T(tz, w1, m_droneParams.km * state.alge(w1)));
    dhdz.push_back(T(tz, w2, - m_droneParams.km * state.alge(w2)));
    dhdz.push_back(T(tz, w3, m_droneParams.km * state.alge(w3)));
    dhdz.push_back(T(tz, w4, - m_droneParams.km * state.alge(w4)));

    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzCurr(double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(23);
    dhdz.push_back(T(eix, eix, 1));
    dhdz.push_back(T(eiy, eiy, 1));
    dhdz.push_back(T(eiz, eiz, 1));
    dhdz.push_back(T(eixdot, eixdot, 1));
    dhdz.push_back(T(edxdot, state_body_vx, -1/timestep));

    dhdz.push_back(T(eiydot, eiydot, 1));
    dhdz.push_back(T(edydot, state_body_vy, -1/timestep));

    dhdz.push_back(T(eizdot, eizdot, 1));
    dhdz.push_back(T(eiphi, eiphi, 1));
    dhdz.push_back(T(eitheta, eitheta, 1));
    dhdz.push_back(T(eipsi, eipsi, 1));
    
    dhdz.push_back(T(delay_1_rollRate, delay_1_rollRate, -lpf_a1));
    dhdz.push_back(T(delay_1_rollRate, delay_2_rollRate, -lpf_a2));
    dhdz.push_back(T(delay_2_rollRate, delay_1_rollRate, -1));
    
    dhdz.push_back(T(delay_1_pitchRate, delay_1_pitchRate, -lpf_a1));
    dhdz.push_back(T(delay_1_pitchRate, delay_2_pitchRate, -lpf_a2));
    dhdz.push_back(T(delay_2_pitchRate, delay_1_pitchRate, -1));

    dhdz.push_back(T(eip, eip, 1));
    dhdz.push_back(T(edp, delay_1_rollRate, -lpf_b0*lpf_a1+lpf_b1));
    dhdz.push_back(T(edp, delay_2_rollRate, -lpf_b0*lpf_a2+lpf_b2));
    dhdz.push_back(T(eiq, eiq, 1));
    dhdz.push_back(T(edq, delay_1_pitchRate, -lpf_b0*lpf_a1+lpf_b1));
    dhdz.push_back(T(edq, delay_2_pitchRate, -lpf_b0*lpf_a2+lpf_b2));
    dhdz.push_back(T(eir, eir, 1));

    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdy()
{
    Eigen::SparseMatrix<double> dhdy(NUM_Z_STATES, NUM_Y_STATES);
    dhdy.insert(epphi, 0) = 1;
    dhdy.insert(eptheta, 1) = 1;
    return dhdy; 
}

