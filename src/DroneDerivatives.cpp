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

    dfdx(3,3) = plant(q)*std::cos(plant(phi))*tan(plant(theta)) 
              - plant(r)*std::sin(plant(phi))*tan(plant(theta));
    dfdx(3,4) = plant(r)*std::cos(plant(phi))*(pow(tan(plant(theta)), 2) + 1) 
              + plant(q)*std::sin(plant(phi))*(pow(tan(plant(theta)), 2) + 1);
    dfdx(3,9) = 1;
    dfdx(3,10) = std::sin(plant(phi))*tan(plant(theta));
    dfdx(3,11) = std::cos(plant(phi))*tan(plant(theta));
    dfdx(4,3) = - plant(r)*std::cos(plant(phi)) - plant(q)*std::sin(plant(phi));
    dfdx(4,10) = std::cos(plant(phi));
    dfdx(4,11) = -std::sin(plant(phi));
    dfdx(5,3) = (plant(q)*std::cos(plant(phi)))/std::cos(plant(theta)) 
              - (plant(r)*std::sin(plant(phi)))/std::cos(plant(theta));
    dfdx(5,4) = (plant(r)*std::cos(plant(phi))*std::sin(plant(theta)))/pow(std::cos(plant(theta)), 2) 
              + (plant(q)*std::sin(plant(phi))*std::sin(plant(theta)))/pow(std::cos(plant(theta)), 2);
    dfdx(5,10) = std::sin(plant(phi))/std::cos(plant(theta));
    dfdx(5,11) = std::cos(plant(phi))/std::cos(plant(theta));

    dfdx(6,3) = (alge(ft)*(std::cos(plant(phi))*std::sin(plant(psi)) - std::cos(plant(psi))*std::sin(plant(phi))*std::sin(plant(theta))))/m_droneParams.mass;
    dfdx(6,4) = (alge(ft)* std::cos(plant(phi))*std::cos(plant(psi))*std::cos(plant(theta)))/m_droneParams.mass;
    dfdx(6,5) = (alge(ft)*(std::cos(plant(psi))*std::sin(plant(phi)) - std::cos(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;

    dfdx(7,3) = -(alge(ft)*(std::cos(plant(phi))*std::cos(plant(psi)) + std::sin(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;
    dfdx(7,4) = (alge(ft)* std::cos(plant(phi))*std::cos(plant(theta))*std::sin(plant(psi)))/m_droneParams.mass;
    dfdx(7,5) = (alge(ft)*(std::sin(plant(phi))*std::sin(plant(psi)) + std::cos(plant(phi))*std::cos(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;

    dfdx(8,3) = -(alge(ft)*std::cos(plant(theta))*std::sin(plant(phi)))/m_droneParams.mass;
    dfdx(8,4) = -(alge(ft)*std::cos(plant(phi))*std::sin(plant(theta)))/m_droneParams.mass;

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

    dfdz.push_back(T(xdot, ft, 1.0/m_droneParams.mass * 
                (std::sin(state.plant(phi)) * std::sin(state.plant(psi)) 
               + std::cos(state.plant(phi))*std::cos(state.plant(psi))*std::sin(state.plant(theta)))));

    dfdz.push_back(T(ydot, ft, 1.0/m_droneParams.mass * 
                (std::cos(state.plant(phi))*std::sin(state.plant(psi))*std::sin(state.plant(theta)) 
               - std::cos(state.plant(psi))*std::sin(state.plant(phi)))));

    dfdz.push_back(T(zdot, ft, 1.0/m_droneParams.mass * 
                (std::cos(state.plant(phi)) * std::cos(state.plant(theta)))));

    dfdz.push_back(T(p, tx, 1.0/m_droneParams.Ix));
    dfdz.push_back(T(q, ty, 1.0/m_droneParams.Iy));
    dfdz.push_back(T(r, tz, 1.0/m_droneParams.Iz));

    Eigen::SparseMatrix<double> dfdz_mat(NUM_PLANT_STATES, NUM_Z_STATES);
    dfdz_mat.setFromTriplets(dfdz.begin(), dfdz.end());
    return dfdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dgdx(SystemState state)
{
    std::vector<T> dgdx;
    dgdx.reserve(6);
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdx.push_back(T(0, xdot, -m_ctrlParams.at(velY).kp*sin(state.plant(psi))));
        dgdx.push_back(T(0, ydot, m_ctrlParams.at(velY).kp*cos(state.plant(psi))));
        dgdx.push_back(T(0, psi, -m_ctrlParams.at(velY).kp*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdx.push_back(T(1, xdot, -m_ctrlParams.at(velX).kp*cos(state.plant(psi))));
        dgdx.push_back(T(1, ydot, -m_ctrlParams.at(velX).kp*sin(state.plant(psi))));
        dgdx.push_back(T(1, psi, -m_ctrlParams.at(velX).kp*(-state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi)))));
    }
    Eigen::SparseMatrix<double> dgdx_mat(NUM_Y_STATES, NUM_PLANT_STATES);
    dgdx_mat.setFromTriplets(dgdx.begin(), dgdx.end());
    return dgdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dgdz(SystemState state)
{
    std::vector<T> dgdz;
    dgdz.reserve(6);
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdz.push_back(T(0, desVelY, -m_ctrlParams.at(velY).kp));
        dgdz.push_back(T(0, eiydot, -m_ctrlParams.at(velY).ki));
        dgdz.push_back(T(0, edydot, -m_ctrlParams.at(velY).kd));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdz.push_back(T(1, desVelX, m_ctrlParams.at(velX).kp));
        dgdz.push_back(T(1, eixdot, m_ctrlParams.at(velX).ki));
        dgdz.push_back(T(1, edxdot, m_ctrlParams.at(velX).kd));
    }
    Eigen::SparseMatrix<double> dgdz_mat(NUM_Y_STATES, NUM_Z_STATES);
    dgdz_mat.setFromTriplets(dgdz.begin(), dgdz.end());
    return dgdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxPlus(SystemState state, double time, double timestep)
{
    std::vector<T> dhdx;
    dhdx.reserve(150);
    dhdx.push_back(T(eix, x, -timestep*cos(state.plant(psi))));
    dhdx.push_back(T(eix, y, -timestep*sin(state.plant(psi))));
    dhdx.push_back(T(eix, psi, timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))));
    dhdx.push_back(T(edx, x, -1/timestep*cos(state.plant(psi))));
    dhdx.push_back(T(edx, y, -1/timestep*sin(state.plant(psi))));
    dhdx.push_back(T(edx, psi, -(state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, x, -m_ctrlParams.at(posX).kp*cos(state.plant(psi))));
    dhdx.push_back(T(desVelX, y, -m_ctrlParams.at(posX).kp*sin(state.plant(psi))));
    dhdx.push_back(T(desVelX, psi, m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) ));

    dhdx.push_back(T(eiy, x, timestep*sin(state.plant(psi))));
    dhdx.push_back(T(eiy, y, -timestep*cos(state.plant(psi))));
    dhdx.push_back(T(eiy, psi, -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)))));
    dhdx.push_back(T(edy, x, 1/timestep*sin(state.plant(psi))));
    dhdx.push_back(T(edy, y, -1/timestep*cos(state.plant(psi))));
    dhdx.push_back(T(edy, psi, (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, x, m_ctrlParams.at(posY).kp*sin(state.plant(psi))));
    dhdx.push_back(T(desVelY, y, -m_ctrlParams.at(posY).kp*cos(state.plant(psi))));
    dhdx.push_back(T(desVelY, psi, -m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) ));

    dhdx.push_back(T(eiz, z, -timestep));
    dhdx.push_back(T(edz, z, -1/timestep));
    dhdx.push_back(T(desVelZ, z, -m_ctrlParams.at(posZ).kp));

    dhdx.push_back(T(eixdot, xdot, -timestep*cos(state.plant(psi))));
    dhdx.push_back(T(eixdot, ydot, -timestep*sin(state.plant(psi))));
    dhdx.push_back(T(eixdot, psi, -timestep*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))));
    dhdx.push_back(T(edxdot, xdot, -1/timestep*cos(state.plant(psi))));
    dhdx.push_back(T(edxdot, ydot, -1/timestep*sin(state.plant(psi))));
    dhdx.push_back(T(edxdot, psi, -(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))/timestep ));

    dhdx.push_back(T(eiydot, xdot, timestep*sin(state.plant(psi))));
    dhdx.push_back(T(eiydot, ydot, -timestep*cos(state.plant(psi))));
    dhdx.push_back(T(eiydot, psi, timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi))) ));
    dhdx.push_back(T(edydot, xdot, 1/timestep*sin(state.plant(psi))));
    dhdx.push_back(T(edydot, ydot, -1/timestep*cos(state.plant(psi))));
    dhdx.push_back(T(edydot, psi, (state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))/timestep ));

    dhdx.push_back(T(eizdot, zdot, -timestep ));
    dhdx.push_back(T(edzdot, zdot, -1/timestep ));

    dhdx.push_back(T(desThrust, zdot, -m_thrustScale*m_ctrlParams.at(velZ).kp ));

    dhdx.push_back(T(eiphi, phi, -180/M_PI*timestep ));
    dhdx.push_back(T(edphi, phi, -180/M_PI/timestep ));
    dhdx.push_back(T(desRollRate, phi, -m_ctrlParams.at(roll).kp*180/M_PI));

    dhdx.push_back(T(eitheta, theta, -180/M_PI*timestep ));
    dhdx.push_back(T(edtheta, theta, -180/M_PI/timestep ));
    dhdx.push_back(T(desPitchRate, theta, -m_ctrlParams.at(pitch).kp*180/M_PI));

    dhdx.push_back(T(eipsi, psi, -180/M_PI*timestep ));
    dhdx.push_back(T(edpsi, psi, -180/M_PI/timestep ));
    dhdx.push_back(T(desYawRate, psi, -m_ctrlParams.at(yaw).kp*180/M_PI));

    dhdx.push_back(T(delay_1_rollRate, p, -180/timestep/M_PI));
    dhdx.push_back(T(delay_1_pitchRate, q, -180/timestep/M_PI));

    dhdx.push_back(T(eip, p, -180/M_PI*timestep ));
    dhdx.push_back(T(edp, p, -180/M_PI/timestep*lpf_b0 ));
    dhdx.push_back(T(desRollOutput, p, -m_ctrlParams.at(rollRate).kp*180/M_PI));

    dhdx.push_back(T(eiq, q, -180/M_PI*timestep ));
    dhdx.push_back(T(edq, q, -180/M_PI/timestep*lpf_b0 ));
    dhdx.push_back(T(desPitchOutput, q, -m_ctrlParams.at(pitchRate).kp*180/M_PI));

    dhdx.push_back(T(eir, r, -180/M_PI*timestep ));
    dhdx.push_back(T(edr, r, -180/M_PI/timestep ));
    dhdx.push_back(T(desYawOutput, r, -m_ctrlParams.at(yawRate).kp*180/M_PI));

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxCurr(SystemState prev, double timestep)
{
    std::vector<T> dhdx;
    dhdx.reserve(150);
    dhdx.push_back(T(edx, x, 1/timestep*cos(prev.plant(psi))));
    dhdx.push_back(T(edx, y, 1/timestep*sin(prev.plant(psi))));
    dhdx.push_back(T(edx, psi, (prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep));

    dhdx.push_back(T(edy, x, -1/timestep*sin(prev.plant(psi))));
    dhdx.push_back(T(edy, y, 1/timestep*cos(prev.plant(psi))));
    dhdx.push_back(T(edy, psi, -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep ));

    dhdx.push_back(T(edz, z, 1/timestep));

    dhdx.push_back(T(edxdot, xdot, 1/timestep*cos(prev.plant(psi))));
    dhdx.push_back(T(edxdot, ydot, 1/timestep*sin(prev.plant(psi))));
    dhdx.push_back(T(edxdot, psi, (prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep ));

    dhdx.push_back(T(edydot, xdot, -1/timestep*sin(prev.plant(psi))));
    dhdx.push_back(T(edydot, ydot, 1/timestep*cos(prev.plant(psi))));
    dhdx.push_back(T(edydot, psi, -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep ));

    dhdx.push_back(T(edzdot, zdot, 1/timestep));

    dhdx.push_back(T(edphi, phi, 180/M_PI/timestep));
    dhdx.push_back(T(edtheta, theta, 180/M_PI/timestep));
    dhdx.push_back(T(edpsi, psi, 180/M_PI/timestep));

    dhdx.push_back(T(delay_1_rollRate, p, 180/M_PI/timestep));
    dhdx.push_back(T(delay_1_pitchRate, q, 180/M_PI/timestep));

    dhdx.push_back(T(edp, p, 180/M_PI/timestep*lpf_b0));
    dhdx.push_back(T(edq, q, 180/M_PI/timestep*lpf_b0));
    dhdx.push_back(T(edr, r, 180/M_PI/timestep));

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzPlus(SystemState state, double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(89);
    
    dhdz.push_back(T(desVelX, eix, m_ctrlParams.at(posX).ki));
    dhdz.push_back(T(desVelX, edx, m_ctrlParams.at(posX).kd));

    dhdz.push_back(T(desVelY, eiy, m_ctrlParams.at(posY).ki));
    dhdz.push_back(T(desVelY, edy, m_ctrlParams.at(posY).kd));

    dhdz.push_back(T(desVelZ, eiz, m_ctrlParams.at(posZ).ki));
    dhdz.push_back(T(desVelZ, edz, m_ctrlParams.at(posZ).kd));

    dhdz.push_back(T(eixdot, desVelX, timestep));
    dhdz.push_back(T(eiydot, desVelY, timestep));
    dhdz.push_back(T(eizdot, desVelZ, timestep));

    dhdz.push_back(T(desThrust, desVelZ, m_ctrlParams.at(velZ).kp*m_thrustScale));
    dhdz.push_back(T(desThrust, eizdot, m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(desThrust, edzdot, m_ctrlParams.at(velZ).kd*m_thrustScale));

    dhdz.push_back(T(desRollRate, eiphi, m_ctrlParams.at(roll).ki));
    dhdz.push_back(T(desRollRate, edphi, m_ctrlParams.at(roll).kd));

    dhdz.push_back(T(desPitchRate, eitheta, m_ctrlParams.at(pitch).ki));
    dhdz.push_back(T(desPitchRate, edtheta, m_ctrlParams.at(pitch).kd));
    
    dhdz.push_back(T(desYawRate, eipsi, m_ctrlParams.at(yaw).ki));
    dhdz.push_back(T(desYawRate, edpsi, m_ctrlParams.at(yaw).kd));
    
    dhdz.push_back(T(eip, desRollRate, timestep));
    dhdz.push_back(T(desRollOutput, desRollRate, m_ctrlParams.at(rollRate).kp));
    dhdz.push_back(T(desRollOutput, eip, m_ctrlParams.at(rollRate).ki));
    dhdz.push_back(T(desRollOutput, edp, m_ctrlParams.at(rollRate).kd));

    dhdz.push_back(T(eiq, desPitchRate, timestep));
    dhdz.push_back(T(desPitchOutput, desPitchRate, m_ctrlParams.at(pitchRate).kp));
    dhdz.push_back(T(desPitchOutput, eiq, m_ctrlParams.at(pitchRate).ki));
    dhdz.push_back(T(desPitchOutput, edq, m_ctrlParams.at(pitchRate).kd));

    dhdz.push_back(T(eir, desYawRate, timestep));
    dhdz.push_back(T(desYawOutput, desYawRate, m_ctrlParams.at(yawRate).kp));
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

    dhdz.push_back(T(tz, w1, m_droneParams.km * 2 * state.alge(w1)));
    dhdz.push_back(T(tz, w2, - m_droneParams.km * 2 * state.alge(w2)));
    dhdz.push_back(T(tz, w3, m_droneParams.km * 2 * state.alge(w3)));
    dhdz.push_back(T(tz, w4, - m_droneParams.km * 2 * state.alge(w4)));

    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzCurr()
{
    std::vector<T> dhdz;
    dhdz.reserve(93);
    dhdz.push_back(T(eix, eix, 1));
    dhdz.push_back(T(eiy, eiy, 1));
    dhdz.push_back(T(eiz, eiz, 1));
    dhdz.push_back(T(eixdot, eixdot, 1));
    dhdz.push_back(T(eiydot, eiydot, 1));
    dhdz.push_back(T(eizdot, eizdot, 1));
    dhdz.push_back(T(eiphi, eiphi, 1));
    dhdz.push_back(T(eitheta, eitheta, 1));
    dhdz.push_back(T(eipsi, eipsi, 1));
    dhdz.push_back(T(eip, eip, 1));
    dhdz.push_back(T(edp, delay_1_rollRate, lpf_b1 - lpf_a1*lpf_b0));
    dhdz.push_back(T(edp, delay_2_rollRate, lpf_b2 - lpf_b0*lpf_a2));
    dhdz.push_back(T(eiq, eiq, 1));
    dhdz.push_back(T(edq, delay_1_pitchRate, lpf_b1 - lpf_a1*lpf_b0));
    dhdz.push_back(T(edq, delay_2_pitchRate, lpf_b2 - lpf_b0*lpf_a2));
    dhdz.push_back(T(eir, eir, 1));
    dhdz.push_back(T(delay_1_rollRate, delay_1_rollRate, -lpf_a1));
    dhdz.push_back(T(delay_1_rollRate, delay_2_rollRate, -lpf_a2));
    dhdz.push_back(T(delay_2_rollRate, delay_1_rollRate, 1));
    dhdz.push_back(T(delay_1_pitchRate, delay_1_pitchRate, -lpf_a1));
    dhdz.push_back(T(delay_1_pitchRate, delay_2_pitchRate, -lpf_a2));
    dhdz.push_back(T(delay_2_pitchRate, delay_1_pitchRate, 1));

    
    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdy(double timestep)
{
    std::vector<T> dhdy;
   
    dhdy.push_back(T(eiphi, 0, timestep));
    dhdy.push_back(T(desRollRate, 0, m_ctrlParams.at(roll).kp));
    dhdy.push_back(T(eitheta, 1, timestep));
    dhdy.push_back(T(desPitchRate, 1, m_ctrlParams.at(pitch).kp ));
    
    Eigen::SparseMatrix<double> dhdy_mat(NUM_Z_STATES, NUM_Y_STATES);
    dhdy_mat.setFromTriplets(dhdy.begin(), dhdy.end()); 
    return dhdy_mat;
}

Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::d2fdx2(SystemState state)
{
    Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> d2fdx2 = Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES>::Zero();
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;

    d2fdx2(39, 3) = - plant(r)*cos(plant(phi))*tan(plant(theta)) - plant(q)*sin(plant(phi))*tan(plant(theta));
    d2fdx2(39, 4) = plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(39, 10) = cos(plant(phi))*tan(plant(theta));
    d2fdx2(39, 11) = -sin(plant(phi))*tan(plant(theta));
    d2fdx2(40, 3) = plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(40, 4) = 2*plant(r)*cos(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1) + 2*plant(q)*sin(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(40, 10) = sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(40, 11) = cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(46, 3) = cos(plant(phi))*tan(plant(theta));
    d2fdx2(46, 4) = sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(47, 3) = -sin(plant(phi))*tan(plant(theta));
    d2fdx2(47, 4) = cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2(51, 3) = plant(r)*sin(plant(phi)) - plant(q)*cos(plant(phi));
    d2fdx2(51, 10) = -sin(plant(phi));
    d2fdx2(51, 11) = -cos(plant(phi));
    d2fdx2(58, 3) = -sin(plant(phi));
    d2fdx2(59, 3) = -cos(plant(phi));
    d2fdx2(63, 3) = - (plant(r)*cos(plant(phi)))/cos(plant(theta)) - (plant(q)*sin(plant(phi)))/cos(plant(theta));
    d2fdx2(63, 4) = (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(63, 10) = cos(plant(phi))/cos(plant(theta));
    d2fdx2(63, 11) = -sin(plant(phi))/cos(plant(theta));
    d2fdx2(64, 3) = (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(64, 4) = (plant(r)*cos(plant(phi)))/cos(plant(theta)) + (plant(q)*sin(plant(phi)))/cos(plant(theta)) + std::pow((2*plant(r)*cos(plant(phi))*std::pow(sin(plant(theta)), 2))/cos(plant(theta)), 3) + std::pow((2*plant(q)*sin(plant(phi))*std::pow(sin(plant(theta)), 2))/cos(plant(theta)), 3);
    d2fdx2(64, 10) = (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(64, 11) = (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(70, 3) = cos(plant(phi))/cos(plant(theta));
    d2fdx2(70, 4) = (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(71, 3) = -sin(plant(phi))/cos(plant(theta));
    d2fdx2(71, 4) = (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2(75, 3) = -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(75, 4) = -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    d2fdx2(75, 5) = (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(76, 3) = -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    d2fdx2(76, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(psi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2(76, 5) = -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2(77, 3) = (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(77, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2(77, 5) = -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(87, 3) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(87, 4) = -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2(87, 5) = (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(88, 3) = -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2(88, 4) = -(alge(ft)*cos(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2(88, 5) = (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2(89, 3) = (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(89, 4) = (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2(89, 5) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2(99, 3) = -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2(99, 4) = (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2(100, 3) = (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2(100, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2(118, 11) = (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix;
    d2fdx2(119, 10) = (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix;
    d2fdx2(129, 11) = -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy;
    d2fdx2(131, 9) = -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy;
    d2fdx2(141, 10) = (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz;
    d2fdx2(142, 9) = (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz;
    return d2fdx2;
}


Eigen::Matrix<double, NUM_Y_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::d2gdx2(SystemState state)
{
    Eigen::Matrix<double, NUM_Y_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> d2gdx2 = Eigen::Matrix<double, NUM_Y_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES>::Zero();
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        d2gdx2(xdot, psi) = -m_ctrlParams.at(velY).kp*cos(state.plant(psi));
        d2gdx2(ydot, psi) = -m_ctrlParams.at(velY).kp*sin(state.plant(psi));
        d2gdx2(psi, xdot) = -m_ctrlParams.at(velY).kp*cos(state.plant(psi));
        d2gdx2(psi, ydot) = -m_ctrlParams.at(velY).kp*sin(state.plant(psi));
        d2gdx2(psi, psi) =  -m_ctrlParams.at(velY).kp*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        d2gdx2(NUM_PLANT_STATES + xdot, psi) =  m_ctrlParams.at(velX).kp*sin(state.plant(psi));
        d2gdx2(NUM_PLANT_STATES + ydot, psi) = -m_ctrlParams.at(velX).kp*cos(state.plant(psi));
        d2gdx2(NUM_PLANT_STATES + psi, xdot) =  m_ctrlParams.at(velX).kp*sin(state.plant(psi));
        d2gdx2(NUM_PLANT_STATES + psi, ydot) = -m_ctrlParams.at(velX).kp*cos(state.plant(psi));
        d2gdx2(NUM_PLANT_STATES + psi, psi) = m_ctrlParams.at(velX).kp*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)));
    }
    return d2gdx2;
}

Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::d2hdx2_plus(SystemState state, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> d2hdx2_plus = Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES>::Zero();
    d2hdx2_plus(eix*NUM_PLANT_STATES+x, psi) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eix*NUM_PLANT_STATES+y, psi) = -timestep*cos(state.plant(psi));
    d2hdx2_plus(eix*NUM_PLANT_STATES+psi, x) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eix*NUM_PLANT_STATES+psi, y) = -timestep*cos(state.plant(psi));
    d2hdx2_plus(eix*NUM_PLANT_STATES+psi, psi) = -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)));

    d2hdx2_plus(edx*NUM_PLANT_STATES+x, psi) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edx*NUM_PLANT_STATES+y, psi) = -1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edx*NUM_PLANT_STATES+psi, x) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edx*NUM_PLANT_STATES+psi, y) = -1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edx*NUM_PLANT_STATES+psi, psi) = (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep;

    d2hdx2_plus(desVelX*NUM_PLANT_STATES+x, psi) = m_ctrlParams.at(posX).kp*sin(state.plant(psi));
    d2hdx2_plus(desVelX*NUM_PLANT_STATES+y, psi) = -m_ctrlParams.at(posX).kp*cos(state.plant(psi));
    d2hdx2_plus(desVelX*NUM_PLANT_STATES+psi, x) = m_ctrlParams.at(posX).kp*sin(state.plant(psi));
    d2hdx2_plus(desVelX*NUM_PLANT_STATES+psi, y) = -m_ctrlParams.at(posX).kp*cos(state.plant(psi));
    d2hdx2_plus(desVelX*NUM_PLANT_STATES+psi, psi) = -m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)));

    d2hdx2_plus(eiy*NUM_PLANT_STATES+x, psi) = timestep*cos(state.plant(psi));
    d2hdx2_plus(eiy*NUM_PLANT_STATES+y, psi) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eiy*NUM_PLANT_STATES+psi, x) = timestep*cos(state.plant(psi));
    d2hdx2_plus(eiy*NUM_PLANT_STATES+psi, y) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eiy*NUM_PLANT_STATES+psi, psi) = -timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -x));

    d2hdx2_plus(edy*NUM_PLANT_STATES+x, psi) = 1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edy*NUM_PLANT_STATES+y, psi) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edy*NUM_PLANT_STATES+psi, x) = 1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edy*NUM_PLANT_STATES+psi, y) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edy*NUM_PLANT_STATES+psi, psi) = (state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi)))/timestep;

    d2hdx2_plus(desVelY*NUM_PLANT_STATES+x, psi) = m_ctrlParams.at(posY).kp*cos(state.plant(psi));
    d2hdx2_plus(desVelY*NUM_PLANT_STATES+y, psi) = m_ctrlParams.at(posY).kp*sin(state.plant(psi));
    d2hdx2_plus(desVelY*NUM_PLANT_STATES+psi, x) = m_ctrlParams.at(posY).kp*cos(state.plant(psi));
    d2hdx2_plus(desVelY*NUM_PLANT_STATES+psi, y) = m_ctrlParams.at(posY).kp*sin(state.plant(psi));
    d2hdx2_plus(desVelY*NUM_PLANT_STATES+psi, psi) = -m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)));
    
    d2hdx2_plus(eixdot*NUM_PLANT_STATES+xdot, psi) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eixdot*NUM_PLANT_STATES+ydot, psi) = -timestep*cos(state.plant(psi));
    d2hdx2_plus(eixdot*NUM_PLANT_STATES+psi, x) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eixdot*NUM_PLANT_STATES+psi, y) = -timestep*cos(state.plant(psi));
    d2hdx2_plus(eixdot*NUM_PLANT_STATES+psi, psi) = timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)));

    d2hdx2_plus(edxdot*NUM_PLANT_STATES+xdot, psi) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edxdot*NUM_PLANT_STATES+ydot, psi) = -1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edxdot*NUM_PLANT_STATES+psi, xdot) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edxdot*NUM_PLANT_STATES+psi, ydot) = -1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edxdot*NUM_PLANT_STATES+psi, psi) = (state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))/timestep;

    d2hdx2_plus(eiydot*NUM_PLANT_STATES+xdot, psi) = timestep*cos(state.plant(psi));
    d2hdx2_plus(eiydot*NUM_PLANT_STATES+ydot, psi) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eiydot*NUM_PLANT_STATES+psi, xdot) = timestep*cos(state.plant(psi));
    d2hdx2_plus(eiydot*NUM_PLANT_STATES+psi, ydot) = timestep*sin(state.plant(psi));
    d2hdx2_plus(eiydot*NUM_PLANT_STATES+psi, psi) = timestep*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)));

    d2hdx2_plus(edydot*NUM_PLANT_STATES+xdot, psi) = 1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edydot*NUM_PLANT_STATES+ydot, psi) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edydot*NUM_PLANT_STATES+psi, xdot) = 1/timestep*cos(state.plant(psi));
    d2hdx2_plus(edydot*NUM_PLANT_STATES+psi, ydot) = 1/timestep*sin(state.plant(psi));
    d2hdx2_plus(edydot*NUM_PLANT_STATES+psi, psi) = (state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))/timestep;
  
    return d2hdx2_plus;
}

Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::d2hdx2_curr(SystemState prev, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES> d2hdx2_curr = Eigen::Matrix<double, NUM_Z_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES>::Zero();
    d2hdx2_curr(edx*NUM_PLANT_STATES + x, psi) = -1/timestep*sin(prev.plant(psi));
    d2hdx2_curr(edx*NUM_PLANT_STATES + y, psi) = 1/timestep*cos(prev.plant(psi));
    d2hdx2_curr(edx*NUM_PLANT_STATES + psi, x) = -sin(prev.plant(psi))/timestep;
    d2hdx2_curr(edx*NUM_PLANT_STATES + psi, y) = cos(prev.plant(psi))/timestep;
    d2hdx2_curr(edx*NUM_PLANT_STATES + psi, psi) = -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep;

    d2hdx2_curr(edy*NUM_PLANT_STATES + x, psi) = -1/timestep*cos(prev.plant(psi));
    d2hdx2_curr(edy*NUM_PLANT_STATES + y, psi) = -1/timestep*sin(prev.plant(psi));
    d2hdx2_curr(edy*NUM_PLANT_STATES + psi, x) = -cos(prev.plant(psi))/timestep;
    d2hdx2_curr(edy*NUM_PLANT_STATES + psi, y) = -sin(prev.plant(psi))/timestep;
    d2hdx2_curr(edy*NUM_PLANT_STATES + psi, psi) = -(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep;

    d2hdx2_curr(edxdot*NUM_PLANT_STATES + xdot, psi) = -1/timestep*sin(prev.plant(psi));
    d2hdx2_curr(edxdot*NUM_PLANT_STATES + ydot, psi) = 1/timestep*cos(prev.plant(psi));
    d2hdx2_curr(edxdot*NUM_PLANT_STATES + psi, x) = -sin(prev.plant(psi))/timestep;
    d2hdx2_curr(edxdot*NUM_PLANT_STATES + psi, y) = cos(prev.plant(psi))/timestep;
    d2hdx2_curr(edxdot*NUM_PLANT_STATES + psi, psi) = -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep;

    d2hdx2_curr(edydot*NUM_PLANT_STATES + xdot, psi) = -1/timestep*cos(prev.plant(psi));
    d2hdx2_curr(edydot*NUM_PLANT_STATES + ydot, psi) = -1/timestep*sin(prev.plant(psi));
    d2hdx2_curr(edydot*NUM_PLANT_STATES + psi, x) = -cos(prev.plant(psi))/timestep;
    d2hdx2_curr(edydot*NUM_PLANT_STATES + psi, y) = -sin(prev.plant(psi))/timestep;
    d2hdx2_curr(edydot*NUM_PLANT_STATES + psi, psi) = -(prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep;

    return d2hdx2_curr;
}

Eigen::Matrix<double, NUM_Z_STATES*NUM_Z_STATES, NUM_Z_STATES> DroneTrajectory::d2hdz2_plus()
{
    Eigen::Matrix<double, NUM_Z_STATES*NUM_Z_STATES, NUM_Z_STATES> d2hdz2_plus = Eigen::Matrix<double, NUM_Z_STATES*NUM_Z_STATES, NUM_Z_STATES>::Zero();
  
    d2hdz2_plus(ft*NUM_Z_STATES + w1, w1) = m_droneParams.kf * 2;
    d2hdz2_plus(ft*NUM_Z_STATES + w2, w2) = m_droneParams.kf * 2;
    d2hdz2_plus(ft*NUM_Z_STATES + w3, w3) = m_droneParams.kf * 2;
    d2hdz2_plus(ft*NUM_Z_STATES + w4, w4) = m_droneParams.kf * 2;

    d2hdz2_plus(tx*NUM_Z_STATES + w1, w1) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(tx*NUM_Z_STATES + w2, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(tx*NUM_Z_STATES + w3, w3) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(tx*NUM_Z_STATES + w4, w4) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);

    d2hdz2_plus(ty*NUM_Z_STATES + w1, w1) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(ty*NUM_Z_STATES + w2, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(ty*NUM_Z_STATES + w3, w3) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    d2hdz2_plus(ty*NUM_Z_STATES + w4, w4) = m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2);
    
    d2hdz2_plus(tz*NUM_Z_STATES + w1, w1) = m_droneParams.km * 2;
    d2hdz2_plus(tz*NUM_Z_STATES + w2, w2) = - m_droneParams.km * 2;
    d2hdz2_plus(tz*NUM_Z_STATES + w3, w3) = m_droneParams.km * 2;
    d2hdz2_plus(tz*NUM_Z_STATES + w4, w4) = - m_droneParams.km * 2;

    return d2hdz2_plus;
}
