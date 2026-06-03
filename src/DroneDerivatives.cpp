#include "DroneTrajectory.h"
typedef Eigen::Triplet<double> T;

Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> DroneTrajectory::dfdx(SystemState state) const
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

Eigen::SparseMatrix<double> DroneTrajectory::dfdz(SystemState state) const
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

Eigen::SparseMatrix<double> DroneTrajectory::dgdx(SystemState state) const
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

Eigen::SparseMatrix<double> DroneTrajectory::dgdz(SystemState state) const
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

Eigen::SparseMatrix<double> DroneTrajectory::dgdp(SystemState state)
{
    std::vector<T> dgdp;
    dgdp.reserve(6);
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdp.push_back(T(0, kpvy, -(state.alge(desVelY) + state.plant(xdot)*std::sin(state.plant(psi)) - state.plant(ydot)*std::cos(state.plant(psi)) ) ));
        dgdp.push_back(T(0, kivy, -state.alge(eiydot)));
        dgdp.push_back(T(0, kdvy, -state.alge(edydot)));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdp.push_back(T(0, kpvz, (state.alge(desVelX) - state.plant(xdot)*std::cos(state.plant(psi)) - state.plant(ydot)*std::sin(state.plant(psi)) ) ));
        dgdp.push_back(T(0, kivz, state.alge(eixdot)));
        dgdp.push_back(T(0, kdvz, state.alge(edxdot)));
    }
    Eigen::SparseMatrix<double> dgdp_mat(NUM_Y_STATES, NUM_PARAMETERS);
    dgdp_mat.setFromTriplets(dgdp.begin(), dgdp.end());
    return dgdp_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxPlus(SystemState state, double time, double timestep) const
{
    std::vector<T> dhdx;
    dhdx.reserve(67);
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

Eigen::SparseMatrix<double> DroneTrajectory::dhdxCurr(SystemState prev, double timestep) const
{
    std::vector<T> dhdx;
    dhdx.reserve(30);
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

Eigen::SparseMatrix<double> DroneTrajectory::dhdzPlus(SystemState state, double timestep) const
{
    std::vector<T> dhdz;
    dhdz.reserve(75);
    
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

Eigen::SparseMatrix<double> DroneTrajectory::dhdzCurr() const
{
    std::vector<T> dhdz;
    dhdz.reserve(26);
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

Eigen::SparseMatrix<double> DroneTrajectory::dhdy(double timestep) const
{
    std::vector<T> dhdy;
    dhdy.reserve(4);
   
    dhdy.push_back(T(eiphi, 0, timestep));
    dhdy.push_back(T(desRollRate, 0, m_ctrlParams.at(roll).kp));
    dhdy.push_back(T(eitheta, 1, timestep));
    dhdy.push_back(T(desPitchRate, 1, m_ctrlParams.at(pitch).kp ));
    
    Eigen::SparseMatrix<double> dhdy_mat(NUM_Z_STATES, NUM_Y_STATES);
    dhdy_mat.setFromTriplets(dhdy.begin(), dhdy.end()); 
    return dhdy_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdp(SystemState state, double time)
{
    std::vector<T> dhdp;
    dhdp.reserve(38);

    dhdp.push_back(T(desVelX, kppx, (m_ref.at(refx)(time) - state.plant(x))*std::cos(state.plant(psi)) 
                                  + (m_ref.at(refy)(time) - state.plant(y))*std::sin(state.plant(psi)) ));
    dhdp.push_back(T(desVelX, kipx, state.alge(eix)));
    dhdp.push_back(T(desVelX, kdpx, state.alge(edx)));

    dhdp.push_back(T(desVelY, kppy, - (m_ref.at(refx)(time) - state.plant(x))*std::sin(state.plant(psi)) 
                                  + (m_ref.at(refy)(time) - state.plant(y))*std::cos(state.plant(psi)) ));
    dhdp.push_back(T(desVelY, kipy, state.alge(eiy)));
    dhdp.push_back(T(desVelY, kdpy, state.alge(edy)));

    dhdp.push_back(T(desVelZ, kppz, m_ref.at(refz)(time) - state.plant(z)));
    dhdp.push_back(T(desVelZ, kipz, state.alge(eiz)));
    dhdp.push_back(T(desVelZ, kdpz, state.alge(edz)));

    dhdp.push_back(T(desThrust, kpvz, m_thrustScale*(state.alge(desVelZ) - state.plant(zdot)) ));
    dhdp.push_back(T(desThrust, kivz, m_thrustScale*state.alge(eizdot) ));
    dhdp.push_back(T(desThrust, kdvz, m_thrustScale*state.alge(edzdot) ));

    dhdp.push_back(T(desRollRate, kpphi, state.alge(desRoll) - 180/M_PI*state.plant(phi) ));
    dhdp.push_back(T(desRollRate, kiphi, state.alge(eiphi) ));
    dhdp.push_back(T(desRollRate, kdphi, state.alge(edphi) ));

    dhdp.push_back(T(desPitchRate, kptheta, state.alge(desPitch) - 180/M_PI*state.plant(theta) ));
    dhdp.push_back(T(desPitchRate, kitheta, state.alge(eitheta) ));
    dhdp.push_back(T(desPitchRate, kdtheta, state.alge(edtheta) ));

    dhdp.push_back(T(desYawRate, kppsi, m_ref.at(refyaw)(time) - 180/M_PI*state.plant(psi) ));
    dhdp.push_back(T(desYawRate, kipsi, state.alge(eipsi) ));
    dhdp.push_back(T(desYawRate, kdpsi, state.alge(edpsi) ));

    dhdp.push_back(T(desRollOutput, kpp, state.alge(desRollRate) - 180/M_PI*state.plant(p) ));
    dhdp.push_back(T(desRollOutput, kip, state.alge(eip) ));
    dhdp.push_back(T(desRollOutput, kdp, state.alge(edp) ));

    dhdp.push_back(T(desPitchOutput, kpq, state.alge(desPitchRate) - 180/M_PI*state.plant(q) ));
    dhdp.push_back(T(desPitchOutput, kiq, state.alge(eiq) ));
    dhdp.push_back(T(desPitchOutput, kdq, state.alge(edq) ));

    dhdp.push_back(T(desYawOutput, kpr, state.alge(desYawRate) - 180/M_PI*state.plant(r) ));
    dhdp.push_back(T(desYawOutput, kir, state.alge(eir) ));
    dhdp.push_back(T(desYawOutput, kdr, state.alge(edr) ));
    
    Eigen::SparseMatrix<double> dhdp_mat(NUM_Z_STATES, NUM_PARAMETERS);
    dhdp_mat.setFromTriplets(dhdp.begin(), dhdp.end()); 
    return dhdp_mat;

}

Eigen::SparseMatrix<double> DroneTrajectory::d2fdx2(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    std::vector<T> d2fdx2;
    d2fdx2.reserve(59);

    d2fdx2.emplace_back(T(39, 3,  - plant(r)*cos(plant(phi))*tan(plant(theta)) - plant(q)*sin(plant(phi))*tan(plant(theta)) ));
    d2fdx2.emplace_back(T(39, 4,  plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(39, 10,  cos(plant(phi))*tan(plant(theta))));
    d2fdx2.emplace_back(T(39, 11,  -sin(plant(phi))*tan(plant(theta))));
    d2fdx2.emplace_back(T(40, 3,  plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(40, 4,  2*plant(r)*cos(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1) + 2*plant(q)*sin(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(40, 10,  sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(40, 11,  cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(46, 3,  cos(plant(phi))*tan(plant(theta))));
    d2fdx2.emplace_back(T(46, 4,  sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(47, 3,  -sin(plant(phi))*tan(plant(theta))));
    d2fdx2.emplace_back(T(47, 4,  cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1)));
    d2fdx2.emplace_back(T(51, 3,  plant(r)*sin(plant(phi)) - plant(q)*cos(plant(phi))));
    d2fdx2.emplace_back(T(51, 10,  -sin(plant(phi))));
    d2fdx2.emplace_back(T(51, 11,  -cos(plant(phi))));
    d2fdx2.emplace_back(T(58, 3,  -sin(plant(phi))));
    d2fdx2.emplace_back(T(59, 3,  -cos(plant(phi))));
    d2fdx2.emplace_back(T(63, 3,  - (plant(r)*cos(plant(phi)))/cos(plant(theta)) - (plant(q)*sin(plant(phi)))/cos(plant(theta))));
    d2fdx2.emplace_back(T(63, 4,  (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(63, 10,  cos(plant(phi))/cos(plant(theta))));
    d2fdx2.emplace_back(T(63, 11,  -sin(plant(phi))/cos(plant(theta))));
    d2fdx2.emplace_back(T(64, 3,  (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(64, 4,  (plant(r)*cos(plant(phi)))/cos(plant(theta)) + (plant(q)*sin(plant(phi)))/cos(plant(theta)) + std::pow((2*plant(r)*cos(plant(phi))*std::pow(sin(plant(theta)), 2))/cos(plant(theta)), 3) + std::pow((2*plant(q)*sin(plant(phi))*std::pow(sin(plant(theta)), 2))/cos(plant(theta)), 3)));
    d2fdx2.emplace_back(T(64, 10,  (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(64, 11,  (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(70, 3,  cos(plant(phi))/cos(plant(theta))));
    d2fdx2.emplace_back(T(70, 4,  (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(71, 3,  -sin(plant(phi))/cos(plant(theta))));
    d2fdx2.emplace_back(T(71, 4,  (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2)));
    d2fdx2.emplace_back(T(75, 3,  -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(75, 4,  -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(75, 5,  (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(76, 3,  -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(76, 4,  -(alge(ft)*cos(plant(phi))*cos(plant(psi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(76, 5,  -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(77, 3,  (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(77, 4,  -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(77, 5,  -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(87, 3,  (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(87, 4,  -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(87, 5,  (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(88, 3,  -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(88, 4,  -(alge(ft)*cos(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(88, 5,  (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(89, 3,  (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(89, 4,  (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(89, 5,  (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass));
    d2fdx2.emplace_back(T(99, 3,  -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(99, 4,  (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(100, 3,  (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(100, 4,  -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass));
    d2fdx2.emplace_back(T(118, 11,  (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix));
    d2fdx2.emplace_back(T(119, 10,  (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix));
    d2fdx2.emplace_back(T(129, 11,  -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy));
    d2fdx2.emplace_back(T(131, 9,  -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy));
    d2fdx2.emplace_back(T(141, 10,  (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz));
    d2fdx2.emplace_back(T(142, 9,  (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz));

    Eigen::SparseMatrix<double> d2fdx2_mat(NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2fdx2_mat.setFromTriplets(d2fdx2.begin(), d2fdx2.end());
    return d2fdx2_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2fdxdz(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    std::vector<T> d2fdxdz;
    d2fdxdz.reserve(8);

    d2fdxdz.emplace_back(T(6*NUM_PLANT_STATES + 3, ft, ((std::cos(plant(phi))*std::sin(plant(psi)) - std::cos(plant(psi))*std::sin(plant(phi))*std::sin(plant(theta))))/m_droneParams.mass));
    d2fdxdz.emplace_back(T(6*NUM_PLANT_STATES + 4, ft,  ( std::cos(plant(phi))*std::cos(plant(psi))*std::cos(plant(theta)))/m_droneParams.mass));
    d2fdxdz.emplace_back(T(6*NUM_PLANT_STATES + 5, ft, ((std::cos(plant(psi))*std::sin(plant(phi)) - std::cos(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass));

    d2fdxdz.emplace_back(T(7*NUM_PLANT_STATES + 3, ft, -((std::cos(plant(phi))*std::cos(plant(psi)) + std::sin(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass));
    d2fdxdz.emplace_back(T(7*NUM_PLANT_STATES + 4, ft, (std::cos(plant(phi))*std::cos(plant(theta))*std::sin(plant(psi)))/m_droneParams.mass));
    d2fdxdz.emplace_back(T(7*NUM_PLANT_STATES + 5, ft, ((std::sin(plant(phi))*std::sin(plant(psi)) + std::cos(plant(phi))*std::cos(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass));

    d2fdxdz.emplace_back(T(8*NUM_PLANT_STATES + 3, ft, -(std::cos(plant(theta))*std::sin(plant(phi)))/m_droneParams.mass));
    d2fdxdz.emplace_back(T(8*NUM_PLANT_STATES + 4, ft, -(std::cos(plant(phi))*std::sin(plant(theta)))/m_droneParams.mass));

    Eigen::SparseMatrix<double> d2fdxdz_mat(NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_Z_STATES);
    d2fdxdz_mat.setFromTriplets(d2fdxdz.begin(), d2fdxdz.end());
    return d2fdxdz_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2fdzdx(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    std::vector<T> d2fdzdx;
    d2fdzdx.reserve(8);

    d2fdzdx.emplace_back(T(xdot*NUM_PLANT_STATES + ft, phi, (cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdzdx.emplace_back(T(xdot*NUM_PLANT_STATES + ft, theta, (cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass));
    d2fdzdx.emplace_back(T(xdot*NUM_PLANT_STATES + ft, psi, (cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass));

    d2fdzdx.emplace_back(T(ydot*NUM_PLANT_STATES + ft, phi, -(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass));
    d2fdzdx.emplace_back(T(ydot*NUM_PLANT_STATES + ft, theta, (cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass));
    d2fdzdx.emplace_back(T(ydot*NUM_PLANT_STATES + ft, psi, (sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta)))/m_droneParams.mass));

    d2fdzdx.emplace_back(T(zdot*NUM_PLANT_STATES + ft, phi, -(cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass));
    d2fdzdx.emplace_back(T(zdot*NUM_PLANT_STATES + ft, theta, -(cos(plant(phi))*sin(plant(theta)))/m_droneParams.mass));
   
    Eigen::SparseMatrix<double> d2fdzdx_mat(NUM_PLANT_STATES*NUM_Z_STATES, NUM_PLANT_STATES);
    d2fdzdx_mat.setFromTriplets(d2fdzdx.begin(), d2fdzdx.end());
    return d2fdzdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2gdx2_mult_dxdwo(SystemState state, Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo)
{
    std::vector<T> d2gdx2;
    d2gdx2.reserve(3*NUM_STATES*2);
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        for (int i = 0; i < NUM_STATES; i++)
        {
            d2gdx2.emplace_back(T(xdot, i, -m_ctrlParams.at(velY).kp*cos(state.plant(psi))*dxdwo(psi, i) ));
            d2gdx2.emplace_back(T(ydot, i, -m_ctrlParams.at(velY).kp*sin(state.plant(psi))*dxdwo(psi, i) ));
            d2gdx2.emplace_back(T(psi, i, -m_ctrlParams.at(velY).kp*cos(state.plant(psi))*dxdwo(xdot, i) ));
            d2gdx2.emplace_back(T(psi, i, -m_ctrlParams.at(velY).kp*sin(state.plant(psi))*dxdwo(ydot, i) ));
            d2gdx2.emplace_back(T(psi, i,  -m_ctrlParams.at(velY).kp*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))*dxdwo(psi, i) ));
        }
        
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        for(int i = 0; i < NUM_STATES; i++)
        {
            d2gdx2.emplace_back(T(NUM_PLANT_STATES + xdot, i,  m_ctrlParams.at(velX).kp*sin(state.plant(psi))*dxdwo(psi, i) ));
            d2gdx2.emplace_back(T(NUM_PLANT_STATES + ydot, i, -m_ctrlParams.at(velX).kp*cos(state.plant(psi))*dxdwo(psi, i) ));
            d2gdx2.emplace_back(T(NUM_PLANT_STATES + psi, i,   m_ctrlParams.at(velX).kp*sin(state.plant(psi))*dxdwo(xdot, i) ));
            d2gdx2.emplace_back(T(NUM_PLANT_STATES + psi, i,  -m_ctrlParams.at(velX).kp*cos(state.plant(psi))*dxdwo(ydot, i) ));
            d2gdx2.emplace_back(T(NUM_PLANT_STATES + psi, i,   m_ctrlParams.at(velX).kp*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))*dxdwo(psi, i) ));
        }
    }

    Eigen::SparseMatrix<double> d2gdx2_mat(NUM_Y_STATES*NUM_PLANT_STATES, NUM_STATES);
    d2gdx2_mat.setFromTriplets(d2gdx2.begin(), d2gdx2.end());
    return d2gdx2_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2hdx2_plus_mult_dxdwo(SystemState state, double time, double timestep, Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(10*NUM_STATES);

    for(int i = 0; i < NUM_STATES; i ++)
    {
        d2hdx2_plus.emplace_back(T(eix*NUM_PLANT_STATES+x, i, timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eix*NUM_PLANT_STATES+y, i, -timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eix*NUM_PLANT_STATES+psi, i, timestep*sin(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(eix*NUM_PLANT_STATES+psi, i, -timestep*cos(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(eix*NUM_PLANT_STATES+psi, i, -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(edx*NUM_PLANT_STATES+x, i, 1/timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edx*NUM_PLANT_STATES+y, i, -1/timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edx*NUM_PLANT_STATES+psi, i, 1/timestep*sin(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(edx*NUM_PLANT_STATES+psi, i, -1/timestep*cos(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(edx*NUM_PLANT_STATES+psi, i, (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(desVelX*NUM_PLANT_STATES+x, i, m_ctrlParams.at(posX).kp*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(desVelX*NUM_PLANT_STATES+y, i, -m_ctrlParams.at(posX).kp*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(desVelX*NUM_PLANT_STATES+psi, i, m_ctrlParams.at(posX).kp*sin(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(desVelX*NUM_PLANT_STATES+psi, i, -m_ctrlParams.at(posX).kp*cos(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(desVelX*NUM_PLANT_STATES+psi, i, -m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(eiy*NUM_PLANT_STATES+x, i, timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eiy*NUM_PLANT_STATES+y, i, timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eiy*NUM_PLANT_STATES+psi, i, timestep*cos(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(eiy*NUM_PLANT_STATES+psi, i, timestep*sin(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(eiy*NUM_PLANT_STATES+psi, i, -timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -x))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(edy*NUM_PLANT_STATES+x, i, 1/timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edy*NUM_PLANT_STATES+y, i, 1/timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edy*NUM_PLANT_STATES+psi, i, 1/timestep*cos(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(edy*NUM_PLANT_STATES+psi, i, 1/timestep*sin(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(edy*NUM_PLANT_STATES+psi, i, (state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(desVelY*NUM_PLANT_STATES+x, i, m_ctrlParams.at(posY).kp*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(desVelY*NUM_PLANT_STATES+y, i, m_ctrlParams.at(posY).kp*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(desVelY*NUM_PLANT_STATES+psi, i, m_ctrlParams.at(posY).kp*cos(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(desVelY*NUM_PLANT_STATES+psi, i, m_ctrlParams.at(posY).kp*sin(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(desVelY*NUM_PLANT_STATES+psi, i, -m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(eixdot*NUM_PLANT_STATES+xdot, i, timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eixdot*NUM_PLANT_STATES+ydot, i, -timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eixdot*NUM_PLANT_STATES+psi, i, timestep*sin(state.plant(psi))*dxdwo(x, i) ));
        d2hdx2_plus.emplace_back(T(eixdot*NUM_PLANT_STATES+psi, i, -timestep*cos(state.plant(psi))*dxdwo(y, i) ));
        d2hdx2_plus.emplace_back(T(eixdot*NUM_PLANT_STATES+psi, i, timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(edxdot*NUM_PLANT_STATES+xdot, i, 1/timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edxdot*NUM_PLANT_STATES+ydot, i, -1/timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edxdot*NUM_PLANT_STATES+psi, i, 1/timestep*sin(state.plant(psi))*dxdwo(xdot, i) ));
        d2hdx2_plus.emplace_back(T(edxdot*NUM_PLANT_STATES+psi, i, -1/timestep*cos(state.plant(psi))*dxdwo(ydot, i) ));
        d2hdx2_plus.emplace_back(T(edxdot*NUM_PLANT_STATES+psi, i, (state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(eiydot*NUM_PLANT_STATES+xdot, i, timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eiydot*NUM_PLANT_STATES+ydot, i, timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(eiydot*NUM_PLANT_STATES+psi, i, timestep*cos(state.plant(psi))*dxdwo(xdot, i) ));
        d2hdx2_plus.emplace_back(T(eiydot*NUM_PLANT_STATES+psi, i, timestep*sin(state.plant(psi))*dxdwo(ydot, i) ));
        d2hdx2_plus.emplace_back(T(eiydot*NUM_PLANT_STATES+psi, i, timestep*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))*dxdwo(psi, i) ));

        d2hdx2_plus.emplace_back(T(edydot*NUM_PLANT_STATES+xdot, i, 1/timestep*cos(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edydot*NUM_PLANT_STATES+ydot, i, 1/timestep*sin(state.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_plus.emplace_back(T(edydot*NUM_PLANT_STATES+psi, i, 1/timestep*cos(state.plant(psi))*dxdwo(xdot, i) ));
        d2hdx2_plus.emplace_back(T(edydot*NUM_PLANT_STATES+psi, i, 1/timestep*sin(state.plant(psi))*dxdwo(ydot, i) ));
        d2hdx2_plus.emplace_back(T(edydot*NUM_PLANT_STATES+psi, i, (state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))/timestep*dxdwo(psi, i) ));
    }
  
    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_Z_STATES*NUM_PLANT_STATES, NUM_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2hdx2_curr_mult_dxdwo(SystemState prev, double timestep, Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo)
{
    std::vector<T> d2hdx2_curr;
    d2hdx2_curr.reserve(4*NUM_STATES);
    
    for (int i = 0; i < NUM_STATES; i++)
    {
        d2hdx2_curr.emplace_back(T(edx*NUM_PLANT_STATES + x, i, -1/timestep*sin(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edx*NUM_PLANT_STATES + y, i, 1/timestep*cos(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edx*NUM_PLANT_STATES + psi, i, -sin(prev.plant(psi))/timestep*dxdwo(x, i) ));
        d2hdx2_curr.emplace_back(T(edx*NUM_PLANT_STATES + psi, i, cos(prev.plant(psi))/timestep*dxdwo(y, i) ));
        d2hdx2_curr.emplace_back(T(edx*NUM_PLANT_STATES + psi, i, -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_curr.emplace_back(T(edy*NUM_PLANT_STATES + x, i, -1/timestep*cos(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edy*NUM_PLANT_STATES + y, i, -1/timestep*sin(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edy*NUM_PLANT_STATES + psi, i, -cos(prev.plant(psi))/timestep*dxdwo(x, i) ));
        d2hdx2_curr.emplace_back(T(edy*NUM_PLANT_STATES + psi, i, -sin(prev.plant(psi))/timestep*dxdwo(y, i) ));
        d2hdx2_curr.emplace_back(T(edy*NUM_PLANT_STATES + psi, i, -(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_curr.emplace_back(T(edxdot*NUM_PLANT_STATES + xdot, i, -1/timestep*sin(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edxdot*NUM_PLANT_STATES + ydot, i, 1/timestep*cos(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edxdot*NUM_PLANT_STATES + psi, i, -sin(prev.plant(psi))/timestep*dxdwo(xdot, i) ));
        d2hdx2_curr.emplace_back(T(edxdot*NUM_PLANT_STATES + psi, i, cos(prev.plant(psi))/timestep*dxdwo(ydot, i) ));
        d2hdx2_curr.emplace_back(T(edxdot*NUM_PLANT_STATES + psi, i, -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep*dxdwo(psi, i) ));

        d2hdx2_curr.emplace_back(T(edydot*NUM_PLANT_STATES + xdot, i, -1/timestep*cos(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edydot*NUM_PLANT_STATES + ydot, i, -1/timestep*sin(prev.plant(psi))*dxdwo(psi, i) ));
        d2hdx2_curr.emplace_back(T(edydot*NUM_PLANT_STATES + psi, i, -cos(prev.plant(psi))/timestep*dxdwo(xdot, i) ));
        d2hdx2_curr.emplace_back(T(edydot*NUM_PLANT_STATES + psi, i, -sin(prev.plant(psi))/timestep*dxdwo(ydot, i) ));
        d2hdx2_curr.emplace_back(T(edydot*NUM_PLANT_STATES + psi, i, -(prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep*dxdwo(psi, i) ));
    }
    Eigen::SparseMatrix<double> d2hdx2_curr_mat(NUM_Z_STATES*NUM_PLANT_STATES, NUM_STATES);
    d2hdx2_curr_mat.setFromTriplets(d2hdx2_curr.begin(), d2hdx2_curr.end());
    return d2hdx2_curr_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2hdz2_plus_mult_dzdwo(Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo)
{
    std::vector<T> d2hdz2_plus;
    d2hdz2_plus.reserve(4*NUM_STATES*4);
    for(int i = 0; i < NUM_STATES; i++) {
        d2hdz2_plus.emplace_back(T(ft*NUM_Z_STATES + w1, i, m_droneParams.kf * 2 * dzdwo(w1, i) ));
        d2hdz2_plus.emplace_back(T(ft*NUM_Z_STATES + w2, i, m_droneParams.kf * 2 * dzdwo(w2, i) ));
        d2hdz2_plus.emplace_back(T(ft*NUM_Z_STATES + w3, i, m_droneParams.kf * 2 * dzdwo(w3, i) ));
        d2hdz2_plus.emplace_back(T(ft*NUM_Z_STATES + w4, i, m_droneParams.kf * 2 * dzdwo(w4, i) )) ;

        d2hdz2_plus.emplace_back(T(tx*NUM_Z_STATES + w1, i, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w1,i) ));
        d2hdz2_plus.emplace_back(T(tx*NUM_Z_STATES + w2, i, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w2,i) ));
        d2hdz2_plus.emplace_back(T(tx*NUM_Z_STATES + w3, i,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w3,i) ));
        d2hdz2_plus.emplace_back(T(tx*NUM_Z_STATES + w4, i,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w4,i) ));

        d2hdz2_plus.emplace_back(T(ty*NUM_Z_STATES + w1, i,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w1, i) ));
        d2hdz2_plus.emplace_back(T(ty*NUM_Z_STATES + w2, i, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w2, i) ));
        d2hdz2_plus.emplace_back(T(ty*NUM_Z_STATES + w3, i, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w3, i) ));
        d2hdz2_plus.emplace_back(T(ty*NUM_Z_STATES + w4, i,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) * dzdwo(w4, i) ));

        d2hdz2_plus.emplace_back(T(tz*NUM_Z_STATES + w1, i,    m_droneParams.km * 2 * dzdwo(w1, i) )) ;
        d2hdz2_plus.emplace_back(T(tz*NUM_Z_STATES + w2, i,  - m_droneParams.km * 2 * dzdwo(w2, i) )) ;
        d2hdz2_plus.emplace_back(T(tz*NUM_Z_STATES + w3, i,    m_droneParams.km * 2 * dzdwo(w3, i) )) ;
        d2hdz2_plus.emplace_back(T(tz*NUM_Z_STATES + w4, i,  - m_droneParams.km * 2 * dzdwo(w4, i) )) ;
    }

    Eigen::SparseMatrix<double> d2hdz2_plus_mat(NUM_Z_STATES*NUM_Z_STATES, NUM_STATES);
    d2hdz2_plus_mat.setFromTriplets(d2hdz2_plus.begin(), d2hdz2_plus.end());
    return d2hdz2_plus_mat;
}


Eigen::Tensor<double, 3> DroneTrajectory::dfdz_mult_d2zdwo2(SystemState state, Eigen::Tensor<double, 3> const& d2zdwo2)
{
    Eigen::Tensor<double, 3> result(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    result.setZero();
   
    result.chip(xdot, 0) += (1.0/m_droneParams.mass * (std::sin(state.plant(phi)) * std::sin(state.plant(psi)) 
                             + std::cos(state.plant(phi))*std::cos(state.plant(psi))*std::sin(state.plant(theta)))) * d2zdwo2.chip(ft, 0);
    result.chip(ydot, 0) += (1.0/m_droneParams.mass * (std::cos(state.plant(phi))*std::sin(state.plant(psi))*std::sin(state.plant(theta)) 
                             - std::cos(state.plant(psi))*std::sin(state.plant(phi)))) * d2zdwo2.chip(ft, 0);
    result.chip(zdot, 0) += (1.0/m_droneParams.mass * (std::cos(state.plant(phi)) * std::cos(state.plant(theta)))) * d2zdwo2.chip(ft, 0);

    result.chip(p, 0) += 1.0/m_droneParams.Ix * d2zdwo2.chip(tx, 0);
    result.chip(q, 0) += 1.0/m_droneParams.Iy * d2zdwo2.chip(ty, 0);
    result.chip(r, 0) += 1.0/m_droneParams.Iz * d2zdwo2.chip(tz, 0);

    return result;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2edx_dx_curr2(SystemState prev, double timestep)
{
    std::vector<T> d2hdx2_curr;
    d2hdx2_curr.reserve(5);
    
    d2hdx2_curr.emplace_back(T(x, psi, -1/timestep*sin(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(y, psi, 1/timestep*cos(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(psi, x, -sin(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, y, cos(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, psi, -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_curr_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_curr_mat.setFromTriplets(d2hdx2_curr.begin(), d2hdx2_curr.end());
    return d2hdx2_curr_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2edy_dx_curr2(SystemState prev, double timestep)
{
    std::vector<T> d2hdx2_curr;
    d2hdx2_curr.reserve(5);
    
    d2hdx2_curr.emplace_back(T(x, psi, -1/timestep*cos(prev.plant(psi))));
    d2hdx2_curr.emplace_back(T(y, psi, -1/timestep*sin(prev.plant(psi))));
    d2hdx2_curr.emplace_back(T(psi, x, -cos(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, y, -sin(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, psi, -(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_curr_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_curr_mat.setFromTriplets(d2hdx2_curr.begin(), d2hdx2_curr.end());
    return d2hdx2_curr_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edxdot_dx_curr2(SystemState prev, double timestep)
{
    std::vector<T> d2hdx2_curr;
    d2hdx2_curr.reserve(5);
    
    d2hdx2_curr.emplace_back(T(xdot, psi, -1/timestep*sin(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(ydot, psi, 1/timestep*cos(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(psi, xdot, -sin(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, ydot, cos(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, psi, -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_curr_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_curr_mat.setFromTriplets(d2hdx2_curr.begin(), d2hdx2_curr.end());
    return d2hdx2_curr_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edydot_dx_curr2(SystemState prev, double timestep)
{
    std::vector<T> d2hdx2_curr;
    d2hdx2_curr.reserve(5);
    
    d2hdx2_curr.emplace_back(T(xdot, psi, -1/timestep*cos(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(ydot, psi, -1/timestep*sin(prev.plant(psi)) ));
    d2hdx2_curr.emplace_back(T(psi, xdot, -cos(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, ydot, -sin(prev.plant(psi))/timestep ));
    d2hdx2_curr.emplace_back(T(psi, psi, -(prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_curr_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_curr_mat.setFromTriplets(d2hdx2_curr.begin(), d2hdx2_curr.end());
    return d2hdx2_curr_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2eix_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);
    
    d2hdx2_plus.emplace_back(T(x, psi, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, -timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, -timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2eiy_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -x)) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edx_dx_plus2(SystemState state, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, -1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, -1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edy_dx_plus2(SystemState state, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, 1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, 1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, (state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2eixdot_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(xdot, psi, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(ydot, psi, -timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, xdot, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, ydot, -timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2eiydot_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(xdot, psi, timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(ydot, psi, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, xdot, timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, ydot, timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, timestep*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edxdot_dx_plus2(SystemState state, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(xdot, psi, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(ydot, psi, -1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, xdot, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, ydot, -1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, (state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2edydot_dx_plus2(SystemState state, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(xdot, psi, 1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(ydot, psi, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, xdot, 1/timestep*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, ydot, 1/timestep*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, (state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))/timestep ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2desVelx_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, m_ctrlParams.at(posX).kp*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, -m_ctrlParams.at(posX).kp*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, m_ctrlParams.at(posX).kp*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, -m_ctrlParams.at(posX).kp*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2desVely_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, m_ctrlParams.at(posY).kp*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, m_ctrlParams.at(posY).kp*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, m_ctrlParams.at(posY).kp*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, m_ctrlParams.at(posY).kp*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2ft_dz_plus2(SystemState state)
{
    std::vector<T> d2hdz2_plus;
    d2hdz2_plus.reserve(4);

    d2hdz2_plus.emplace_back(T(w1, w1, m_droneParams.kf * 2  ));
    d2hdz2_plus.emplace_back(T(w2, w2, m_droneParams.kf * 2  ));
    d2hdz2_plus.emplace_back(T(w3, w3, m_droneParams.kf * 2  ));
    d2hdz2_plus.emplace_back(T(w4, w4, m_droneParams.kf * 2  )) ;

    Eigen::SparseMatrix<double> d2hdz2_plus_mat(NUM_Z_STATES, NUM_Z_STATES);
    d2hdz2_plus_mat.setFromTriplets(d2hdz2_plus.begin(), d2hdz2_plus.end());
    return d2hdz2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2tx_dz_plus2(SystemState state)
{
    std::vector<T> d2hdz2_plus;
    d2hdz2_plus.reserve(4);

    d2hdz2_plus.emplace_back(T(w1, w1, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));
    d2hdz2_plus.emplace_back(T(w2, w2, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ));
    d2hdz2_plus.emplace_back(T(w3, w3,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));
    d2hdz2_plus.emplace_back(T(w4, w4,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));

    Eigen::SparseMatrix<double> d2hdz2_plus_mat(NUM_Z_STATES, NUM_Z_STATES);
    d2hdz2_plus_mat.setFromTriplets(d2hdz2_plus.begin(), d2hdz2_plus.end());
    return d2hdz2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2ty_dz_plus2(SystemState state)
{
    std::vector<T> d2hdz2_plus;
    d2hdz2_plus.reserve(4);

    d2hdz2_plus.emplace_back(T(w1, w1,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));
    d2hdz2_plus.emplace_back(T(w2, w2, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));
    d2hdz2_plus.emplace_back(T(w3, w3, - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));
    d2hdz2_plus.emplace_back(T(w4, w4,   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2)  ));

    Eigen::SparseMatrix<double> d2hdz2_plus_mat(NUM_Z_STATES, NUM_Z_STATES);
    d2hdz2_plus_mat.setFromTriplets(d2hdz2_plus.begin(), d2hdz2_plus.end());
    return d2hdz2_plus_mat;        
}
Eigen::SparseMatrix<double> DroneTrajectory::d2tz_dz_plus2(SystemState state)
{
    std::vector<T> d2hdz2_plus;
    d2hdz2_plus.reserve(4);

    d2hdz2_plus.emplace_back(T(w1, w1,    m_droneParams.km * 2   )) ;
    d2hdz2_plus.emplace_back(T(w2, w2,  - m_droneParams.km * 2   )) ;
    d2hdz2_plus.emplace_back(T(w3, w3,    m_droneParams.km * 2  )) ;
    d2hdz2_plus.emplace_back(T(w4, w4,  - m_droneParams.km * 2  )) ;

    Eigen::SparseMatrix<double> d2hdz2_plus_mat(NUM_Z_STATES, NUM_Z_STATES);
    d2hdz2_plus_mat.setFromTriplets(d2hdz2_plus.begin(), d2hdz2_plus.end());
    return d2hdz2_plus_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2desRolldx2(SystemState state) const
{
    std::vector<T> d2gdx2;
    d2gdx2.reserve(6);

    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        d2gdx2.push_back(T(xdot, psi, -m_ctrlParams.at(velY).kp*cos(state.plant(psi))));
        d2gdx2.push_back(T(ydot, psi, -m_ctrlParams.at(velY).kp*sin(state.plant(psi))));
        d2gdx2.push_back(T(psi, xdot, -m_ctrlParams.at(velY).kp*cos(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, ydot, -m_ctrlParams.at(velY).kp*sin(state.plant(psi))));
        d2gdx2.push_back(T(psi, psi, -m_ctrlParams.at(velY).kp*(-state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi)))));
    }
    
    Eigen::SparseMatrix<double> d2gdx2_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2gdx2_mat.setFromTriplets(d2gdx2.begin(), d2gdx2.end());
    return d2gdx2_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::d2desPitchdx2(SystemState state) const
{
    std::vector<T> d2gdx2;
    d2gdx2.reserve(6);
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        d2gdx2.push_back(T(xdot, psi, m_ctrlParams.at(velX).kp*sin(state.plant(psi))));
        d2gdx2.push_back(T(ydot, psi, -m_ctrlParams.at(velX).kp*cos(state.plant(psi))));
        d2gdx2.push_back(T(psi, xdot, m_ctrlParams.at(velX).kp*sin(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, ydot, -m_ctrlParams.at(velX).kp*cos(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, psi, -m_ctrlParams.at(velX).kp*(-state.plant(xdot)*cos(state.plant(psi)) - state.plant(ydot)*sin(state.plant(psi)))));
    }
    Eigen::SparseMatrix<double> d2gdx2_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2gdx2_mat.setFromTriplets(d2gdx2.begin(), d2gdx2.end());
    return d2gdx2_mat;
}