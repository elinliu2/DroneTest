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
        dgdx.push_back(T(0, xdot, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*sin(state.plant(psi))));
        dgdx.push_back(T(0, ydot, m_ctrlParams.at(velY).kp*m_sf(kpvy)*cos(state.plant(psi))));
        dgdx.push_back(T(0, psi, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdx.push_back(T(1, xdot, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*cos(state.plant(psi))));
        dgdx.push_back(T(1, ydot, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*sin(state.plant(psi))));
        dgdx.push_back(T(1, psi, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*(-state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi)))));
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
        dgdz.push_back(T(0, desVelY, -m_ctrlParams.at(velY).kp*m_sf(kpvy)));
        dgdz.push_back(T(0, eiydot, -m_ctrlParams.at(velY).ki*m_sf(kivy)));
        dgdz.push_back(T(0, edydot, -m_ctrlParams.at(velY).kd*m_sf(kdvy)));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdz.push_back(T(1, desVelX, m_ctrlParams.at(velX).kp*m_sf(kpvx)));
        dgdz.push_back(T(1, eixdot, m_ctrlParams.at(velX).ki*m_sf(kivx)));
        dgdz.push_back(T(1, edxdot, m_ctrlParams.at(velX).kd*m_sf(kdvx)));
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
        dgdp.push_back(T(0, kpvy, (-(state.alge(desVelY) + state.plant(xdot)*std::sin(state.plant(psi)) - state.plant(ydot)*std::cos(state.plant(psi)) ))*m_ctrlParams.at(velY).kp ));
        dgdp.push_back(T(0, kivy, -state.alge(eiydot)*m_ctrlParams.at(velY).ki));
        dgdp.push_back(T(0, kdvy, -state.alge(edydot)*m_ctrlParams.at(velY).kd));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdp.push_back(T(1, kpvx, (state.alge(desVelX) - state.plant(xdot)*std::cos(state.plant(psi)) - state.plant(ydot)*std::sin(state.plant(psi)) )*m_ctrlParams.at(velX).kp  ));
        dgdp.push_back(T(1, kivx, state.alge(eixdot)*m_ctrlParams.at(velX).ki));
        dgdp.push_back(T(1, kdvx, state.alge(edxdot)*m_ctrlParams.at(velX).kd));
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
    dhdx.push_back(T(desVelX, x, -m_ctrlParams.at(posX).kp*m_sf(kppx)*cos(state.plant(psi))));
    dhdx.push_back(T(desVelX, y, -m_ctrlParams.at(posX).kp*m_sf(kppx)*sin(state.plant(psi))));
    dhdx.push_back(T(desVelX, psi, m_ctrlParams.at(posX).kp*m_sf(kppx)*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) ));

    dhdx.push_back(T(eiy, x, timestep*sin(state.plant(psi))));
    dhdx.push_back(T(eiy, y, -timestep*cos(state.plant(psi))));
    dhdx.push_back(T(eiy, psi, -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)))));
    dhdx.push_back(T(edy, x, 1/timestep*sin(state.plant(psi))));
    dhdx.push_back(T(edy, y, -1/timestep*cos(state.plant(psi))));
    dhdx.push_back(T(edy, psi, (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, x, m_ctrlParams.at(posY).kp*m_sf(kppy)*sin(state.plant(psi))));
    dhdx.push_back(T(desVelY, y, -m_ctrlParams.at(posY).kp*m_sf(kppy)*cos(state.plant(psi))));
    dhdx.push_back(T(desVelY, psi, -m_ctrlParams.at(posY).kp*m_sf(kppy)*(cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) ));

    dhdx.push_back(T(eiz, z, -timestep));
    dhdx.push_back(T(edz, z, -1/timestep));
    dhdx.push_back(T(desVelZ, z, -m_ctrlParams.at(posZ).kp*m_sf(kppz)));

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

    dhdx.push_back(T(desThrust, zdot, -m_thrustScale*m_ctrlParams.at(velZ).kp*m_sf(kpvz) ));

    dhdx.push_back(T(eiphi, phi, -180/M_PI*timestep ));
    dhdx.push_back(T(edphi, phi, -180/M_PI/timestep ));
    dhdx.push_back(T(desRollRate, phi, -m_ctrlParams.at(roll).kp*m_sf(kpphi)*180/M_PI));

    dhdx.push_back(T(eitheta, theta, -180/M_PI*timestep ));
    dhdx.push_back(T(edtheta, theta, -180/M_PI/timestep ));
    dhdx.push_back(T(desPitchRate, theta, -m_ctrlParams.at(pitch).kp*m_sf(kptheta)*180/M_PI));

    dhdx.push_back(T(eipsi, psi, -180/M_PI*timestep ));
    dhdx.push_back(T(edpsi, psi, -180/M_PI/timestep ));
    dhdx.push_back(T(desYawRate, psi, -m_ctrlParams.at(yaw).kp*m_sf(kppsi)*180/M_PI));

    dhdx.push_back(T(delay_1_rollRate, p, -180/timestep/M_PI));
    dhdx.push_back(T(delay_1_pitchRate, q, -180/timestep/M_PI));

    dhdx.push_back(T(eip, p, -180/M_PI*timestep ));
    dhdx.push_back(T(edp, p, -180/M_PI/timestep*lpf_b0 ));
    dhdx.push_back(T(desRollOutput, p, -m_ctrlParams.at(rollRate).kp*m_sf(kpp)*180/M_PI));

    dhdx.push_back(T(eiq, q, -180/M_PI*timestep ));
    dhdx.push_back(T(edq, q, -180/M_PI/timestep*lpf_b0 ));
    dhdx.push_back(T(desPitchOutput, q, -m_ctrlParams.at(pitchRate).kp*m_sf(kpq)*180/M_PI));

    dhdx.push_back(T(eir, r, -180/M_PI*timestep ));
    dhdx.push_back(T(edr, r, -180/M_PI/timestep ));
    dhdx.push_back(T(desYawOutput, r, -m_ctrlParams.at(yawRate).kp*m_sf(kpr)*180/M_PI));

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
    
    dhdz.push_back(T(desVelX, eix, m_ctrlParams.at(posX).ki*m_sf(kipx)));
    dhdz.push_back(T(desVelX, edx, m_ctrlParams.at(posX).kd*m_sf(kdpx)));

    dhdz.push_back(T(desVelY, eiy, m_ctrlParams.at(posY).ki*m_sf(kipy)));
    dhdz.push_back(T(desVelY, edy, m_ctrlParams.at(posY).kd*m_sf(kdpy)));

    dhdz.push_back(T(desVelZ, eiz, m_ctrlParams.at(posZ).ki*m_sf(kipz)));
    dhdz.push_back(T(desVelZ, edz, m_ctrlParams.at(posZ).kd*m_sf(kdpz)));

    dhdz.push_back(T(eixdot, desVelX, timestep));
    dhdz.push_back(T(eiydot, desVelY, timestep));
    dhdz.push_back(T(eizdot, desVelZ, timestep));

    dhdz.push_back(T(desThrust, desVelZ, m_ctrlParams.at(velZ).kp*m_sf(kpvz)*m_thrustScale));
    dhdz.push_back(T(desThrust, eizdot, m_ctrlParams.at(velZ).ki*m_sf(kivz)*m_thrustScale));
    dhdz.push_back(T(desThrust, edzdot, m_ctrlParams.at(velZ).kd*m_sf(kdvz)*m_thrustScale));

    dhdz.push_back(T(desRollRate, eiphi, m_ctrlParams.at(roll).ki*m_sf(kiphi)));
    dhdz.push_back(T(desRollRate, edphi, m_ctrlParams.at(roll).kd*m_sf(kdphi)));

    dhdz.push_back(T(desPitchRate, eitheta, m_ctrlParams.at(pitch).ki*m_sf(kitheta)));
    dhdz.push_back(T(desPitchRate, edtheta, m_ctrlParams.at(pitch).kd*m_sf(kdtheta)));
    
    dhdz.push_back(T(desYawRate, eipsi, m_ctrlParams.at(yaw).ki*m_sf(kipsi)));
    dhdz.push_back(T(desYawRate, edpsi, m_ctrlParams.at(yaw).kd*m_sf(kdpsi)));
    
    dhdz.push_back(T(eip, desRollRate, timestep));
    dhdz.push_back(T(desRollOutput, desRollRate, m_ctrlParams.at(rollRate).kp*m_sf(kpp)));
    dhdz.push_back(T(desRollOutput, eip, m_ctrlParams.at(rollRate).ki*m_sf(kip)));
    dhdz.push_back(T(desRollOutput, edp, m_ctrlParams.at(rollRate).kd*m_sf(kdp)));

    dhdz.push_back(T(eiq, desPitchRate, timestep));
    dhdz.push_back(T(desPitchOutput, desPitchRate, m_ctrlParams.at(pitchRate).kp*m_sf(kpq)));
    dhdz.push_back(T(desPitchOutput, eiq, m_ctrlParams.at(pitchRate).ki*m_sf(kiq)));
    dhdz.push_back(T(desPitchOutput, edq, m_ctrlParams.at(pitchRate).kd*m_sf(kdq)));

    dhdz.push_back(T(eir, desYawRate, timestep));
    dhdz.push_back(T(desYawOutput, desYawRate, m_ctrlParams.at(yawRate).kp*m_sf(kpr)));
    dhdz.push_back(T(desYawOutput, eir, m_ctrlParams.at(yawRate).ki*m_sf(kir)));
    dhdz.push_back(T(desYawOutput, edr, m_ctrlParams.at(yawRate).kd*m_sf(kdr)));

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
    dhdy.push_back(T(desRollRate, 0, m_ctrlParams.at(roll).kp*m_sf(kpphi)));
    dhdy.push_back(T(eitheta, 1, timestep));
    dhdy.push_back(T(desPitchRate, 1, m_ctrlParams.at(pitch).kp*m_sf(kptheta)));
    
    Eigen::SparseMatrix<double> dhdy_mat(NUM_Z_STATES, NUM_Y_STATES);
    dhdy_mat.setFromTriplets(dhdy.begin(), dhdy.end()); 
    return dhdy_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdp(SystemState state, double time)
{
    std::vector<T> dhdp;
    dhdp.reserve(38);

    dhdp.push_back(T(desVelX, kppx, ((m_ref.at(refx)(time) - state.plant(x))*std::cos(state.plant(psi)) 
                                  + (m_ref.at(refy)(time) - state.plant(y))*std::sin(state.plant(psi)))*m_ctrlParams.at(posX).kp ));
    dhdp.push_back(T(desVelX, kipx, state.alge(eix)*m_ctrlParams.at(posX).ki));
    dhdp.push_back(T(desVelX, kdpx, state.alge(edx)*m_ctrlParams.at(posX).kd));

    dhdp.push_back(T(desVelY, kppy, (-(m_ref.at(refx)(time) - state.plant(x))*std::sin(state.plant(psi)) 
                                  + (m_ref.at(refy)(time) - state.plant(y))*std::cos(state.plant(psi)))*m_ctrlParams.at(posY).kp ));
    dhdp.push_back(T(desVelY, kipy, state.alge(eiy)*m_ctrlParams.at(posY).ki));
    dhdp.push_back(T(desVelY, kdpy, state.alge(edy)*m_ctrlParams.at(posY).kd));

    dhdp.push_back(T(desVelZ, kppz, (m_ref.at(refz)(time) - state.plant(z))*m_ctrlParams.at(posZ).kp));
    dhdp.push_back(T(desVelZ, kipz, state.alge(eiz)*m_ctrlParams.at(posZ).ki));
    dhdp.push_back(T(desVelZ, kdpz, state.alge(edz)*m_ctrlParams.at(posZ).kd));

    dhdp.push_back(T(desThrust, kpvz, m_thrustScale*(state.alge(desVelZ) - state.plant(zdot))*m_ctrlParams.at(velZ).kp ));
    dhdp.push_back(T(desThrust, kivz, m_thrustScale*state.alge(eizdot)*m_ctrlParams.at(velZ).ki ));
    dhdp.push_back(T(desThrust, kdvz, m_thrustScale*state.alge(edzdot)*m_ctrlParams.at(velZ).kd ));

    dhdp.push_back(T(desRollRate, kpphi, (state.alge(desRoll) - 180/M_PI*state.plant(phi))*m_ctrlParams.at(roll).kp ));
    dhdp.push_back(T(desRollRate, kiphi, state.alge(eiphi)*m_ctrlParams.at(roll).ki ));
    dhdp.push_back(T(desRollRate, kdphi, state.alge(edphi)*m_ctrlParams.at(roll).kd ));

    dhdp.push_back(T(desPitchRate, kptheta, (state.alge(desPitch) - 180/M_PI*state.plant(theta))*m_ctrlParams.at(pitch).kp ));
    dhdp.push_back(T(desPitchRate, kitheta, state.alge(eitheta)*m_ctrlParams.at(pitch).ki ));
    dhdp.push_back(T(desPitchRate, kdtheta, state.alge(edtheta)*m_ctrlParams.at(pitch).kd ));

    dhdp.push_back(T(desYawRate, kppsi, (m_ref.at(refyaw)(time) - 180/M_PI*state.plant(psi))*m_ctrlParams.at(yaw).kp ));
    dhdp.push_back(T(desYawRate, kipsi, state.alge(eipsi)*m_ctrlParams.at(yaw).ki ));
    dhdp.push_back(T(desYawRate, kdpsi, state.alge(edpsi)*m_ctrlParams.at(yaw).kd ));

    dhdp.push_back(T(desRollOutput, kpp, (state.alge(desRollRate) - 180/M_PI*state.plant(p))*m_ctrlParams.at(rollRate).kp ));
    dhdp.push_back(T(desRollOutput, kip, state.alge(eip)*m_ctrlParams.at(rollRate).ki ));
    dhdp.push_back(T(desRollOutput, kdp, state.alge(edp)*m_ctrlParams.at(rollRate).kd ));

    dhdp.push_back(T(desPitchOutput, kpq, (state.alge(desPitchRate) - 180/M_PI*state.plant(q))*m_ctrlParams.at(pitchRate).kp ));
    dhdp.push_back(T(desPitchOutput, kiq, state.alge(eiq)*m_ctrlParams.at(pitchRate).ki ));
    dhdp.push_back(T(desPitchOutput, kdq, state.alge(edq)*m_ctrlParams.at(pitchRate).kd ));

    dhdp.push_back(T(desYawOutput, kpr, (state.alge(desYawRate) - 180/M_PI*state.plant(r))*m_ctrlParams.at(yawRate).kp ));
    dhdp.push_back(T(desYawOutput, kir, state.alge(eir)*m_ctrlParams.at(yawRate).ki ));
    dhdp.push_back(T(desYawOutput, kdr, state.alge(edr)*m_ctrlParams.at(yawRate).kd ));
    
    Eigen::SparseMatrix<double> dhdp_mat(NUM_Z_STATES, NUM_PARAMETERS);
    dhdp_mat.setFromTriplets(dhdp.begin(), dhdp.end()); 
    return dhdp_mat;

}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2fdx2(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_tensor(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2fdx2_tensor.setZero();

    d2fdx2_tensor(3, 3, 3) = - plant(r)*cos(plant(phi))*tan(plant(theta)) - plant(q)*sin(plant(phi))*tan(plant(theta));
    d2fdx2_tensor(3, 3, 4) = plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 3, 10) = cos(plant(phi))*tan(plant(theta));
    d2fdx2_tensor(3, 3, 11) = -sin(plant(phi))*tan(plant(theta));
    d2fdx2_tensor(3, 4, 3) = plant(q)*cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1) - plant(r)*sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 4, 4) = 2*plant(r)*cos(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1) + 2*plant(q)*sin(plant(phi))*tan(plant(theta))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 4, 10) = sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 4, 11) = cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 10, 3) = cos(plant(phi))*tan(plant(theta));
    d2fdx2_tensor(3, 10, 4) = sin(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(3, 11, 3) = -sin(plant(phi))*tan(plant(theta));
    d2fdx2_tensor(3, 11, 4) = cos(plant(phi))*(std::pow(tan(plant(theta)), 2) + 1);
    d2fdx2_tensor(4, 3, 3) = plant(r)*sin(plant(phi)) - plant(q)*cos(plant(phi));
    d2fdx2_tensor(4, 3, 10) = -sin(plant(phi));
    d2fdx2_tensor(4, 3, 11) = -cos(plant(phi));
    d2fdx2_tensor(4, 10, 3) = -sin(plant(phi));
    d2fdx2_tensor(4, 11, 3) = -cos(plant(phi));
    d2fdx2_tensor(5, 3, 3) = - (plant(r)*cos(plant(phi)))/cos(plant(theta)) - (plant(q)*sin(plant(phi)))/cos(plant(theta));
    d2fdx2_tensor(5, 3, 4) = (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(5, 3, 10) = cos(plant(phi))/cos(plant(theta));
    d2fdx2_tensor(5, 3, 11) = -sin(plant(phi))/cos(plant(theta));
    d2fdx2_tensor(5, 4, 3) = (plant(q)*cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2) - (plant(r)*sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(5, 4, 4) = (plant(r)*cos(plant(phi)))/cos(plant(theta)) + (plant(q)*sin(plant(phi)))/cos(plant(theta)) + (2*plant(r)*cos(plant(phi))*std::pow(sin(plant(theta)), 2))/std::pow(cos(plant(theta)), 3) + (2*plant(q)*sin(plant(phi))*std::pow(sin(plant(theta)), 2))/std::pow(cos(plant(theta)), 3);
    d2fdx2_tensor(5, 4, 10) = (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(5, 4, 11) = (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(5, 10, 3) = cos(plant(phi))/cos(plant(theta));
    d2fdx2_tensor(5, 10, 4) = (sin(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(5, 11, 3) = -sin(plant(phi))/cos(plant(theta));
    d2fdx2_tensor(5, 11, 4) = (cos(plant(phi))*sin(plant(theta)))/std::pow(cos(plant(theta)), 2);
    d2fdx2_tensor(6, 3, 3) = -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(6, 3, 4) = -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    d2fdx2_tensor(6, 3, 5) = (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(6, 4, 3) = -(alge(ft)*cos(plant(psi))*cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    d2fdx2_tensor(6, 4, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(psi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(6, 4, 5) = -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2_tensor(6, 5, 3) = (alge(ft)*(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(6, 5, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2_tensor(6, 5, 5) = -(alge(ft)*(sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(7, 3, 3) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(7, 3, 4) = -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2_tensor(7, 3, 5) = (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(7, 4, 3) = -(alge(ft)*cos(plant(theta))*sin(plant(phi))*sin(plant(psi)))/m_droneParams.mass;
    d2fdx2_tensor(7, 4, 4) = -(alge(ft)*cos(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(7, 4, 5) = (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(7, 5, 3) = (alge(ft)*(cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(7, 5, 4) = (alge(ft)*cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(7, 5, 5) = (alge(ft)*(cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta))))/m_droneParams.mass;
    d2fdx2_tensor(8, 3, 3) = -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(8, 3, 4) = (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(8, 4, 3) = (alge(ft)*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(8, 4, 4) = -(alge(ft)*cos(plant(phi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdx2_tensor(9, 10, 11) = (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix;
    d2fdx2_tensor(9, 11, 10) = (m_droneParams.Iy - m_droneParams.Iz)/m_droneParams.Ix;
    d2fdx2_tensor(10, 9, 11) = -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy;
    d2fdx2_tensor(10, 11, 9) = -(m_droneParams.Ix - m_droneParams.Iz)/m_droneParams.Iy;
    d2fdx2_tensor(11, 9, 10) = (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz;
    d2fdx2_tensor(11, 10, 9) = (m_droneParams.Ix - m_droneParams.Iy)/m_droneParams.Iz;

    return d2fdx2_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2fdxdz(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_tensor(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES);
    d2fdxdz_tensor.setZero();

    d2fdxdz_tensor(6 , 3, ft) = ((std::cos(plant(phi))*std::sin(plant(psi)) - std::cos(plant(psi))*std::sin(plant(phi))*std::sin(plant(theta))))/m_droneParams.mass;
    d2fdxdz_tensor(6 , 4, ft) =  ( std::cos(plant(phi))*std::cos(plant(psi))*std::cos(plant(theta)))/m_droneParams.mass;
    d2fdxdz_tensor(6 , 5, ft) = ((std::cos(plant(psi))*std::sin(plant(phi)) - std::cos(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;

    d2fdxdz_tensor(7 , 3, ft) = -((std::cos(plant(phi))*std::cos(plant(psi)) + std::sin(plant(phi))*std::sin(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;
    d2fdxdz_tensor(7 , 4, ft) = (std::cos(plant(phi))*std::cos(plant(theta))*std::sin(plant(psi)))/m_droneParams.mass;
    d2fdxdz_tensor(7 , 5, ft) = ((std::sin(plant(phi))*std::sin(plant(psi)) + std::cos(plant(phi))*std::cos(plant(psi))*std::sin(plant(theta))))/m_droneParams.mass;

    d2fdxdz_tensor(8 , 3, ft) = -(std::cos(plant(theta))*std::sin(plant(phi)))/m_droneParams.mass;
    d2fdxdz_tensor(8 , 4, ft) = -(std::cos(plant(phi))*std::sin(plant(theta)))/m_droneParams.mass;

    return d2fdxdz_tensor;
}
Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2fdzdx(SystemState state)
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant = state.plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge = state.alge;
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_tensor(NUM_PLANT_STATES, NUM_Z_STATES, NUM_PLANT_STATES);
    d2fdzdx_tensor.setZero();

    d2fdzdx_tensor(xdot , ft, phi)= (cos(plant(phi))*sin(plant(psi)) - cos(plant(psi))*sin(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdzdx_tensor(xdot , ft, theta)= (cos(plant(phi))*cos(plant(psi))*cos(plant(theta)))/m_droneParams.mass;
    d2fdzdx_tensor(xdot , ft, psi)= (cos(plant(psi))*sin(plant(phi)) - cos(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass;

    d2fdzdx_tensor(ydot , ft, phi)= -(cos(plant(phi))*cos(plant(psi)) + sin(plant(phi))*sin(plant(psi))*sin(plant(theta)))/m_droneParams.mass;
    d2fdzdx_tensor(ydot , ft, theta)= (cos(plant(phi))*cos(plant(theta))*sin(plant(psi)))/m_droneParams.mass;
    d2fdzdx_tensor(ydot , ft, psi)= (sin(plant(phi))*sin(plant(psi)) + cos(plant(phi))*cos(plant(psi))*sin(plant(theta)))/m_droneParams.mass;

    d2fdzdx_tensor(zdot , ft, phi)= -(cos(plant(theta))*sin(plant(phi)))/m_droneParams.mass;
    d2fdzdx_tensor(zdot , ft, theta)= -(cos(plant(phi))*sin(plant(theta)))/m_droneParams.mass;
   
    return d2fdzdx_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2gdx2(SystemState state)
{
    
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdx2_tensor(NUM_Y_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2gdx2_tensor.setZero();

    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
            d2gdx2_tensor(0, xdot, psi)= -m_ctrlParams.at(velY).kp*m_sf(kpvy)*cos(state.plant(psi)) ;
            d2gdx2_tensor(0, ydot, psi)= -m_ctrlParams.at(velY).kp*m_sf(kpvy)*sin(state.plant(psi)) ;
            d2gdx2_tensor(0, psi, xdot)= -m_ctrlParams.at(velY).kp*m_sf(kpvy)*cos(state.plant(psi)) ;
            d2gdx2_tensor(0, psi, ydot)= -m_ctrlParams.at(velY).kp*m_sf(kpvy)*sin(state.plant(psi)) ;
            d2gdx2_tensor(0, psi, psi)=  -m_ctrlParams.at(velY).kp*m_sf(kpvy)*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi))) ;
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
            d2gdx2_tensor(1, xdot, psi)=  m_ctrlParams.at(velX).kp*m_sf(kpvx)*sin(state.plant(psi));
            d2gdx2_tensor(1, ydot, psi)= -m_ctrlParams.at(velX).kp*m_sf(kpvx)*cos(state.plant(psi));
            d2gdx2_tensor(1, psi, xdot)=   m_ctrlParams.at(velX).kp*m_sf(kpvx)*sin(state.plant(psi));
            d2gdx2_tensor(1, psi, ydot)=  -m_ctrlParams.at(velX).kp*m_sf(kpvx)*cos(state.plant(psi));
            d2gdx2_tensor(1, psi, psi)=   m_ctrlParams.at(velX).kp*m_sf(kpvx)*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)));
        
    }

    return d2gdx2_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdx2_plus(SystemState state, double time, double timestep )
{
    
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_tensor(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_tensor.setZero();

    d2hdx2_tensor(eix, x, psi)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eix, y, psi)= -timestep*cos(state.plant(psi));
    d2hdx2_tensor(eix, psi, x)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eix, psi, y)= -timestep*cos(state.plant(psi));
    d2hdx2_tensor(eix, psi, psi)= -timestep*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)));

    d2hdx2_tensor(edx, x, psi)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edx, y, psi)= -1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edx, psi, x)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edx, psi, y)= -1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edx, psi, psi)= (state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi)))/timestep;

    d2hdx2_tensor(desVelX,x, psi)= m_ctrlParams.at(posX).kp*m_sf(kppx)*sin(state.plant(psi));
    d2hdx2_tensor(desVelX,y, psi)= -m_ctrlParams.at(posX).kp*m_sf(kppx)*cos(state.plant(psi));
    d2hdx2_tensor(desVelX,psi, x)= m_ctrlParams.at(posX).kp*m_sf(kppx)*sin(state.plant(psi));
    d2hdx2_tensor(desVelX,psi, y)= -m_ctrlParams.at(posX).kp*m_sf(kppx)*cos(state.plant(psi));
    d2hdx2_tensor(desVelX,psi, psi)= -m_ctrlParams.at(posX).kp*m_sf(kppx)*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)));

    d2hdx2_tensor(eiy,x, psi)= timestep*cos(state.plant(psi));
    d2hdx2_tensor(eiy,y, psi)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eiy,psi, x)= timestep*cos(state.plant(psi));
    d2hdx2_tensor(eiy,psi, y)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eiy,psi, psi)= -timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)));

    d2hdx2_tensor(edy,x, psi)= 1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edy,y, psi)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edy,psi, x)= 1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edy,psi, y)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edy,psi, psi)= (state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi)))/timestep;

    d2hdx2_tensor(desVelY,x, psi)= m_ctrlParams.at(posY).kp*m_sf(kppy)*cos(state.plant(psi));
    d2hdx2_tensor(desVelY,y, psi)= m_ctrlParams.at(posY).kp*m_sf(kppy)*sin(state.plant(psi));
    d2hdx2_tensor(desVelY,psi, x)= m_ctrlParams.at(posY).kp*m_sf(kppy)*cos(state.plant(psi));
    d2hdx2_tensor(desVelY,psi, y)= m_ctrlParams.at(posY).kp*m_sf(kppy)*sin(state.plant(psi));
    d2hdx2_tensor(desVelY,psi, psi)= -m_ctrlParams.at(posY).kp*m_sf(kppy)*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)));

    d2hdx2_tensor(eixdot,xdot, psi)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eixdot,ydot, psi)= -timestep*cos(state.plant(psi));
    d2hdx2_tensor(eixdot,psi, xdot)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eixdot,psi, ydot)= -timestep*cos(state.plant(psi));
    d2hdx2_tensor(eixdot,psi, psi)= timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)));

    d2hdx2_tensor(edxdot,xdot, psi)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edxdot,ydot, psi)= -1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edxdot,psi, xdot)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edxdot,psi, ydot)= -1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edxdot,psi, psi)= (state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)))/timestep;

    d2hdx2_tensor(eiydot,xdot, psi)= timestep*cos(state.plant(psi));
    d2hdx2_tensor(eiydot,ydot, psi)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eiydot,psi, xdot)= timestep*cos(state.plant(psi));
    d2hdx2_tensor(eiydot,psi, ydot)= timestep*sin(state.plant(psi));
    d2hdx2_tensor(eiydot,psi, psi)= timestep*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)));

    d2hdx2_tensor(edydot,xdot, psi)= 1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edydot,ydot, psi)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edydot,psi, xdot)= 1/timestep*cos(state.plant(psi));
    d2hdx2_tensor(edydot,psi, ydot)= 1/timestep*sin(state.plant(psi));
    d2hdx2_tensor(edydot,psi, psi)= (state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi)))/timestep;

    return d2hdx2_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdx2_curr(SystemState prev, double timestep)
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_tensor(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_tensor.setZero();
    
    d2hdx2_tensor(edx,x, psi)= -1/timestep*sin(prev.plant(psi));
    d2hdx2_tensor(edx,y, psi)= 1/timestep*cos(prev.plant(psi));
    d2hdx2_tensor(edx,psi, x)= -sin(prev.plant(psi))/timestep;
    d2hdx2_tensor(edx,psi, y)= cos(prev.plant(psi))/timestep;
    d2hdx2_tensor(edx,psi, psi)= -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep;

    d2hdx2_tensor(edy,x, psi)= -1/timestep*cos(prev.plant(psi));
    d2hdx2_tensor(edy,y, psi)= -1/timestep*sin(prev.plant(psi));
    d2hdx2_tensor(edy,psi, x)= -cos(prev.plant(psi))/timestep;
    d2hdx2_tensor(edy,psi, y)= -sin(prev.plant(psi))/timestep;
    d2hdx2_tensor(edy,psi, psi)= -(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep;

    d2hdx2_tensor(edxdot,xdot, psi)= -1/timestep*sin(prev.plant(psi));
    d2hdx2_tensor(edxdot,ydot, psi)= 1/timestep*cos(prev.plant(psi));
    d2hdx2_tensor(edxdot,psi, xdot)= -sin(prev.plant(psi))/timestep;
    d2hdx2_tensor(edxdot,psi, ydot)= cos(prev.plant(psi))/timestep;
    d2hdx2_tensor(edxdot,psi, psi)= -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep;

    d2hdx2_tensor(edydot,xdot, psi)= -1/timestep*cos(prev.plant(psi));
    d2hdx2_tensor(edydot,ydot, psi)= -1/timestep*sin(prev.plant(psi));
    d2hdx2_tensor(edydot,psi, xdot)= -cos(prev.plant(psi))/timestep;
    d2hdx2_tensor(edydot,psi, ydot)= -sin(prev.plant(psi))/timestep;
    d2hdx2_tensor(edydot,psi, psi)= -(prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep;

    return d2hdx2_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdz2()
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdz2_plus_tensor(NUM_Z_STATES, NUM_Z_STATES, NUM_Z_STATES);
    d2hdz2_plus_tensor.setZero();

    d2hdz2_plus_tensor(ft, w1, w1)= m_droneParams.kf * 2 ;
    d2hdz2_plus_tensor(ft, w2, w2)= m_droneParams.kf * 2 ;
    d2hdz2_plus_tensor(ft, w3, w3)= m_droneParams.kf * 2 ;
    d2hdz2_plus_tensor(ft, w4, w4)= m_droneParams.kf * 2 ;

    d2hdz2_plus_tensor(tx, w1, w1)= - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(tx, w2, w2)= - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(tx, w3, w3)=   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(tx, w4, w4)=   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;

    d2hdz2_plus_tensor(ty, w1, w1)=   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(ty, w2, w2)= - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(ty, w3, w3)= - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdz2_plus_tensor(ty, w4, w4)=   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;

    d2hdz2_plus_tensor(tz, w1, w1)=    m_droneParams.km * 2  ;
    d2hdz2_plus_tensor(tz, w2, w2)=  - m_droneParams.km * 2  ;
    d2hdz2_plus_tensor(tz, w3, w3)=    m_droneParams.km * 2  ;
    d2hdz2_plus_tensor(tz, w4, w4)=  - m_droneParams.km * 2  ;

    return d2hdz2_plus_tensor;
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
    d2hdx2_plus.emplace_back(T(psi, psi, -timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x))) ));

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

    d2hdx2_plus.emplace_back(T(x, psi, m_ctrlParams.at(posX).kp*m_sf(kppx)*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, -m_ctrlParams.at(posX).kp*m_sf(kppx)*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, m_ctrlParams.at(posX).kp*m_sf(kppx)*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, -m_ctrlParams.at(posX).kp*m_sf(kppx)*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -m_ctrlParams.at(posX).kp*m_sf(kppx)*(cos(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) ));

    Eigen::SparseMatrix<double> d2hdx2_plus_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2hdx2_plus_mat.setFromTriplets(d2hdx2_plus.begin(), d2hdx2_plus.end());
    return d2hdx2_plus_mat;
}
Eigen::SparseMatrix<double> DroneTrajectory::d2desVely_dx_plus2(SystemState state, double time, double timestep)
{
    std::vector<T> d2hdx2_plus;
    d2hdx2_plus.reserve(5);

    d2hdx2_plus.emplace_back(T(x, psi, m_ctrlParams.at(posY).kp*m_sf(kppy)*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(y, psi, m_ctrlParams.at(posY).kp*m_sf(kppy)*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, x, m_ctrlParams.at(posY).kp*m_sf(kppy)*cos(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, y, m_ctrlParams.at(posY).kp*m_sf(kppy)*sin(state.plant(psi)) ));
    d2hdx2_plus.emplace_back(T(psi, psi, -m_ctrlParams.at(posY).kp*m_sf(kppy)*(cos(state.plant(psi))*(m_ref.at(refy)(time) -state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) -state.plant(x))) ));

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
        d2gdx2.push_back(T(xdot, psi, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*cos(state.plant(psi))));
        d2gdx2.push_back(T(ydot, psi, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*sin(state.plant(psi))));
        d2gdx2.push_back(T(psi, xdot, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*cos(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, ydot, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*sin(state.plant(psi))));
        d2gdx2.push_back(T(psi, psi, -m_ctrlParams.at(velY).kp*m_sf(kpvy)*(-state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi)))));
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
        d2gdx2.push_back(T(xdot, psi, m_ctrlParams.at(velX).kp*m_sf(kpvx)*sin(state.plant(psi))));
        d2gdx2.push_back(T(ydot, psi, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*cos(state.plant(psi))));
        d2gdx2.push_back(T(psi, xdot, m_ctrlParams.at(velX).kp*m_sf(kpvx)*sin(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, ydot, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*cos(state.plant(psi)) ));
        d2gdx2.push_back(T(psi, psi, -m_ctrlParams.at(velX).kp*m_sf(kpvx)*(-state.plant(xdot)*cos(state.plant(psi)) - state.plant(ydot)*sin(state.plant(psi)))));
    }
    Eigen::SparseMatrix<double> d2gdx2_mat(NUM_PLANT_STATES, NUM_PLANT_STATES);
    d2gdx2_mat.setFromTriplets(d2gdx2.begin(), d2gdx2.end());
    return d2gdx2_mat;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdxplus_dp(SystemState state, double time)
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdxplus_dp_tensor(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PARAMETERS);
    d2hdxplus_dp_tensor.setZero();

    d2hdxplus_dp_tensor(desVelX, x, kppx) = -cos(state.plant(psi))*m_ctrlParams.at(posX).kp;
    d2hdxplus_dp_tensor(desVelX, y, kppx) =  -sin(state.plant(psi))*m_ctrlParams.at(posX).kp;
    d2hdxplus_dp_tensor(desVelX, psi, kppx) = (cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))*m_ctrlParams.at(posX).kp;
    d2hdxplus_dp_tensor(desVelY, x, kppy) = sin(state.plant(psi))*m_ctrlParams.at(posY).kp;
    d2hdxplus_dp_tensor(desVelY, y, kppy) = -cos(state.plant(psi))*m_ctrlParams.at(posY).kp;
    d2hdxplus_dp_tensor(desVelY, psi, kppy) = (-(cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))))*m_ctrlParams.at(posY).kp ;
    d2hdxplus_dp_tensor(desVelZ, z, kppz) = -1*m_ctrlParams.at(posZ).kp;
    d2hdxplus_dp_tensor(desThrust, zdot, kpvz) = -m_thrustScale*m_ctrlParams.at(velZ).kp;
    d2hdxplus_dp_tensor(desRollRate, phi, kpphi) = -180/M_PI*m_ctrlParams.at(roll).kp;
    d2hdxplus_dp_tensor(desPitchRate, theta, kptheta) = -180/M_PI*m_ctrlParams.at(pitch).kp;
    d2hdxplus_dp_tensor(desYawRate, psi, kppsi) = -180/M_PI*m_ctrlParams.at(yaw).kp;
    d2hdxplus_dp_tensor(desRollOutput, p, kpp) = -180/M_PI*m_ctrlParams.at(rollRate).kp;
    d2hdxplus_dp_tensor(desPitchOutput, q, kpq) = -180/M_PI*m_ctrlParams.at(pitchRate).kp;
    d2hdxplus_dp_tensor(desYawOutput, r, kpr) = -180/M_PI*m_ctrlParams.at(yawRate).kp;

    return d2hdxplus_dp_tensor;
}
Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdzplus_dp()
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdzplus_dp_tensor(NUM_Z_STATES, NUM_Z_STATES, NUM_PARAMETERS);
    d2hdzplus_dp_tensor.setZero();

    d2hdzplus_dp_tensor(desVelX, eix, kipx) = 1*m_ctrlParams.at(posX).ki;
    d2hdzplus_dp_tensor(desVelX, edx, kdpx) = 1*m_ctrlParams.at(posX).kd;
    d2hdzplus_dp_tensor(desVelY, eiy, kipy) = 1*m_ctrlParams.at(posY).ki;
    d2hdzplus_dp_tensor(desVelY, edy, kdpy) = 1*m_ctrlParams.at(posY).kd;
    d2hdzplus_dp_tensor(desVelZ, eiz, kipz) = 1*m_ctrlParams.at(posZ).ki;
    d2hdzplus_dp_tensor(desVelZ, edz, kdpz) = 1*m_ctrlParams.at(posZ).kd;
    d2hdzplus_dp_tensor(desThrust, desVelZ, kpvz) = m_thrustScale*m_ctrlParams.at(velZ).kp;
    d2hdzplus_dp_tensor(desThrust, eizdot, kivz) = m_thrustScale*m_ctrlParams.at(velZ).ki;
    d2hdzplus_dp_tensor(desThrust, edzdot, kdvz) = m_thrustScale*m_ctrlParams.at(velZ).kd;
    d2hdzplus_dp_tensor(desRollRate, eiphi, kiphi) = 1*m_ctrlParams.at(roll).ki;
    d2hdzplus_dp_tensor(desRollRate, edphi, kdphi) = 1*m_ctrlParams.at(roll).kd;
    d2hdzplus_dp_tensor(desPitchRate, eitheta, kitheta) = 1*m_ctrlParams.at(pitch).ki;
    d2hdzplus_dp_tensor(desPitchRate, edtheta, kdtheta) = 1*m_ctrlParams.at(pitch).kd;
    d2hdzplus_dp_tensor(desYawRate, eipsi, kipsi) = 1*m_ctrlParams.at(yaw).ki;
    d2hdzplus_dp_tensor(desYawRate, edpsi, kdpsi) = 1*m_ctrlParams.at(yaw).kd;
    d2hdzplus_dp_tensor(desRollOutput, desRollRate, kpp ) = 1*m_ctrlParams.at(rollRate).kp;
    d2hdzplus_dp_tensor(desRollOutput, eip, kip) = 1*m_ctrlParams.at(rollRate).ki;
    d2hdzplus_dp_tensor(desRollOutput, edp, kdp) = 1*m_ctrlParams.at(rollRate).kd;
    d2hdzplus_dp_tensor(desPitchOutput, desPitchRate, kpq) = 1*m_ctrlParams.at(pitchRate).kp;
    d2hdzplus_dp_tensor(desPitchOutput, eiq, kiq) = 1*m_ctrlParams.at(pitchRate).ki;
    d2hdzplus_dp_tensor(desPitchOutput, edq, kdq) = 1*m_ctrlParams.at(pitchRate).kd;
    d2hdzplus_dp_tensor(desYawOutput, desYawRate, kpr) = 1*m_ctrlParams.at(yawRate).kp;
    d2hdzplus_dp_tensor(desYawOutput, eir, kir) = 1*m_ctrlParams.at(yawRate).ki;
    d2hdzplus_dp_tensor(desYawOutput, edr, kdr) = 1*m_ctrlParams.at(yawRate).kd;

    return d2hdzplus_dp_tensor;

}
Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2gdxplus_dp(SystemState state)
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdxplus_dp_tensor(NUM_Y_STATES, NUM_PLANT_STATES, NUM_PARAMETERS);
    d2gdxplus_dp_tensor.setZero();

    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        d2gdxplus_dp_tensor( 0, xdot, kpvy) = -sin(state.plant(psi))*m_ctrlParams.at(velY).kp;
        d2gdxplus_dp_tensor( 0, ydot, kpvy) = cos(state.plant(psi))*m_ctrlParams.at(velY).kp;
        d2gdxplus_dp_tensor( 0, psi, kpvy) = (-(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi))))*m_ctrlParams.at(velY).kp;
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        d2gdxplus_dp_tensor( 1, xdot, kpvx) = -cos(state.plant(psi))*m_ctrlParams.at(velX).kp;
        d2gdxplus_dp_tensor( 1, ydot, kpvx) = -sin(state.plant(psi))*m_ctrlParams.at(velX).kp;
        d2gdxplus_dp_tensor( 1, psi, kpvx) = (-(-state.plant(xdot)*sin(state.plant(psi)) + state.plant(ydot)*cos(state.plant(psi))))*m_ctrlParams.at(velX).kp;
    }

    return d2gdxplus_dp_tensor;
}
Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2gdzplus_dp(SystemState state)
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdzplus_dp_tensor(NUM_Y_STATES, NUM_Z_STATES, NUM_PARAMETERS);
    d2gdzplus_dp_tensor.setZero();
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        d2gdzplus_dp_tensor(0, desVelY, kpvy) = -1*m_ctrlParams.at(velY).kp;
        d2gdzplus_dp_tensor(0, eiydot, kivy) = -1*m_ctrlParams.at(velY).ki;
        d2gdzplus_dp_tensor(0, edydot,kdvy ) = -1*m_ctrlParams.at(velY).kd;
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        d2gdzplus_dp_tensor(1, desVelX, kpvx) = 1*m_ctrlParams.at(velX).kp;
        d2gdzplus_dp_tensor(1, eixdot, kivx) = 1*m_ctrlParams.at(velX).ki;
        d2gdzplus_dp_tensor(1, edxdot, kdvx) = 1*m_ctrlParams.at(velX).kd;
    }
    return d2gdzplus_dp_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdydp()
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdyplus_dp_tensor(NUM_Z_STATES, NUM_Y_STATES, NUM_PARAMETERS);
    d2hdyplus_dp_tensor.setZero();
    d2hdyplus_dp_tensor(desRollRate, 0, kpphi) = 1*m_ctrlParams.at(roll).kp;
    d2hdyplus_dp_tensor(desPitchRate, 1, kptheta ) = 1*m_ctrlParams.at(pitch).kp;
    return d2hdyplus_dp_tensor;
}

Eigen::Tensor<double, 3, Eigen::ColMajor> DroneTrajectory::d2hdzplus2(SystemState state)
{
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdzplus2_tensor(NUM_Z_STATES, NUM_Z_STATES, NUM_Z_STATES);
    d2hdzplus2_tensor.setZero();
    d2hdzplus2_tensor(ft, w1, w1) = m_droneParams.kf * 2 ;
    d2hdzplus2_tensor(ft, w2, w2) = m_droneParams.kf * 2 ;
    d2hdzplus2_tensor(ft, w3, w3) = m_droneParams.kf * 2 ;
    d2hdzplus2_tensor(ft, w4, w4) = m_droneParams.kf * 2 ;

    d2hdzplus2_tensor(tx, w1, w1) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(tx, w2, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(tx, w3, w3) =   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(tx, w4, w4) =   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;

    d2hdzplus2_tensor(ty, w1, w1) =   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(ty, w2, w2) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(ty, w3, w3) = - m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;
    d2hdzplus2_tensor(ty, w4, w4) =   m_droneParams.kf * m_droneParams.length * 2 * 1/sqrt(2) ;

    d2hdzplus2_tensor(tz, w1, w1) =    m_droneParams.km * 2  ;
    d2hdzplus2_tensor(tz, w2, w2) =  - m_droneParams.km * 2  ;
    d2hdzplus2_tensor(tz, w3, w3) =    m_droneParams.km * 2  ;
    d2hdzplus2_tensor(tz, w4, w4) =  - m_droneParams.km * 2  ;

    return d2hdzplus2_tensor;
}