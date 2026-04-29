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

Eigen::SparseMatrix<double> DroneTrajectory::dgdz(SystemState state, double timestep)
{
    std::vector<T> dgdz;
    dgdz.reserve(8);
    if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
        dgdz.push_back(T(0, eiydot, m_ctrlParams.at(velY).ki));
        dgdz.push_back(T(0, desVelY, m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep));
        dgdz.push_back(T(0, eiy, m_ctrlParams.at(posY).ki*(m_ctrlParams.at(velY).kp + timestep*m_ctrlParams.at(velY).ki)));
        dgdz.push_back(T(0, edy, m_ctrlParams.at(posY).kd*(m_ctrlParams.at(velY).kp + timestep*m_ctrlParams.at(velY).ki)));
    }
    if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
        dgdz.push_back(T(1, eixdot, m_ctrlParams.at(velX).ki));
        dgdz.push_back(T(1, desVelX, m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep));
        dgdz.push_back(T(1, eix, m_ctrlParams.at(posX).ki*(m_ctrlParams.at(velX).kp + timestep*m_ctrlParams.at(velX).ki)));
        dgdz.push_back(T(1, edx, m_ctrlParams.at(posX).kd*(m_ctrlParams.at(velX).kp + timestep*m_ctrlParams.at(velX).ki)));
    }
    Eigen::SparseMatrix<double> dgdz_mat(NUM_Y_STATES, NUM_Z_STATES);
    dgdz_mat.setFromTriplets(dgdz.begin(), dgdz.end());
    return dgdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxPlus(SystemState state, double time, double timestep)
{
    std::vector<T> dhdx;
    dhdx.reserve(150);

    // TODO: what if ref is not constant?
    dhdx.push_back(T(eix, 0, -timestep*std::cos(state.plant(psi))));
    dhdx.push_back(T(eix, 1, -timestep*std::sin(state.plant(psi))));
    dhdx.push_back(T(eix, 5, timestep*(std::cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - std::sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))));
    dhdx.push_back(T(edx, 0, -std::cos(state.plant(psi))/timestep));
    dhdx.push_back(T(edx, 1, -std::sin(state.plant(psi))/timestep));
    dhdx.push_back(T(edx, 5, -(state.plant(y)*std::cos(state.plant(psi)) - state.plant(x)*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 0, - m_ctrlParams.at(posX).kp*std::cos(state.plant(psi)) - m_ctrlParams.at(posX).ki*timestep*std::cos(state.plant(psi)) - (m_ctrlParams.at(posX).kd*std::cos(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 1, - m_ctrlParams.at(posX).kp*std::sin(state.plant(psi)) - m_ctrlParams.at(posX).ki*timestep*std::sin(state.plant(psi)) - (m_ctrlParams.at(posX).kd*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 5, m_ctrlParams.at(posX).kp*(std::cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - std::sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) - (m_ctrlParams.at(posX).kd*(state.plant(y)*std::cos(state.plant(psi)) - state.plant(x)*std::sin(state.plant(psi))))/timestep + m_ctrlParams.at(posX).ki*timestep*(std::cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - std::sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))));
    dhdx.push_back(T(eiy, 1, -timestep*(std::cos(state.plant(psi)) - std::sin(state.plant(psi)))));
    dhdx.push_back(T(eiy, 5, -timestep*(std::cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(y)) + std::sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)))));
    dhdx.push_back(T(edy, 0, std::sin(state.plant(psi))/timestep));
    dhdx.push_back(T(edy, 1, -std::cos(state.plant(psi))/timestep));
    dhdx.push_back(T(edy, 5, (state.plant(x)*std::cos(state.plant(psi)) + state.plant(y)*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, 0, m_ctrlParams.at(posY).kp*std::sin(state.plant(psi)) + (m_ctrlParams.at(posY).kd*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, 1, - m_ctrlParams.at(posY).kp*std::cos(state.plant(psi)) - (m_ctrlParams.at(posY).kd*std::cos(state.plant(psi)))/timestep - m_ctrlParams.at(posY).ki*timestep*(std::cos(state.plant(psi)) - std::sin(state.plant(psi)))));
    dhdx.push_back(T(desVelY, 5, (m_ctrlParams.at(posY).kd*(state.plant(x)*std::cos(state.plant(psi)) + state.plant(y)*std::sin(state.plant(psi))))/timestep - m_ctrlParams.at(posY).kp*(std::cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + std::sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) - m_ctrlParams.at(posY).ki*timestep*(std::cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(y)) + std::sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)))));
    dhdx.push_back(T(eiz, 2, -timestep));
    dhdx.push_back(T(edz, 2, -1/timestep));
    dhdx.push_back(T(desVelZ, 2, - m_ctrlParams.at(posZ).kp - m_ctrlParams.at(posZ).ki*timestep - m_ctrlParams.at(posZ).kd/timestep));
    dhdx.push_back(T(eixdot, 0, -timestep*(m_ctrlParams.at(posX).kp*std::cos(state.plant(psi)) + m_ctrlParams.at(posX).ki*timestep*std::cos(state.plant(psi)) + (m_ctrlParams.at(posX).kd*std::cos(state.plant(psi)))/timestep)));
    dhdx.push_back(T(eixdot, 1, -timestep*(m_ctrlParams.at(posX).kp*std::sin(state.plant(psi)) + m_ctrlParams.at(posX).ki*timestep*std::sin(state.plant(psi)) + (m_ctrlParams.at(posX).kd*std::sin(state.plant(psi)))/timestep)));
    dhdx.push_back(T(eixdot, 5, timestep*(state.plant(xdot)*std::sin(state.plant(psi)) - state.plant(ydot)*std::cos(state.plant(psi)) + m_ctrlParams.at(posX).kp*(std::cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - std::sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) - (m_ctrlParams.at(posX).kd*(state.plant(y)*std::cos(state.plant(psi)) - state.plant(x)*std::sin(state.plant(psi))))/timestep + m_ctrlParams.at(posX).ki*timestep*(std::cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - std::sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))))));
    dhdx.push_back(T(eixdot, 6, -timestep*std::cos(state.plant(psi))));
    dhdx.push_back(T(eixdot, 7, -timestep*std::sin(state.plant(psi))));
    dhdx.push_back(T(edxdot, 5, -(state.plant(ydot)*std::cos(state.plant(psi)) - state.plant(xdot)*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(edxdot, 6, -std::cos(state.plant(psi))/timestep));
    dhdx.push_back(T(edxdot, 7, -std::sin(state.plant(psi))/timestep));
    dhdx.push_back(T(eiydot, 0, timestep*(m_ctrlParams.at(posY).kp*std::sin(state.plant(psi)) + (m_ctrlParams.at(posY).kd*std::sin(state.plant(psi)))/timestep)));
    dhdx.push_back(T(eiydot, 1, -timestep*(m_ctrlParams.at(posY).kp*std::cos(state.plant(psi)) + (m_ctrlParams.at(posY).kd*std::cos(state.plant(psi)))/timestep + m_ctrlParams.at(posY).ki*timestep*(std::cos(state.plant(psi)) - std::sin(state.plant(psi))))));
    dhdx.push_back(T(eiydot, 5, timestep*(state.plant(xdot)*std::cos(state.plant(psi)) + state.plant(ydot)*std::sin(state.plant(psi)) - m_ctrlParams.at(posY).kp*(std::cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)) + std::sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))) + (m_ctrlParams.at(posY).kd*(state.plant(x)*std::cos(state.plant(psi)) + state.plant(y)*std::sin(state.plant(psi))))/timestep - m_ctrlParams.at(posY).ki*timestep*(std::cos(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(y)) + std::sin(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y))))));
    dhdx.push_back(T(eiydot, 6, timestep*std::sin(state.plant(psi))));
    dhdx.push_back(T(eiydot, 7, -timestep*std::cos(state.plant(psi))));
    dhdx.push_back(T(edydot, 5, (state.plant(xdot)*std::cos(state.plant(psi)) + state.plant(ydot)*std::sin(state.plant(psi)))/timestep));
    dhdx.push_back(T(edydot, 6, std::sin(state.plant(psi))/timestep));
    dhdx.push_back(T(edydot, 7, -std::cos(state.plant(psi))/timestep));
    dhdx.push_back(T(eizdot, 2, -timestep*(m_ctrlParams.at(posZ).kp + m_ctrlParams.at(posZ).ki*timestep + m_ctrlParams.at(posZ).kd/timestep)));
    dhdx.push_back(T(eizdot, 8, -timestep));
    dhdx.push_back(T(edzdot, 8, -1/timestep));
    dhdx.push_back(T(desThrust, 2, -m_thrustScale*(m_ctrlParams.at(velZ).kp*(m_ctrlParams.at(posZ).kp + m_ctrlParams.at(posZ).ki*timestep + m_ctrlParams.at(posZ).kd/timestep) + m_ctrlParams.at(velZ).ki*timestep*(m_ctrlParams.at(posZ).kp + m_ctrlParams.at(posZ).ki*timestep + m_ctrlParams.at(posZ).kd/timestep))));
    dhdx.push_back(T(desThrust, 8, -m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep + m_ctrlParams.at(velZ).kd/timestep)));
    dhdx.push_back(T(eiphi, 3, -(180*timestep)/M_PI));
    dhdx.push_back(T(edphi, 3, -180/(M_PI*timestep)));
    dhdx.push_back(T(desRollRate, 3, - (180*m_ctrlParams.at(roll).kp)/M_PI - (180*m_ctrlParams.at(roll).kd)/(M_PI*timestep) - (180*m_ctrlParams.at(roll).ki*timestep)/M_PI));
    dhdx.push_back(T(eitheta, 4, -(180*timestep)/M_PI));
    dhdx.push_back(T(edtheta, 4, -180/(M_PI*timestep)));
    dhdx.push_back(T(desPitchRate, 4, - (180*m_ctrlParams.at(pitch).kp)/M_PI - (180*m_ctrlParams.at(pitch).kd)/(M_PI*timestep) - (180*m_ctrlParams.at(pitch).ki*timestep)/M_PI));
    dhdx.push_back(T(eipsi, 5, -(180*timestep)/M_PI));
    dhdx.push_back(T(edpsi, 5, -180/(M_PI*timestep)));
    dhdx.push_back(T(desYawRate, 5, - (180*m_ctrlParams.at(yaw).kp)/M_PI - (180*m_ctrlParams.at(yaw).kd)/(M_PI*timestep) - (180*m_ctrlParams.at(yaw).ki*timestep)/M_PI));
    dhdx.push_back(T(eip, 3, -timestep*((180*m_ctrlParams.at(roll).kp)/M_PI + (180*m_ctrlParams.at(roll).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(roll).ki*timestep)/M_PI)));
    dhdx.push_back(T(eip, 9, -(180*timestep)/M_PI));
    dhdx.push_back(T(edp, 9, -(180*lpf_b0)/(M_PI*timestep)));
    dhdx.push_back(T(desRollOutput, 3, - m_ctrlParams.at(rollRate).kp*((180*m_ctrlParams.at(roll).kp)/M_PI + (180*m_ctrlParams.at(roll).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(roll).ki*timestep)/M_PI) - m_ctrlParams.at(rollRate).ki*timestep*((180*m_ctrlParams.at(roll).kp)/M_PI + (180*m_ctrlParams.at(roll).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(roll).ki*timestep)/M_PI)));
    dhdx.push_back(T(desRollOutput, 9, - (180*m_ctrlParams.at(rollRate).kp)/M_PI - (180*m_ctrlParams.at(rollRate).ki*timestep)/M_PI - (180*lpf_b0*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(eiq, 4, -timestep*((180*m_ctrlParams.at(pitch).kp)/M_PI + (180*m_ctrlParams.at(pitch).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(pitch).ki*timestep)/M_PI)));
    dhdx.push_back(T(eiq, 10, -(180*timestep)/M_PI));
    dhdx.push_back(T(edq, 10, -(180*lpf_b0)/(M_PI*timestep)));
    dhdx.push_back(T(desPitchOutput, 4, - m_ctrlParams.at(pitchRate).kp*((180*m_ctrlParams.at(pitch).kp)/M_PI + (180*m_ctrlParams.at(pitch).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(pitch).ki*timestep)/M_PI) - m_ctrlParams.at(pitchRate).ki*timestep*((180*m_ctrlParams.at(pitch).kp)/M_PI + (180*m_ctrlParams.at(pitch).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(pitch).ki*timestep)/M_PI)));
    dhdx.push_back(T(desPitchOutput, 10, - (180*m_ctrlParams.at(pitchRate).kp)/M_PI - (180*m_ctrlParams.at(pitchRate).ki*timestep)/M_PI - (180*lpf_b0*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(eir, 5, -timestep*((180*m_ctrlParams.at(yaw).kp)/M_PI + (180*m_ctrlParams.at(yaw).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(yaw).ki*timestep)/M_PI)));
    dhdx.push_back(T(eir, 11, -(180*timestep)/M_PI));
    dhdx.push_back(T(edr, 11, -180/(M_PI*timestep)));
    dhdx.push_back(T(desYawOutput, 5, - m_ctrlParams.at(yawRate).kp*((180*m_ctrlParams.at(yaw).kp)/M_PI + (180*m_ctrlParams.at(yaw).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(yaw).ki*timestep)/M_PI) - m_ctrlParams.at(yawRate).ki*timestep*((180*m_ctrlParams.at(yaw).kp)/M_PI + (180*m_ctrlParams.at(yaw).kd)/(M_PI*timestep) + (180*m_ctrlParams.at(yaw).ki*timestep)/M_PI)));
    dhdx.push_back(T(desYawOutput, 11, - (180*m_ctrlParams.at(yawRate).kp)/M_PI - (180*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep) - (180*m_ctrlParams.at(yawRate).ki*timestep)/M_PI));
    dhdx.push_back(T(delay_1_rollRate, 9, -180/(M_PI*timestep)));
    dhdx.push_back(T(delay_1_pitchRate, 10, -180/(M_PI*timestep)));
    dhdx.push_back(T(w1, 2, -(m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)*(m_ctrlParams.at(posZ).kd + m_ctrlParams.at(posZ).kp*timestep + m_ctrlParams.at(posZ).ki*std::pow(timestep, 2)))/timestep));
    dhdx.push_back(T(w1, 3, (90*m_alpha*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(roll).kd + m_ctrlParams.at(roll).kp*timestep + m_ctrlParams.at(roll).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 4, -(90*m_alpha*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(pitch).kd + m_ctrlParams.at(pitch).kp*timestep + m_ctrlParams.at(pitch).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 8, -m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep + m_ctrlParams.at(velZ).kd/timestep)));
    dhdx.push_back(T(w1, 9, (90*m_alpha*(lpf_b0*m_ctrlParams.at(rollRate).kd + m_ctrlParams.at(rollRate).kp*timestep + m_ctrlParams.at(rollRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 10, -(90*m_alpha*(lpf_b0*m_ctrlParams.at(pitchRate).kd + m_ctrlParams.at(pitchRate).kp*timestep + m_ctrlParams.at(pitchRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 11, -(180*m_alpha*(m_ctrlParams.at(yawRate).kd + m_ctrlParams.at(yawRate).kp*timestep + m_ctrlParams.at(yawRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 2, -(m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)*(m_ctrlParams.at(posZ).kd + m_ctrlParams.at(posZ).kp*timestep + m_ctrlParams.at(posZ).ki*std::pow(timestep, 2)))/timestep));
    dhdx.push_back(T(w2, 3, (90*m_alpha*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(roll).kd + m_ctrlParams.at(roll).kp*timestep + m_ctrlParams.at(roll).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 4, (90*m_alpha*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(pitch).kd + m_ctrlParams.at(pitch).kp*timestep + m_ctrlParams.at(pitch).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 5, (180*m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)*(m_ctrlParams.at(yaw).kd + m_ctrlParams.at(yaw).kp*timestep + m_ctrlParams.at(yaw).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 8, -m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep + m_ctrlParams.at(velZ).kd/timestep)));
    dhdx.push_back(T(w2, 9, (90*m_alpha*(lpf_b0*m_ctrlParams.at(rollRate).kd + m_ctrlParams.at(rollRate).kp*timestep + m_ctrlParams.at(rollRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 10, (90*m_alpha*(lpf_b0*m_ctrlParams.at(pitchRate).kd + m_ctrlParams.at(pitchRate).kp*timestep + m_ctrlParams.at(pitchRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 11, (180*m_alpha*(m_ctrlParams.at(yawRate).kd + m_ctrlParams.at(yawRate).kp*timestep + m_ctrlParams.at(yawRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 2, -(m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)*(m_ctrlParams.at(posZ).kd + m_ctrlParams.at(posZ).kp*timestep + m_ctrlParams.at(posZ).ki*std::pow(timestep, 2)))/timestep));
    dhdx.push_back(T(w3, 3, -(90*m_alpha*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(roll).kd + m_ctrlParams.at(roll).kp*timestep + m_ctrlParams.at(roll).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 4, (90*m_alpha*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(pitch).kd + m_ctrlParams.at(pitch).kp*timestep + m_ctrlParams.at(pitch).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 8, -m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep + m_ctrlParams.at(velZ).kd/timestep)));
    dhdx.push_back(T(w3, 9, -(90*m_alpha*(lpf_b0*m_ctrlParams.at(rollRate).kd + m_ctrlParams.at(rollRate).kp*timestep + m_ctrlParams.at(rollRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 10, (90*m_alpha*(lpf_b0*m_ctrlParams.at(pitchRate).kd + m_ctrlParams.at(pitchRate).kp*timestep + m_ctrlParams.at(pitchRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 11, -(180*m_alpha*(m_ctrlParams.at(yawRate).kd + m_ctrlParams.at(yawRate).kp*timestep + m_ctrlParams.at(yawRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 2, -(m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)*(m_ctrlParams.at(posZ).kd + m_ctrlParams.at(posZ).kp*timestep + m_ctrlParams.at(posZ).ki*std::pow(timestep, 2)))/timestep));
    dhdx.push_back(T(w4, 3, -(90*m_alpha*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(roll).kd + m_ctrlParams.at(roll).kp*timestep + m_ctrlParams.at(roll).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 4, -(90*m_alpha*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(pitch).kd + m_ctrlParams.at(pitch).kp*timestep + m_ctrlParams.at(pitch).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 8, -m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep + m_ctrlParams.at(velZ).kd/timestep)));
    dhdx.push_back(T(w4, 9, -(90*m_alpha*(lpf_b0*m_ctrlParams.at(rollRate).kd + m_ctrlParams.at(rollRate).kp*timestep + m_ctrlParams.at(rollRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 10, -(90*m_alpha*(lpf_b0*m_ctrlParams.at(pitchRate).kd + m_ctrlParams.at(pitchRate).kp*timestep + m_ctrlParams.at(pitchRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 11, (180*m_alpha*(m_ctrlParams.at(yawRate).kd + m_ctrlParams.at(yawRate).kp*timestep + m_ctrlParams.at(yawRate).ki*std::pow(timestep, 2)))/(M_PI*timestep)));

    if (state.alge(desRoll) > -20  && state.alge(desRoll) < 20)
    {
        dhdx.push_back(T(eiphi, 0, -sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiphi, 1, cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiphi, 5, -timestep*(m_ctrlParams.at(velY).kp*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)) - m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(refx - state.plant(x)) + sin(state.plant(psi))*(refy - state.plant(y))) + (m_ctrlParams.at(posY).kd*(state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi))))/timestep - m_ctrlParams.at(posY).ki*timestep*(cos(state.plant(psi))*(refx - state.plant(x)) + sin(state.plant(psi))*(refy - state.plant(y)))) + (m_ctrlParams.at(velY).kd*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi))))/timestep + m_ctrlParams.at(velY).ki*timestep*(state.plant(xdot)*cos(state.plant(psi)) + state.plant(ydot)*sin(state.plant(psi)) - m_ctrlParams.at(posY).kp*(cos(state.plant(psi))*(refx - state.plant(x)) + sin(state.plant(psi))*(refy - state.plant(y))) + (m_ctrlParams.at(posY).kd*(state.plant(x)*cos(state.plant(psi)) + state.plant(y)*sin(state.plant(psi))))/timestep - m_ctrlParams.at(posY).ki*timestep*(cos(state.plant(psi))*(refx - state.plant(x)) + sin(state.plant(psi))*(refy - state.plant(y)))))));
        dhdx.push_back(T(eiphi, 6, -sin(state.plant(psi))*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiphi, 7, cos(state.plant(psi))*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(desRollRate, 0, -(sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollRate, 1, (cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollRate, 5, -((m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kd*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(ydot)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refy*timestep*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refx*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 3)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refy*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refx*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(x)*cos(state.plant(psi))))/timestep));
        dhdx.push_back(T(desRollRate, 6, -(sin(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollRate, 7, (cos(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(eip, 0, -sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eip, 1, cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eip, 5, -(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kd*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(ydot)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refy*timestep*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refx*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 3)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refy*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refx*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(x)*cos(state.plant(psi)))));
        dhdx.push_back(T(eip, 6, -sin(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eip, 7, cos(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(desRollOutput, 0, -(sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollOutput, 1, (cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollOutput, 5, -((m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(velY).kd*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kd*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*timestep*state.plant(ydot)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(xdot)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)*state.plant(ydot)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refy*timestep*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refx*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refx*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(x)*cos(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 3)*sin(state.plant(psi)) - m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*refy*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*refy*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*std::pow(timestep, 3)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kp*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki*std::pow(timestep, 2)*state.plant(y)*sin(state.plant(psi)) - m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*refx*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep*state.plant(x)*cos(state.plant(psi)) + m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).kp*timestep*state.plant(x)*cos(state.plant(psi))))/timestep));
        dhdx.push_back(T(desRollOutput, 6, -(sin(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desRollOutput, 7, (cos(state.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(velY).kd + m_ctrlParams.at(velY).kp*timestep + m_ctrlParams.at(velY).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(w1, 0, (m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w1, 1, -(m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w2, 0, (m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w2, 1, -(m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w3, 0, -(m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w3, 1, (m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w4, 0, -(m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w4, 1, (m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(posY).kd + m_ctrlParams.at(posY).kp*timestep + m_ctrlParams.at(posY).ki*std::pow(timestep, 2)))/(2*timestep)));
    }
    if (state.alge(desPitch) > -20  && state.alge(desPitch) < 20)
    {
        dhdx.push_back(T(eitheta, 0, -cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eitheta, 1, -sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eitheta, 5, timestep*(m_ctrlParams.at(velX).kp*(state.plant(xdot)*sin(state.plant(psi)) - state.plant(ydot)*cos(state.plant(psi)) + m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) - (m_ctrlParams.at(posX).kd*(state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi))))/timestep + m_ctrlParams.at(posX).ki*timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))) - (m_ctrlParams.at(velX).kd*(state.plant(ydot)*cos(state.plant(psi)) - state.plant(xdot)*sin(state.plant(psi))))/timestep + m_ctrlParams.at(velX).ki*timestep*(state.plant(xdot)*sin(state.plant(psi)) - state.plant(ydot)*cos(state.plant(psi)) + m_ctrlParams.at(posX).kp*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x))) - (m_ctrlParams.at(posX).kd*(state.plant(y)*cos(state.plant(psi)) - state.plant(x)*sin(state.plant(psi))))/timestep + m_ctrlParams.at(posX).ki*timestep*(cos(state.plant(psi))*(m_ref.at(refy)(time) - state.plant(y)) - sin(state.plant(psi))*(m_ref.at(refx)(time) - state.plant(x)))))));
        dhdx.push_back(T(eitheta, 6, -cos(state.plant(psi))*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eitheta, 7, -sin(state.plant(psi))*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(desPitchRate, 0, -(cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchRate, 1, -(sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchRate, 5, -((m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kd*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*timestep*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*timestep*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*timestep*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 3)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(y)*cos(state.plant(psi))))/timestep));
        dhdx.push_back(T(desPitchRate, 6, -(cos(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchRate, 7, -(sin(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(eiq, 0, -cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiq, 1, -sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiq, 5, -(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kd*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*timestep*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*timestep*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*timestep*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 3)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(y)*cos(state.plant(psi)))));
        dhdx.push_back(T(eiq, 6, -cos(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(eiq, 7, -sin(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2))));
        dhdx.push_back(T(desPitchOutput, 0, -(cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchOutput, 1, -(sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchOutput, 5, -((m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(velX).kd*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kd*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*timestep*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*timestep*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(ydot)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*std::pow(timestep, 2)*state.plant(xdot)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*timestep*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 3)*cos(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refy)(time)*std::pow(timestep, 2)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 3)*sin(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*m_ref.at(refx)(time)*std::pow(timestep, 2)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*std::pow(timestep, 3)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kp*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki*std::pow(timestep, 2)*state.plant(x)*sin(state.plant(psi)) - m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*m_ref.at(refy)(time)*timestep*cos(state.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep*state.plant(y)*cos(state.plant(psi)) + m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).kp*timestep*state.plant(y)*cos(state.plant(psi))))/timestep));
        dhdx.push_back(T(desPitchOutput, 6, -(cos(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(desPitchOutput, 7, -(sin(state.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(velX).kd + m_ctrlParams.at(velX).kp*timestep + m_ctrlParams.at(velX).ki*std::pow(timestep, 2)))/timestep));
        dhdx.push_back(T(w1, 0, -(m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w1, 1, -(m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w2, 0, (m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w2, 1, (m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w3, 0, (m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w3, 1, (m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w4, 0, -(m_alpha*cos(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
        dhdx.push_back(T(w4, 1, -(m_alpha*sin(state.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(posX).kd + m_ctrlParams.at(posX).kp*timestep + m_ctrlParams.at(posX).ki*std::pow(timestep, 2)))/(2*timestep)));
    }
    if (state.alge(desRoll) == -20  || state.alge(desRoll) == 20 || state.alge(desPitch) == -20  || state.alge(desPitch) == 20) 
    {
        dhdx.push_back(T(w1, 5, -(180*m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)*(m_ctrlParams.at(yaw).kd + m_ctrlParams.at(yaw).kp*timestep + m_ctrlParams.at(yaw).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
        dhdx.push_back(T(w2, 5, (180*m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)*(m_ctrlParams.at(yaw).kd + m_ctrlParams.at(yaw).kp*timestep + m_ctrlParams.at(yaw).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
        dhdx.push_back(T(w3, 5, -(180*m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)*(m_ctrlParams.at(yaw).kd + m_ctrlParams.at(yaw).kp*timestep + m_ctrlParams.at(yaw).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
        dhdx.push_back(T(w4, 5, (180*m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)*(m_ctrlParams.at(yaw).kd + m_ctrlParams.at(yaw).kp*timestep + m_ctrlParams.at(yaw).ki*std::pow(timestep, 2)))/(M_PI*timestep)));
    }
    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());

    if ((state.alge(desRoll) > -20  && state.alge(desRoll) < 20) || (state.alge(desPitch) > -20  && state.alge(desPitch) < 20))
    {
        dhdx.push_back(T(w1, 5, m_alpha*(-0.5*dhdx_mat.coeff(desRollOutput, 5) + 0.5*dhdx_mat.coeff(desPitchOutput, 5) + dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w2, 5, m_alpha*(-0.5*dhdx_mat.coeff(desRollOutput, 5) - 0.5*dhdx_mat.coeff(desPitchOutput, 5) - dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w3, 5, m_alpha*(+0.5*dhdx_mat.coeff(desRollOutput, 5) - 0.5*dhdx_mat.coeff(desPitchOutput, 5) + dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w4, 5, m_alpha*(+0.5*dhdx_mat.coeff(desRollOutput, 5) + 0.5*dhdx_mat.coeff(desPitchOutput, 5) - dhdx_mat.coeff(desYawOutput, 5))));
    }

    Eigen::Vector<double, NUM_PLANT_STATES> ft_vec = m_droneParams.kf * 2 * (state.alge(w1)*dhdx_mat.row(w1) 
                                                                        + state.alge(w2)*dhdx_mat.row(w2)
                                                                        + state.alge(w3)*dhdx_mat.row(w3)
                                                                        + state.alge(w4)*dhdx_mat.row(w4));
                                             
    Eigen::Vector<double, NUM_PLANT_STATES> tx_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (- state.alge(w1)*dhdx_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdx_mat.row(w2) 
                                                                                                               + state.alge(w3)*dhdx_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdx_mat.row(w4));
    Eigen::Vector<double, NUM_PLANT_STATES> ty_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (  state.alge(w1)*dhdx_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdx_mat.row(w2) 
                                                                                                               - state.alge(w3)*dhdx_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdx_mat.row(w4));
    Eigen::Vector<double, NUM_PLANT_STATES> tz_vec = m_droneParams.km * 2 * (state.alge(w1)*dhdx_mat.row(w1) 
                                                                       - state.alge(w2)*dhdx_mat.row(w2) 
                                                                       + state.alge(w3)*dhdx_mat.row(w3) 
                                                                       - state.alge(w4)*dhdx_mat.row(w4));

    // Append  entries to the original triplets
    for (int col = 0; col < NUM_PLANT_STATES; ++col) {
        if (ft_vec(col) != 0.0) {
            dhdx.emplace_back(ft, col, ft_vec(col));
        }
        if (tx_vec(col) != 0.0) {
            dhdx.emplace_back(tx, col, tx_vec(col));
        }
        if (ty_vec(col) != 0.0) {
            dhdx.emplace_back(ty, col, ty_vec(col));
        }
        if (tz_vec(col) != 0.0) {
            dhdx.emplace_back(tz, col, tz_vec(col));
        }
    }

    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxCurr(SystemState curr, SystemState prev, double timestep)
{
    // TODO: what if ref is not constant?
    std::vector<T> dhdx;
    dhdx.reserve(10);
    
    dhdx.push_back(T(edx, 0, cos(prev.plant(psi))/timestep));
    dhdx.push_back(T(edx, 1, sin(prev.plant(psi))/timestep));
    dhdx.push_back(T(edx, 5, (prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 0, (m_ctrlParams.at(posX).kd*cos(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 1, (m_ctrlParams.at(posX).kd*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelX, 5, (m_ctrlParams.at(posX).kd*(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi))))/timestep));
    dhdx.push_back(T(edy, 0, -sin(prev.plant(psi))/timestep));
    dhdx.push_back(T(edy, 1, cos(prev.plant(psi))/timestep));
    dhdx.push_back(T(edy, 5, -(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, 0, -(m_ctrlParams.at(posY).kd*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, 1, (m_ctrlParams.at(posY).kd*cos(prev.plant(psi)))/timestep));
    dhdx.push_back(T(desVelY, 5, -(m_ctrlParams.at(posY).kd*(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi))))/timestep));
    dhdx.push_back(T(edz, 2, 1/timestep));
    dhdx.push_back(T(desVelZ, 2, m_ctrlParams.at(posZ).kd/timestep));
    dhdx.push_back(T(eixdot, 0, m_ctrlParams.at(posX).kd*cos(prev.plant(psi))));
    dhdx.push_back(T(eixdot, 1, m_ctrlParams.at(posX).kd*sin(prev.plant(psi))));
    dhdx.push_back(T(eixdot, 5, m_ctrlParams.at(posX).kd*(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi)))));
    dhdx.push_back(T(edxdot, 5, (prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(edxdot, 6, cos(prev.plant(psi))/timestep));
    dhdx.push_back(T(edxdot, 7, sin(prev.plant(psi))/timestep));
    dhdx.push_back(T(eiydot, 0, -m_ctrlParams.at(posY).kd*sin(prev.plant(psi))));
    dhdx.push_back(T(eiydot, 1, m_ctrlParams.at(posY).kd*cos(prev.plant(psi))));
    dhdx.push_back(T(eiydot, 5, -m_ctrlParams.at(posY).kd*(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi)))));
    dhdx.push_back(T(edydot, 5, -(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi)))/timestep));
    dhdx.push_back(T(edydot, 6, -sin(prev.plant(psi))/timestep));
    dhdx.push_back(T(edydot, 7, cos(prev.plant(psi))/timestep));
    dhdx.push_back(T(eizdot, 2, m_ctrlParams.at(posZ).kd));
    dhdx.push_back(T(edzdot, 8, 1/timestep));
    dhdx.push_back(T(desThrust, 2, (m_ctrlParams.at(posZ).kd*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep))/timestep));
    dhdx.push_back(T(desThrust, 8, (m_ctrlParams.at(velZ).kd*m_thrustScale)/timestep));
    dhdx.push_back(T(edphi, 3, 180/(M_PI*timestep)));
    dhdx.push_back(T(desRollRate, 3, (180*m_ctrlParams.at(roll).kd)/(M_PI*timestep)));
    dhdx.push_back(T(edtheta, 4, 180/(M_PI*timestep)));
    dhdx.push_back(T(desPitchRate, 4, (180*m_ctrlParams.at(pitch).kd)/(M_PI*timestep)));
    dhdx.push_back(T(edpsi, 5, 180/(M_PI*timestep)));
    dhdx.push_back(T(desYawRate, 5, (180*m_ctrlParams.at(yaw).kd)/(M_PI*timestep)));
    dhdx.push_back(T(eip, 3, (180*m_ctrlParams.at(roll).kd)/M_PI));
    dhdx.push_back(T(edp, 9, (180*lpf_b0)/(M_PI*timestep)));
    dhdx.push_back(T(desRollOutput, 3, (180*m_ctrlParams.at(roll).kd*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(desRollOutput, 9, (180*lpf_b0*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(eiq, 4, (180*m_ctrlParams.at(pitch).kd)/M_PI));
    dhdx.push_back(T(edq, 10, (180*lpf_b0)/(M_PI*timestep)));
    dhdx.push_back(T(desPitchOutput, 4, (180*m_ctrlParams.at(pitch).kd*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(desPitchOutput, 10, (180*lpf_b0*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(eir, 5, (180*m_ctrlParams.at(yaw).kd)/M_PI));
    dhdx.push_back(T(edr, 11, 180/(M_PI*timestep)));
    dhdx.push_back(T(desYawOutput, 5, (180*m_ctrlParams.at(yaw).kd*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(desYawOutput, 11, (180*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(delay_1_rollRate, 9, 180/(M_PI*timestep)));
    dhdx.push_back(T(delay_1_pitchRate, 10, 180/(M_PI*timestep)));
    dhdx.push_back(T(w1, 2, (m_alpha*m_ctrlParams.at(posZ).kd*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep))/timestep));
    dhdx.push_back(T(w1, 3, -(90*m_alpha*m_ctrlParams.at(roll).kd*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 4, (90*m_alpha*m_ctrlParams.at(pitch).kd*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w1, 8, (m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale)/timestep));
    dhdx.push_back(T(w1, 9, -(90*lpf_b0*m_alpha*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w1, 10, (90*lpf_b0*m_alpha*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w1, 11, (180*m_alpha*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w2, 2, (m_alpha*m_ctrlParams.at(posZ).kd*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep))/timestep));
    dhdx.push_back(T(w2, 3, -(90*m_alpha*m_ctrlParams.at(roll).kd*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 4, -(90*m_alpha*m_ctrlParams.at(pitch).kd*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w2, 8, (m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale)/timestep));
    dhdx.push_back(T(w2, 9, -(90*lpf_b0*m_alpha*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w2, 10, -(90*lpf_b0*m_alpha*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w2, 11, -(180*m_alpha*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w3, 2, (m_alpha*m_ctrlParams.at(posZ).kd*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep))/timestep));
    dhdx.push_back(T(w3, 3, (90*m_alpha*m_ctrlParams.at(roll).kd*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 4, -(90*m_alpha*m_ctrlParams.at(pitch).kd*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w3, 8, (m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale)/timestep));
    dhdx.push_back(T(w3, 9, (90*lpf_b0*m_alpha*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w3, 10, -(90*lpf_b0*m_alpha*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w3, 11, (180*m_alpha*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w4, 2, (m_alpha*m_ctrlParams.at(posZ).kd*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep))/timestep));
    dhdx.push_back(T(w4, 3, (90*m_alpha*m_ctrlParams.at(roll).kd*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 4, (90*m_alpha*m_ctrlParams.at(pitch).kd*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(M_PI*timestep)));
    dhdx.push_back(T(w4, 8, (m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale)/timestep));
    dhdx.push_back(T(w4, 9, (90*lpf_b0*m_alpha*m_ctrlParams.at(rollRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w4, 10, (90*lpf_b0*m_alpha*m_ctrlParams.at(pitchRate).kd)/(M_PI*timestep)));
    dhdx.push_back(T(w4, 11, -(180*m_alpha*m_ctrlParams.at(yawRate).kd)/(M_PI*timestep)));

    if (curr.alge(desRoll) > -20 && curr.alge(desRoll) < 20)
    {
        dhdx.push_back(T(eiphi, 0, m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)));
        dhdx.push_back(T(eiphi, 1, -m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)));
        dhdx.push_back(T(eiphi, 5, timestep*((m_ctrlParams.at(velY).kd*(prev.plant(xdot)*cos(prev.plant(psi)) + prev.plant(ydot)*sin(prev.plant(psi))))/timestep + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi))) + (m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*(prev.plant(x)*cos(prev.plant(psi)) + prev.plant(y)*sin(prev.plant(psi))))/timestep)));
        dhdx.push_back(T(eiphi, 6, m_ctrlParams.at(velY).kd*sin(prev.plant(psi))));
        dhdx.push_back(T(eiphi, 7, -m_ctrlParams.at(velY).kd*cos(prev.plant(psi))));
        dhdx.push_back(T(desRollRate, 0, (m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep))/timestep));
        dhdx.push_back(T(desRollRate, 1, -(m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep))/timestep));
        dhdx.push_back(T(desRollRate, 5, ((m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd*prev.plant(xdot)*cos(prev.plant(psi)) + m_ctrlParams.at(velY).kd*prev.plant(ydot)*sin(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(x)*cos(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(y)*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(y)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(x)*timestep*cos(prev.plant(psi))))/timestep));
        dhdx.push_back(T(desRollRate, 6, (m_ctrlParams.at(velY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep))/timestep));
        dhdx.push_back(T(desRollRate, 7, -(m_ctrlParams.at(velY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep))/timestep));
        dhdx.push_back(T(eip, 0, m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(eip, 1, -m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(eip, 0, m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(eip, 1, -m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(eip, 5, (m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(velY).kd*prev.plant(xdot)*cos(prev.plant(psi)) + m_ctrlParams.at(velY).kd*prev.plant(ydot)*sin(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(x)*cos(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(y)*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(y)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(x)*timestep*cos(prev.plant(psi)))));
        dhdx.push_back(T(eip, 6, m_ctrlParams.at(velY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(eip, 7, -m_ctrlParams.at(velY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)));
        dhdx.push_back(T(desRollOutput, 0, (m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/timestep));
        dhdx.push_back(T(desRollOutput, 1, -(m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/timestep));
        dhdx.push_back(T(desRollOutput, 5, ((m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep)*(m_ctrlParams.at(velY).kd*prev.plant(xdot)*cos(prev.plant(psi)) + m_ctrlParams.at(velY).kd*prev.plant(ydot)*sin(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(x)*cos(prev.plant(psi)) + m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp*prev.plant(y)*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(y)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*prev.plant(x)*timestep*cos(prev.plant(psi))))/timestep));
        dhdx.push_back(T(desRollOutput, 6, (m_ctrlParams.at(velY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/timestep));
        dhdx.push_back(T(desRollOutput, 7, -(m_ctrlParams.at(velY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/timestep));
        dhdx.push_back(T(w1, 0, -(m_alpha*m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w1, 1, (m_alpha*m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 0, -(m_alpha*m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 1, (m_alpha*m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 0, (m_alpha*m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 1, -(m_alpha*m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 0, (m_alpha*m_ctrlParams.at(posY).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 1, -(m_alpha*m_ctrlParams.at(posY).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)*(m_ctrlParams.at(roll).kp + m_ctrlParams.at(roll).ki*timestep)*(m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep))/(2*timestep)));

    } 
    if (curr.alge(desPitch) > -20 && curr.alge(desPitch) < 20)
    {
        dhdx.push_back(T(eitheta, 0, m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)));
        dhdx.push_back(T(eitheta, 1, m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)));
        dhdx.push_back(T(eitheta, 5, timestep*((m_ctrlParams.at(velX).kd*(prev.plant(ydot)*cos(prev.plant(psi)) - prev.plant(xdot)*sin(prev.plant(psi))))/timestep + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi))) + (m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*(prev.plant(y)*cos(prev.plant(psi)) - prev.plant(x)*sin(prev.plant(psi))))/timestep)));
        dhdx.push_back(T(eitheta, 6, m_ctrlParams.at(velX).kd*cos(prev.plant(psi))));
        dhdx.push_back(T(eitheta, 7, m_ctrlParams.at(velX).kd*sin(prev.plant(psi))));
        dhdx.push_back(T(desPitchRate, 0, (m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchRate, 1, (m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchRate, 5, ((m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd*prev.plant(ydot)*cos(prev.plant(psi)) - m_ctrlParams.at(velX).kd*prev.plant(xdot)*sin(prev.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(y)*cos(prev.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(x)*sin(prev.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(x)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(y)*timestep*cos(prev.plant(psi))))/timestep));
        dhdx.push_back(T(desPitchRate, 6, (m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchRate, 7, (m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep))/timestep));
        dhdx.push_back(T(eiq, 0, m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)));
        dhdx.push_back(T(eiq, 1, m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)));
        dhdx.push_back(T(eiq, 5, (m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(velX).kd*prev.plant(ydot)*cos(prev.plant(psi)) - m_ctrlParams.at(velX).kd*prev.plant(xdot)*sin(prev.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(y)*cos(prev.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(x)*sin(prev.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(x)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(y)*timestep*cos(prev.plant(psi)))));
        dhdx.push_back(T(eiq, 6, m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)));
        dhdx.push_back(T(eiq, 7, m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)));
        dhdx.push_back(T(desPitchOutput, 0, (m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchOutput, 1, (m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchOutput, 5, ((m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep)*(m_ctrlParams.at(velX).kd*prev.plant(ydot)*cos(prev.plant(psi)) - m_ctrlParams.at(velX).kd*prev.plant(xdot)*sin(prev.plant(psi)) + m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(y)*cos(prev.plant(psi)) - m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp*prev.plant(x)*sin(prev.plant(psi)) - m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(x)*timestep*sin(prev.plant(psi)) + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*prev.plant(y)*timestep*cos(prev.plant(psi))))/timestep));
        dhdx.push_back(T(desPitchOutput, 6, (m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/timestep));
        dhdx.push_back(T(desPitchOutput, 7, (m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/timestep));
        dhdx.push_back(T(w1, 0, (m_alpha*m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w1, 1, (m_alpha*m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w1, 6, (m_alpha*m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w1, 7, (m_alpha*m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 0, -(m_alpha*m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 1, -(m_alpha*m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 6, -(m_alpha*m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w2, 7, -(m_alpha*m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 0, -(m_alpha*m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 1, -(m_alpha*m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 6, -(m_alpha*m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w3, 7, -(m_alpha*m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 0, (m_alpha*m_ctrlParams.at(posX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 1, (m_alpha*m_ctrlParams.at(posX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 6, (m_alpha*m_ctrlParams.at(velX).kd*cos(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
        dhdx.push_back(T(w4, 7, (m_alpha*m_ctrlParams.at(velX).kd*sin(prev.plant(psi))*(m_ctrlParams.at(pitch).kp + m_ctrlParams.at(pitch).ki*timestep)*(m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep))/(2*timestep)));
    }
    if (curr.alge(desRoll) == -20 || curr.alge(desRoll) == 20 || curr.alge(desPitch) == -20 || curr.alge(desPitch) == 20) 
    {
        dhdx.push_back(T(w1, 5, (180*m_alpha*m_ctrlParams.at(yaw).kd*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep))/(M_PI*timestep)));
        dhdx.push_back(T(w2, 5, -(180*m_alpha*m_ctrlParams.at(yaw).kd*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep))/(M_PI*timestep)));
        dhdx.push_back(T(w3, 5, (180*m_alpha*m_ctrlParams.at(yaw).kd*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep))/(M_PI*timestep)));
        dhdx.push_back(T(w4, 5, -(180*m_alpha*m_ctrlParams.at(yaw).kd*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep))/(M_PI*timestep)));
    }

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());

    if ((curr.alge(desRoll) > -20  && curr.alge(desRoll) < 20) || (curr.alge(desPitch) > -20  && curr.alge(desPitch) < 20))
    {
        dhdx.push_back(T(w1, 5, m_alpha*(-0.5*dhdx_mat.coeff(desRollOutput, 5) + 0.5*dhdx_mat.coeff(desPitchOutput, 5) + dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w2, 5, m_alpha*(-0.5*dhdx_mat.coeff(desRollOutput, 5) - 0.5*dhdx_mat.coeff(desPitchOutput, 5) - dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w3, 5, m_alpha*(+0.5*dhdx_mat.coeff(desRollOutput, 5) - 0.5*dhdx_mat.coeff(desPitchOutput, 5) + dhdx_mat.coeff(desYawOutput, 5))));
        dhdx.push_back(T(w4, 5, m_alpha*(+0.5*dhdx_mat.coeff(desRollOutput, 5) + 0.5*dhdx_mat.coeff(desPitchOutput, 5) - dhdx_mat.coeff(desYawOutput, 5))));
    }

    Eigen::Vector<double, NUM_PLANT_STATES> ft_vec = m_droneParams.kf * 2 * (curr.alge(w1)*dhdx_mat.row(w1) 
                                                                        + curr.alge(w2)*dhdx_mat.row(w2)
                                                                        + curr.alge(w3)*dhdx_mat.row(w3)
                                                                        + curr.alge(w4)*dhdx_mat.row(w4));
                                             
    Eigen::Vector<double, NUM_PLANT_STATES> tx_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (- curr.alge(w1)*dhdx_mat.row(w1) 
                                                                                                               - curr.alge(w2)*dhdx_mat.row(w2) 
                                                                                                               + curr.alge(w3)*dhdx_mat.row(w3) 
                                                                                                               + curr.alge(w4)*dhdx_mat.row(w4));
    Eigen::Vector<double, NUM_PLANT_STATES> ty_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (  curr.alge(w1)*dhdx_mat.row(w1) 
                                                                                                               - curr.alge(w2)*dhdx_mat.row(w2) 
                                                                                                               - curr.alge(w3)*dhdx_mat.row(w3) 
                                                                                                               + curr.alge(w4)*dhdx_mat.row(w4));
    Eigen::Vector<double, NUM_PLANT_STATES> tz_vec = m_droneParams.km * 2 * (curr.alge(w1)*dhdx_mat.row(w1) 
                                                                       - curr.alge(w2)*dhdx_mat.row(w2) 
                                                                       + curr.alge(w3)*dhdx_mat.row(w3) 
                                                                       - curr.alge(w4)*dhdx_mat.row(w4));

    // Append  entries to the original triplets
    for (int col = 0; col < NUM_PLANT_STATES; ++col) {
        if (ft_vec(col) != 0.0) {
            dhdx.emplace_back(ft, col, ft_vec(col));
        }
        if (tx_vec(col) != 0.0) {
            dhdx.emplace_back(tx, col, tx_vec(col));
        }
        if (ty_vec(col) != 0.0) {
            dhdx.emplace_back(ty, col, ty_vec(col));
        }
        if (tz_vec(col) != 0.0) {
            dhdx.emplace_back(tz, col, tz_vec(col));
        }
    }

    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzPlus(SystemState state, double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(300);
    // derivatives wrt eix
    dhdz.push_back(T(desVelX,  eix, m_ctrlParams.at(posX).ki));
    dhdz.push_back(T(eixdot, eix, m_ctrlParams.at(posX).ki*timestep));
    // derivatives wrt edx
    dhdz.push_back(T(desVelX, edx, m_ctrlParams.at(posX).kd));
    dhdz.push_back(T(eixdot, edx, m_ctrlParams.at(posX).kd*timestep));
    // derivatives wrt desVelX
    dhdz.push_back(T(eixdot, desVelX, timestep));
    // derivatives wrt eitheta
    dhdz.push_back(T(desPitchRate, eitheta, m_ctrlParams.at(pitch).ki));
    dhdz.push_back(T(eiq, eitheta, m_ctrlParams.at(pitch).ki*timestep));
    dhdz.push_back(T(desPitchOutput, eitheta, m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*m_ctrlParams.at(pitch).ki*timestep));
    dhdz.push_back(T(w1, eitheta, m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitchRate).ki*m_ctrlParams.at(pitch).ki*timestep)/2)));
    dhdz.push_back(T(w2, eitheta, -m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitchRate).ki*m_ctrlParams.at(pitch).ki*timestep)/2)));
    dhdz.push_back(T(w3, eitheta, -m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitchRate).ki*m_ctrlParams.at(pitch).ki*timestep)/2)));
    dhdz.push_back(T(w4, eitheta, m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitchRate).ki*m_ctrlParams.at(pitch).ki*timestep)/2)));
    // derivatives wrt edtheta
    dhdz.push_back(T(desPitchRate, edtheta, m_ctrlParams.at(pitch).kd));
    dhdz.push_back(T(eiq, edtheta, m_ctrlParams.at(pitch).kd*timestep));
    dhdz.push_back(T(desPitchOutput, edtheta, m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).ki*timestep));
    dhdz.push_back(T(w1, edtheta, m_alpha*((m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, edtheta, -m_alpha*((m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, edtheta, -m_alpha*((m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, edtheta, m_alpha*((m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).kd*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    // derivatives wrt desPitchRate
    dhdz.push_back(T(eiq, desPitchRate, timestep));
    dhdz.push_back(T(desPitchOutput, desPitchRate, m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitchRate).ki*timestep));
    dhdz.push_back(T(w1, desPitchRate, m_alpha*(m_ctrlParams.at(pitchRate).kp/2 + (m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, desPitchRate, -m_alpha*(m_ctrlParams.at(pitchRate).kp/2 + (m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, desPitchRate, -m_alpha*(m_ctrlParams.at(pitchRate).kp/2 + (m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, desPitchRate, m_alpha*(m_ctrlParams.at(pitchRate).kp/2 + (m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    // derivatives wrt delay_1_pitchRate
    dhdz.push_back(T(edq, delay_1_pitchRate, lpf_b0));
    dhdz.push_back(T(desPitchOutput, delay_1_pitchRate, m_ctrlParams.at(pitchRate).kd*lpf_b0));
    dhdz.push_back(T(w1, delay_1_pitchRate, (m_ctrlParams.at(pitchRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w2, delay_1_pitchRate, -(m_ctrlParams.at(pitchRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w3, delay_1_pitchRate, -(m_ctrlParams.at(pitchRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w4, delay_1_pitchRate, (m_ctrlParams.at(pitchRate).kd*lpf_b0*m_alpha)/2));
    // derivatives wrt delay_2_pitchRate
    // derivatives wrt eiq
    dhdz.push_back(T(desPitchOutput, eiq, m_ctrlParams.at(pitchRate).ki));
    dhdz.push_back(T(w1, eiq, (m_ctrlParams.at(pitchRate).ki*m_alpha)/2));
    dhdz.push_back(T(w2, eiq, -(m_ctrlParams.at(pitchRate).ki*m_alpha)/2));
    dhdz.push_back(T(w3, eiq, -(m_ctrlParams.at(pitchRate).ki*m_alpha)/2));
    dhdz.push_back(T(w4, eiq, (m_ctrlParams.at(pitchRate).ki*m_alpha)/2));
    // derivatives wrt edq
    dhdz.push_back(T(desPitchOutput, edq, m_ctrlParams.at(pitchRate).kd));
    dhdz.push_back(T(w1, edq, (m_ctrlParams.at(pitchRate).kd*m_alpha)/2));
    dhdz.push_back(T(w2, edq, -(m_ctrlParams.at(pitchRate).kd*m_alpha)/2));
    dhdz.push_back(T(w3, edq, -(m_ctrlParams.at(pitchRate).kd*m_alpha)/2));
    dhdz.push_back(T(w4, edq, (m_ctrlParams.at(pitchRate).kd*m_alpha)/2));
    // derivatives wrt desPitchOutput
    dhdz.push_back(T(w1, desPitchOutput, m_alpha/2));
    dhdz.push_back(T(w2, desPitchOutput, -m_alpha/2));
    dhdz.push_back(T(w3, desPitchOutput, -m_alpha/2));
    dhdz.push_back(T(w4, desPitchOutput, m_alpha/2));
    
    // derivatives wrt eiy
    dhdz.push_back(T(desVelY, eiy, m_ctrlParams.at(posY).ki));
    dhdz.push_back(T(eiydot, eiy, m_ctrlParams.at(posY).ki*timestep));
    // derivatives wrt edy
    dhdz.push_back(T(desVelY, edy, m_ctrlParams.at(posY).kd));
    dhdz.push_back(T(eiydot, edy, m_ctrlParams.at(posY).kd*timestep));
    // derivatives wrt desVelY
    dhdz.push_back(T(eiydot, desVelY, timestep));
    // derivatives wrt eiphi
    dhdz.push_back(T(desRollRate, eiphi, m_ctrlParams.at(roll).ki));
    dhdz.push_back(T(eip, eiphi, m_ctrlParams.at(roll).ki*timestep));
    dhdz.push_back(T(desRollOutput, eiphi, m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*m_ctrlParams.at(roll).ki*timestep));
    dhdz.push_back(T(w1, eiphi, -m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(rollRate).ki*m_ctrlParams.at(roll).ki*timestep)/2)));
    dhdz.push_back(T(w2, eiphi, -m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(rollRate).ki*m_ctrlParams.at(roll).ki*timestep)/2)));
    dhdz.push_back(T(w3, eiphi, m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(rollRate).ki*m_ctrlParams.at(roll).ki*timestep)/2)));
    dhdz.push_back(T(w4, eiphi, m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(rollRate).ki*m_ctrlParams.at(roll).ki*timestep)/2)));
    // derivatives wrt edphi
    dhdz.push_back(T(desRollRate, edphi, m_ctrlParams.at(roll).kd));
    dhdz.push_back(T(eip, edphi, m_ctrlParams.at(roll).kd*timestep));
    dhdz.push_back(T(desRollOutput, edphi, m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).ki*timestep));
    dhdz.push_back(T(w1, edphi, -m_alpha*((m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, edphi, -m_alpha*((m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, edphi, m_alpha*((m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, edphi, m_alpha*((m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).kd*m_ctrlParams.at(rollRate).ki*timestep)/2))); 
    // derivatives wrt desRollRate
    dhdz.push_back(T(eip, desRollRate, timestep));
    dhdz.push_back(T(desRollOutput, desRollRate, m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(rollRate).ki*timestep));
    dhdz.push_back(T(w1, desRollRate, -m_alpha*(m_ctrlParams.at(rollRate).kp/2 + (m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, desRollRate, -m_alpha*(m_ctrlParams.at(rollRate).kp/2 + (m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, desRollRate, m_alpha*(m_ctrlParams.at(rollRate).kp/2 + (m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, desRollRate, m_alpha*(m_ctrlParams.at(rollRate).kp/2 + (m_ctrlParams.at(rollRate).ki*timestep)/2)));
    // derivatives wrt delay_1_rollRate
    dhdz.push_back(T(edp, delay_1_rollRate, lpf_b0));
    dhdz.push_back(T(desRollOutput, delay_1_rollRate, m_ctrlParams.at(rollRate).kd*lpf_b0));
    dhdz.push_back(T(w1, delay_1_rollRate, -(m_ctrlParams.at(rollRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w2, delay_1_rollRate, -(m_ctrlParams.at(rollRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w3, delay_1_rollRate, (m_ctrlParams.at(rollRate).kd*lpf_b0*m_alpha)/2));
    dhdz.push_back(T(w4, delay_1_rollRate, (m_ctrlParams.at(rollRate).kd*lpf_b0*m_alpha)/2));
    // derivatives wrt delay_2_rollRate
    // derivatives wrt eip
    dhdz.push_back(T(desRollOutput, eip, m_ctrlParams.at(rollRate).ki));
    dhdz.push_back(T(w1, eip, -(m_ctrlParams.at(rollRate).ki*m_alpha)/2));
    dhdz.push_back(T(w2, eip, -(m_ctrlParams.at(rollRate).ki*m_alpha)/2));
    dhdz.push_back(T(w3, eip, (m_ctrlParams.at(rollRate).ki*m_alpha)/2));
    dhdz.push_back(T(w4, eip, (m_ctrlParams.at(rollRate).ki*m_alpha)/2));
    // derivatives wrt edp
    dhdz.push_back(T(desRollOutput, edp, m_ctrlParams.at(rollRate).kd));
    dhdz.push_back(T(w1, edp, -(m_ctrlParams.at(rollRate).kd*m_alpha)/2));
    dhdz.push_back(T(w2, edp, -(m_ctrlParams.at(rollRate).kd*m_alpha)/2));
    dhdz.push_back(T(w3, edp, (m_ctrlParams.at(rollRate).kd*m_alpha)/2));
    dhdz.push_back(T(w4, edp, (m_ctrlParams.at(rollRate).kd*m_alpha)/2));
    // derivatives wrt desRollOutput
    dhdz.push_back(T(w1, desRollOutput, -m_alpha/2));
    dhdz.push_back(T(w2, desRollOutput, -m_alpha/2));
    dhdz.push_back(T(w3, desRollOutput, m_alpha/2));
    dhdz.push_back(T(w4, desRollOutput, m_alpha/2));

    // derivatives wrt eiz
    dhdz.push_back(T(desVelZ, eiz, m_ctrlParams.at(posZ).ki));
    dhdz.push_back(T(eizdot, eiz, m_ctrlParams.at(posZ).ki*timestep));
    dhdz.push_back(T(desThrust, eiz, m_thrustScale*(m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w1, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w2, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w3, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w4, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).ki*m_ctrlParams.at(velZ).ki*timestep)));  
    // derivatives wrt edz
    dhdz.push_back(T(desVelZ,  edz, m_ctrlParams.at(posZ).kd));
    dhdz.push_back(T(eizdot, edz, m_ctrlParams.at(posZ).kd*timestep));
    dhdz.push_back(T(desThrust, edz, m_thrustScale*(m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w1, edz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w2, edz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w3, edz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w4, edz, m_alpha*m_thrustScale*(m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).kp + m_ctrlParams.at(posZ).kd*m_ctrlParams.at(velZ).ki*timestep))); 
    // derivatives wrt desVelZ
    dhdz.push_back(T(eizdot, desVelZ, timestep));
    dhdz.push_back(T(desThrust, desVelZ, m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w1, desVelZ, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w2, desVelZ, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w3, desVelZ, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)));
    dhdz.push_back(T(w4, desVelZ, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp + m_ctrlParams.at(velZ).ki*timestep)));
    // derivatives wrt to eizdot
    dhdz.push_back(T(desThrust, eizdot, m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(w1, eizdot, m_ctrlParams.at(velZ).ki*m_alpha*m_thrustScale));
    dhdz.push_back(T(w2, eizdot, m_ctrlParams.at(velZ).ki*m_alpha*m_thrustScale));
    dhdz.push_back(T(w3, eizdot, m_ctrlParams.at(velZ).ki*m_alpha*m_thrustScale));
    dhdz.push_back(T(w4, eizdot, m_ctrlParams.at(velZ).ki*m_alpha*m_thrustScale));
    // derivatives wrt edzdot
    dhdz.push_back(T(desThrust, edzdot, m_ctrlParams.at(velZ).kd*m_thrustScale));
    dhdz.push_back(T(w1, edzdot, m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale));
    dhdz.push_back(T(w2, edzdot, m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale));
    dhdz.push_back(T(w3, edzdot, m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale));
    dhdz.push_back(T(w4, edzdot, m_ctrlParams.at(velZ).kd*m_alpha*m_thrustScale));
    // derivatives wrt desThrust
    dhdz.push_back(T(w1, desThrust, m_alpha));
    dhdz.push_back(T(w2, desThrust, m_alpha));
    dhdz.push_back(T(w3, desThrust, m_alpha));
    dhdz.push_back(T(w4, desThrust, m_alpha));
    // derivatives wrt eipsi
    dhdz.push_back(T(desYawRate, eipsi, m_ctrlParams.at(yaw).ki));
    dhdz.push_back(T(eir, eipsi, m_ctrlParams.at(yaw).ki*timestep));
    dhdz.push_back(T(desYawOutput, eipsi, m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep));
    dhdz.push_back(T(w1, eipsi, m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w2, eipsi, -m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w3, eipsi, m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w4, eipsi, -m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    // derivatives wrt edpsi
    dhdz.push_back(T(desYawRate, edpsi, m_ctrlParams.at(yaw).kd));
    dhdz.push_back(T(eir, edpsi, m_ctrlParams.at(yaw).kd*timestep));
    dhdz.push_back(T(desYawOutput, edpsi, m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).ki*timestep));
    dhdz.push_back(T(w1, edpsi, m_alpha*(m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w2, edpsi, -m_alpha*(m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w3, edpsi, m_alpha*(m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w4, edpsi, -m_alpha*(m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).kd*m_ctrlParams.at(yawRate).ki*timestep)));
    // derivatives wrt desYawRate
    dhdz.push_back(T(eir, desYawRate, timestep));
    dhdz.push_back(T(desYawOutput, desYawRate, m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep));
    dhdz.push_back(T(w1, desYawRate, m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w2, desYawRate, -m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w3, desYawRate, m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w4, desYawRate, -m_alpha*(m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yawRate).ki*timestep)));
    // eir 
    dhdz.push_back(T(desYawOutput, eir, m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(w1, eir, m_ctrlParams.at(yawRate).ki*m_alpha));
    dhdz.push_back(T(w2, eir, -m_ctrlParams.at(yawRate).ki*m_alpha));
    dhdz.push_back(T(w3, eir, m_ctrlParams.at(yawRate).ki*m_alpha));
    dhdz.push_back(T(w4, eir, -m_ctrlParams.at(yawRate).ki*m_alpha));
    // derivatives wrt edr
    dhdz.push_back(T(desYawOutput, edr, m_ctrlParams.at(yawRate).kd));
    dhdz.push_back(T(w1, edr, m_ctrlParams.at(yawRate).kd*m_alpha));
    dhdz.push_back(T(w2, edr, -m_ctrlParams.at(yawRate).kd*m_alpha));
    dhdz.push_back(T(w3, edr, m_ctrlParams.at(yawRate).kd*m_alpha));
    dhdz.push_back(T(w4, edr, -m_ctrlParams.at(yawRate).kd*m_alpha));
    // derivatives wrt desYawOutput
    dhdz.push_back(T(w1, desYawOutput, m_alpha));
    dhdz.push_back(T(w2, desYawOutput, -m_alpha));
    dhdz.push_back(T(w3, desYawOutput, m_alpha));
    dhdz.push_back(T(w4, desYawOutput, -m_alpha));

    if (state.alge(desPitch) > -20 && state.alge(desPitch) < 20) {
        m_logger << "desPitch not saturated " << std:: endl;
        // derivatives wrt eix
        dhdz.push_back(T(eitheta, eix, timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)));
        dhdz.push_back(T(desPitchRate, eix, m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)));
        dhdz.push_back(T(eiq, eix, timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep))));
        dhdz.push_back(T(desPitchOutput, eix, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep))));
        dhdz.push_back(T(w1, eix, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w2, eix, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w3, eix, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w4, eix, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        // derivatives wrt edx
        dhdz.push_back(T(eitheta, edx, timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)));
        dhdz.push_back(T(desPitchRate, edx, m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)));
        dhdz.push_back(T(eiq, edx, timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep))));
        dhdz.push_back(T(desPitchOutput, edx, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep))));
        dhdz.push_back(T(w1, edx, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2)));
        dhdz.push_back(T(w2, edx, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2)));
        dhdz.push_back(T(w3, edx, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2)));
        dhdz.push_back(T(w4, edx, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(posX).kd*m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).kd*timestep)))/2)));
        // derivatives wrt desVelX
        dhdz.push_back(T(eitheta, desVelX, timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)));
        dhdz.push_back(T(desPitchRate, desVelX, m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)));
        dhdz.push_back(T(eiq, desVelX, timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep))));
        dhdz.push_back(T(desPitchOutput, desVelX, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep))));
        dhdz.push_back(T(w1, desVelX, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2)));
        dhdz.push_back(T(w2, desVelX, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2)));
        dhdz.push_back(T(w3, desVelX, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2)));
        dhdz.push_back(T(w4, desVelX, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp + m_ctrlParams.at(velX).ki*timestep)))/2)));
        // derivatives wrt eixdot
        dhdz.push_back(T(eitheta, eixdot, m_ctrlParams.at(velX).ki*timestep));
        dhdz.push_back(T(desPitchRate, eixdot, m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep));
        dhdz.push_back(T(eiq, eixdot, timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(desPitchOutput, eixdot, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(w1, eixdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w2, eixdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w3, eixdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w4, eixdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        // derivatives wrt edxdot
        dhdz.push_back(T(eitheta, edxdot, m_ctrlParams.at(velX).kd*timestep));
        dhdz.push_back(T(desPitchRate, edxdot, m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep));
        dhdz.push_back(T(eiq, edxdot, timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(desPitchOutput, edxdot, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(w1, edxdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w2, edxdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w3, edxdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w4, edxdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).kd*m_ctrlParams.at(pitch).ki*timestep))/2)));
    }
    if (state.alge(desRoll) > -20 && state.alge(desRoll) < 20) {
        m_logger << "desRoll not saturated " << std:: endl;
        // derivatives wrt eiy
        dhdz.push_back(T(eiphi, eiy, -timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)));
        dhdz.push_back(T(desRollRate, eiy, - m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) - m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)));
        dhdz.push_back(T(eip, eiy, -timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep))));
        dhdz.push_back(T(desRollOutput, eiy, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep))));
        dhdz.push_back(T(w1, eiy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w2, eiy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w3, eiy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w4, eiy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        // derivatives wrt edy
        dhdz.push_back(T(eiphi, edy, -timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)));
        dhdz.push_back(T(desRollRate, edy, - m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) - m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)));
        dhdz.push_back(T(eip, edy, -timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep))));
        dhdz.push_back(T(desRollOutput, edy, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep))));
        dhdz.push_back(T(w1, edy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2)));
        dhdz.push_back(T(w2, edy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2)));
        dhdz.push_back(T(w3, edy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2)));
        dhdz.push_back(T(w4, edy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(posY).kd*m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).kd*timestep)))/2)));
        // derivatives wrt desVelY
        dhdz.push_back(T(eiphi, desVelY, -timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)));
        dhdz.push_back(T(desRollRate, desVelY, - m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) - m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)));
        dhdz.push_back(T(eip, desVelY, -timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep))));
        dhdz.push_back(T(desRollOutput, desVelY, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep))));
        dhdz.push_back(T(w1, desVelY, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2)));
        dhdz.push_back(T(w2, desVelY, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2)));
        dhdz.push_back(T(w3, desVelY, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2)));
        dhdz.push_back(T(w4, desVelY, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp + m_ctrlParams.at(velY).ki*timestep)))/2)));
        // derivatives wrt edydot
        dhdz.push_back(T(eiphi, eiydot, -m_ctrlParams.at(velY).ki*timestep));
        dhdz.push_back(T(desRollRate, eiydot, - m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp - m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep));
        dhdz.push_back(T(eip, eiydot, -timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(desRollOutput, eiydot, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(w1, eiydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w2, eiydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w3, eiydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w4, eiydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        // derivatives wrt edydot
        dhdz.push_back(T(eiphi, edydot, -m_ctrlParams.at(velY).kd*timestep));
        dhdz.push_back(T(desRollRate, edydot, - m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp - m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep));
        dhdz.push_back(T(eip, edydot, -timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(desRollOutput, edydot, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(w1, edydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w2, edydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w3, edydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w4, edydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).kd*m_ctrlParams.at(roll).ki*timestep))/2)));
    }

    // w1 ... w4
    dhdz.push_back(T(ft, w1, m_droneParams.kf * 2 * state.alge(w1)));
    dhdz.push_back(T(tx, w1, - m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w1)));
    dhdz.push_back(T(ty, w1, m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w1)));
    dhdz.push_back(T(tz, w1, m_droneParams.km * 2 * state.alge(w1)));

    dhdz.push_back(T(ft, w2, m_droneParams.kf * 2 * state.alge(w2)));
    dhdz.push_back(T(tx, w2, - m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w2)));
    dhdz.push_back(T(ty, w2, - m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w2)));
    dhdz.push_back(T(tz, w2, - m_droneParams.km * 2 * state.alge(w2)));

    dhdz.push_back(T(ft, w3, m_droneParams.kf * 2 * state.alge(w3)));
    dhdz.push_back(T(tx, w3, m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w3)));
    dhdz.push_back(T(ty, w3, - m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w3)));
    dhdz.push_back(T(tz, w3, - m_droneParams.km * 2 * state.alge(w3)));

    dhdz.push_back(T(ft, w4, m_droneParams.kf * 2 * state.alge(w4)));
    dhdz.push_back(T(tx, w4, m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w4)));
    dhdz.push_back(T(ty, w4, m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * state.alge(w4)));
    dhdz.push_back(T(tz, w4, - m_droneParams.km * 2 * state.alge(w4)));
    
    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());

    Eigen::Vector<double, NUM_Z_STATES> ft_vec = m_droneParams.kf * 2 * (state.alge(w1)*dhdz_mat.row(w1) 
                                                                       + state.alge(w2)*dhdz_mat.row(w2)
                                                                       + state.alge(w3)*dhdz_mat.row(w3)
                                                                       + state.alge(w4)*dhdz_mat.row(w4));
                                             
    Eigen::Vector<double, NUM_Z_STATES> tx_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (- state.alge(w1)*dhdz_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdz_mat.row(w2) 
                                                                                                               + state.alge(w3)*dhdz_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdz_mat.row(w4));
    Eigen::Vector<double, NUM_Z_STATES> ty_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (  state.alge(w1)*dhdz_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdz_mat.row(w2) 
                                                                                                               - state.alge(w3)*dhdz_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdz_mat.row(w4));
    Eigen::Vector<double, NUM_Z_STATES> tz_vec = m_droneParams.km * 2 * (state.alge(w1)*dhdz_mat.row(w1) 
                                                                       - state.alge(w2)*dhdz_mat.row(w2) 
                                                                       + state.alge(w3)*dhdz_mat.row(w3) 
                                                                       - state.alge(w4)*dhdz_mat.row(w4));

    // Append  entries to the original triplets
    for (int col = 0; col < NUM_Z_STATES; ++col) {
        if (ft_vec(col) != 0.0) {
            dhdz.emplace_back(ft, col, ft_vec(col));
        }
        if (tx_vec(col) != 0.0) {
            dhdz.emplace_back(tx, col, tx_vec(col));
        }
        if (ty_vec(col) != 0.0) {
            dhdz.emplace_back(ty, col, ty_vec(col));
        }
        if (tz_vec(col) != 0.0) {
            dhdz.emplace_back(tz, col, tz_vec(col));
        }
        // also identity matrix
        dhdz.emplace_back(col, col, 1);
    }

    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzCurr(SystemState state, double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(93);
    dhdz.push_back(T(eix, eix, 1));
    dhdz.push_back(T(desVelX, eix, m_ctrlParams.at(posX).ki));
    dhdz.push_back(T(eiy, eiy, 1));
    dhdz.push_back(T(desVelY, eiy, m_ctrlParams.at(posY).ki));
    dhdz.push_back(T(eiz, eiz, 1));
    dhdz.push_back(T(desVelZ, eiz, m_ctrlParams.at(posZ).ki));
    dhdz.push_back(T(eixdot, eix, m_ctrlParams.at(posX).ki*timestep));
    dhdz.push_back(T(eixdot, eixdot, 1));
    dhdz.push_back(T(eiydot, eiy, m_ctrlParams.at(posY).ki*timestep));
    dhdz.push_back(T(eiydot, eiydot, 1));
    dhdz.push_back(T(eizdot, eiz, m_ctrlParams.at(posZ).ki*timestep));
    dhdz.push_back(T(eizdot, eizdot, 1));
    dhdz.push_back(T(desThrust, eiz, m_thrustScale*(m_ctrlParams.at(velZ).kp*m_ctrlParams.at(posZ).ki + m_ctrlParams.at(velZ).ki*m_ctrlParams.at(posZ).ki*timestep)));
    dhdz.push_back(T(desThrust, eizdot, m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(eiphi, eiphi, 1));
    dhdz.push_back(T(desRollRate, eiphi, m_ctrlParams.at(roll).ki));
    dhdz.push_back(T(eitheta, eitheta, 1));
    dhdz.push_back(T(desPitchRate, eitheta, m_ctrlParams.at(pitch).ki));
    dhdz.push_back(T(eipsi, eipsi, 1));
    dhdz.push_back(T(desYawRate, eipsi, m_ctrlParams.at(yaw).ki));
    dhdz.push_back(T(eip, eiphi, m_ctrlParams.at(roll).ki*timestep));
    dhdz.push_back(T(eip, eip, 1));
    dhdz.push_back(T(edp, delay_1_rollRate, lpf_b1 - lpf_a1*lpf_b0));
    dhdz.push_back(T(edp, delay_2_rollRate, lpf_b2 - lpf_b0*lpf_a2));
    dhdz.push_back(T(desRollOutput, eiphi, m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp + m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).ki*timestep));
    dhdz.push_back(T(desRollOutput, eip, m_ctrlParams.at(rollRate).ki));
    dhdz.push_back(T(desRollOutput, delay_1_rollRate, m_ctrlParams.at(rollRate).kd*(lpf_b1 - lpf_a1*lpf_b0)));
    dhdz.push_back(T(desRollOutput, delay_2_rollRate, m_ctrlParams.at(rollRate).kd*(lpf_b2 - lpf_b0*lpf_a2)));
    dhdz.push_back(T(eiq, eitheta, m_ctrlParams.at(pitch).ki*timestep));
    dhdz.push_back(T(eiq, eiq, 1));
    dhdz.push_back(T(edq, delay_1_pitchRate, lpf_b1 - lpf_a1*lpf_b0));
    dhdz.push_back(T(edq, delay_2_pitchRate, lpf_b2 - lpf_b0*lpf_a2));
    dhdz.push_back(T(desPitchOutput, eitheta, m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp + m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).ki*timestep));
    dhdz.push_back(T(desPitchOutput, eiq, m_ctrlParams.at(pitchRate).ki));
    dhdz.push_back(T(desPitchOutput, delay_1_pitchRate, m_ctrlParams.at(pitchRate).kd*(lpf_b1 - lpf_a1*lpf_b0)));
    dhdz.push_back(T(desPitchOutput, delay_2_pitchRate, m_ctrlParams.at(pitchRate).kd*(lpf_b2 - lpf_b0*lpf_a2)));
    dhdz.push_back(T(eir, eipsi, m_ctrlParams.at(yaw).ki*timestep));
    dhdz.push_back(T(eir, eir, 1));
    dhdz.push_back(T(desYawOutput, eipsi, m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep));
    dhdz.push_back(T(desYawOutput, eir, m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(delay_1_rollRate, delay_1_rollRate, -lpf_a1));
    dhdz.push_back(T(delay_1_rollRate, delay_2_rollRate, -lpf_a2));
    dhdz.push_back(T(delay_2_rollRate, delay_1_rollRate, 1));
    dhdz.push_back(T(delay_1_pitchRate, delay_1_pitchRate, -lpf_a1));
    dhdz.push_back(T(delay_1_pitchRate, delay_2_pitchRate, -lpf_a2));
    dhdz.push_back(T(delay_2_pitchRate, delay_1_pitchRate, 1));
    dhdz.push_back(T(w1, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp*m_ctrlParams.at(posZ).ki + m_ctrlParams.at(velZ).ki*m_ctrlParams.at(posZ).ki*timestep)));
    dhdz.push_back(T(w1, eizdot, m_alpha*m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(w1, eiphi, -m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w1, eitheta, m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w1, eipsi, m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w1, eip, -(m_alpha*m_ctrlParams.at(rollRate).ki)/2));
    dhdz.push_back(T(w1, eiq, (m_alpha*m_ctrlParams.at(pitchRate).ki)/2));
    dhdz.push_back(T(w1, eir, m_alpha*m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(w1, delay_1_rollRate, -(m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w1, delay_2_rollRate, -(m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w1, delay_1_pitchRate, (m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w1, delay_2_pitchRate, (m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w2, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp*m_ctrlParams.at(posZ).ki + m_ctrlParams.at(velZ).ki*m_ctrlParams.at(posZ).ki*timestep)));
    dhdz.push_back(T(w2, eizdot, m_alpha*m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(w2, eiphi, -m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, eitheta, -m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w2, eipsi, -m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w2, eip, -(m_alpha*m_ctrlParams.at(rollRate).ki)/2));
    dhdz.push_back(T(w2, eiq, -(m_alpha*m_ctrlParams.at(pitchRate).ki)/2));
    dhdz.push_back(T(w2, eir, -m_alpha*m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(w2, delay_1_rollRate, -(m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w2, delay_2_rollRate, -(m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w2, delay_1_pitchRate, -(m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w2, delay_2_pitchRate, -(m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w3, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp*m_ctrlParams.at(posZ).ki + m_ctrlParams.at(velZ).ki*m_ctrlParams.at(posZ).ki*timestep)));
    dhdz.push_back(T(w3, eizdot, m_alpha*m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(w3, eiphi, m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, eitheta, -m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w3, eipsi, m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w3, eip, (m_alpha*m_ctrlParams.at(rollRate).ki)/2));
    dhdz.push_back(T(w3, eiq, -(m_alpha*m_ctrlParams.at(pitchRate).ki)/2));
    dhdz.push_back(T(w3, eir, m_alpha*m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(w3, delay_1_rollRate, (m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w3, delay_2_rollRate, (m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w3, delay_1_pitchRate, -(m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w3, delay_2_pitchRate, -(m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w4, eiz, m_alpha*m_thrustScale*(m_ctrlParams.at(velZ).kp*m_ctrlParams.at(posZ).ki + m_ctrlParams.at(velZ).ki*m_ctrlParams.at(posZ).ki*timestep)));
    dhdz.push_back(T(w4, eizdot, m_alpha*m_ctrlParams.at(velZ).ki*m_thrustScale));
    dhdz.push_back(T(w4, eiphi, m_alpha*((m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).kp)/2 + (m_ctrlParams.at(roll).ki*m_ctrlParams.at(rollRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, eitheta, m_alpha*((m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).kp)/2 + (m_ctrlParams.at(pitch).ki*m_ctrlParams.at(pitchRate).ki*timestep)/2)));
    dhdz.push_back(T(w4, eipsi, -m_alpha*(m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).kp + m_ctrlParams.at(yaw).ki*m_ctrlParams.at(yawRate).ki*timestep)));
    dhdz.push_back(T(w4, eip, (m_alpha*m_ctrlParams.at(rollRate).ki)/2));
    dhdz.push_back(T(w4, eiq, (m_alpha*m_ctrlParams.at(pitchRate).ki)/2));
    dhdz.push_back(T(w4, eir, -m_alpha*m_ctrlParams.at(yawRate).ki));
    dhdz.push_back(T(w4, delay_1_rollRate, (m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w4, delay_2_rollRate, (m_alpha*m_ctrlParams.at(rollRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));
    dhdz.push_back(T(w4, delay_1_pitchRate, (m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b1 - lpf_a1*lpf_b0))/2));
    dhdz.push_back(T(w4, delay_2_pitchRate, (m_alpha*m_ctrlParams.at(pitchRate).kd*(lpf_b2 - lpf_b0*lpf_a2))/2));

    if (state.alge(desRoll) > -20 && state.alge(desRoll) < 20) {
        dhdz.push_back(T(eiphi, eiy, -timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)));
        dhdz.push_back(T(eiphi, eiydot, -m_ctrlParams.at(velY).ki*timestep));
        dhdz.push_back(T(desRollRate, eiy, - m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) - m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)));
        dhdz.push_back(T(desRollRate, eiydot, - m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp - m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep));
        dhdz.push_back(T(eip, eiy, -timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep))));
        dhdz.push_back(T(eip, eiydot, -timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(desRollOutput, eiy, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep))));
        dhdz.push_back(T(desRollOutput, eiydot, - m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep) - m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep)));
        dhdz.push_back(T(w1, eiy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w1, eiydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w2, eiy, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w2, eiydot, m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w3, eiy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w3, eiydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
        dhdz.push_back(T(w4, eiy, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(roll).kp*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep) + m_ctrlParams.at(roll).ki*timestep*(m_ctrlParams.at(velY).kp*m_ctrlParams.at(posY).ki + m_ctrlParams.at(velY).ki*m_ctrlParams.at(posY).ki*timestep)))/2)));
        dhdz.push_back(T(w4, eiydot, -m_alpha*((m_ctrlParams.at(rollRate).kp*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2 + (m_ctrlParams.at(rollRate).ki*timestep*(m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).kp + m_ctrlParams.at(velY).ki*m_ctrlParams.at(roll).ki*timestep))/2)));
    }
    if (state.alge(desPitch) > -20 && state.alge(desPitch) < 20)
    {
        dhdz.push_back(T(eitheta, eix, timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)));
        dhdz.push_back(T(eitheta, eixdot, m_ctrlParams.at(velX).ki*timestep));
        dhdz.push_back(T(desPitchRate, eix, m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)));
        dhdz.push_back(T(desPitchRate, eixdot, m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep));
        dhdz.push_back(T(eiq, eix, timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep))));
        dhdz.push_back(T(eiq, eixdot, timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(desPitchOutput, eix, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep))));
        dhdz.push_back(T(desPitchOutput, eixdot, m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep) + m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep)));
        dhdz.push_back(T(w1, eix, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w1, eixdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w2, eix, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w2, eixdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w3, eix, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w3, eixdot, -m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
        dhdz.push_back(T(w4, eix, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(pitch).kp*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep) + m_ctrlParams.at(pitch).ki*timestep*(m_ctrlParams.at(velX).kp*m_ctrlParams.at(posX).ki + m_ctrlParams.at(velX).ki*m_ctrlParams.at(posX).ki*timestep)))/2)));
        dhdz.push_back(T(w4, eixdot, m_alpha*((m_ctrlParams.at(pitchRate).kp*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2 + (m_ctrlParams.at(pitchRate).ki*timestep*(m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).kp + m_ctrlParams.at(velX).ki*m_ctrlParams.at(pitch).ki*timestep))/2)));
    }
    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());

    Eigen::Vector<double, NUM_Z_STATES> ft_vec = m_droneParams.kf * 2 * (state.alge(w1)*dhdz_mat.row(w1) 
                                                                       + state.alge(w2)*dhdz_mat.row(w2)
                                                                       + state.alge(w3)*dhdz_mat.row(w3)
                                                                       + state.alge(w4)*dhdz_mat.row(w4));
                                             
    Eigen::Vector<double, NUM_Z_STATES> tx_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (- state.alge(w1)*dhdz_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdz_mat.row(w2) 
                                                                                                               + state.alge(w3)*dhdz_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdz_mat.row(w4));
    Eigen::Vector<double, NUM_Z_STATES> ty_vec = m_droneParams.kf * m_droneParams.length * 2.0/std::sqrt(2) * (  state.alge(w1)*dhdz_mat.row(w1) 
                                                                                                               - state.alge(w2)*dhdz_mat.row(w2) 
                                                                                                               - state.alge(w3)*dhdz_mat.row(w3) 
                                                                                                               + state.alge(w4)*dhdz_mat.row(w4));
    Eigen::Vector<double, NUM_Z_STATES> tz_vec = m_droneParams.km * 2 * (state.alge(w1)*dhdz_mat.row(w1) 
                                                                       - state.alge(w2)*dhdz_mat.row(w2) 
                                                                       + state.alge(w3)*dhdz_mat.row(w3) 
                                                                       - state.alge(w4)*dhdz_mat.row(w4));

    // Append  entries to the original triplets
    for (int col = 0; col < NUM_Z_STATES; ++col) {
        if (ft_vec(col) != 0.0) {
            dhdz.emplace_back(ft, col, ft_vec(col));
        }
        if (tx_vec(col) != 0.0) {
            dhdz.emplace_back(tx, col, tx_vec(col));
        }
        if (ty_vec(col) != 0.0) {
            dhdz.emplace_back(ty, col, ty_vec(col));
        }
        if (tz_vec(col) != 0.0) {
            dhdz.emplace_back(tz, col, tz_vec(col));
        }
    }

    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdy()
{
    Eigen::SparseMatrix<double> dhdy(NUM_Z_STATES, NUM_Y_STATES);
    return dhdy; 
}

