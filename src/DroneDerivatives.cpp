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

Eigen::SparseMatrix<double> DroneTrajectory::dgdz(SystemState state)
{
    std::vector<T> dgdz;
    // if(state.alge(desRoll) > -20 && state.alge(desRoll) < 20){
    //     dgdz.reserve(3);
    //     dgdz.push_back(T(1, epydot, m_ctrlParams.at(velY).kp));
         dgdz.push_back(T(1, eiydot, m_ctrlParams.at(velY).ki*state.alge(desRoll)));
    //     dgdz.push_back(T(1, edydot, m_ctrlParams.at(velY).kd));
    // }
    // if(state.alge(desPitch) > -20 && state.alge(desPitch) < 20){
    //     dgdz.reserve(3);
    //     dgdz.push_back(T(1, epxdot, m_ctrlParams.at(velX).kp));
    //     dgdz.push_back(T(1, eixdot, m_ctrlParams.at(velX).ki));
    //     dgdz.push_back(T(1, edxdot, m_ctrlParams.at(velX).kd));
    // }
    Eigen::SparseMatrix<double> dgdz_mat(NUM_Y_STATES, NUM_Z_STATES);
    dgdz_mat.setFromTriplets(dgdz.begin(), dgdz.end());
    return dgdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdxPlus(SystemState state, double timestep)
{
    std::vector<T> dhdx;
    dhdx.reserve(30);

    // TODO: what if ref is not constant?
    dhdx.push_back(T(edz, z, 1/timestep*state.plant(z)));

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

    Eigen::SparseMatrix<double> dhdx_mat(NUM_Z_STATES, NUM_PLANT_STATES);
    dhdx_mat.setFromTriplets(dhdx.begin(), dhdx.end());
    return dhdx_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzPlus(SystemState state, double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(89);
    dhdz.push_back(T(eix, eix, timestep*state.plant(eix)));

    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdzCurr(double timestep)
{
    std::vector<T> dhdz;
    dhdz.reserve(23);
    dhdz.push_back(T(eix, eix, timestep));

    Eigen::SparseMatrix<double> dhdz_mat(NUM_Z_STATES, NUM_Z_STATES);
    dhdz_mat.setFromTriplets(dhdz.begin(), dhdz.end());
    return dhdz_mat;
}

Eigen::SparseMatrix<double> DroneTrajectory::dhdy()
{
    Eigen::SparseMatrix<double> dhdy(NUM_Z_STATES, NUM_Y_STATES);
    return dhdy; 
}

