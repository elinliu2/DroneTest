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
    (void) state;
    return dgdz;
}

