#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>
#include <cmath>

#include "load_experiment_traj.cpp"
#include "test_trajectory_from_file.cpp"


Eigen::Vector<double, NUM_DRONE_PARAMS> DroneTrajectory::getDroneParams()
{
    return { m_droneParams.Ix,m_droneParams.Iy, m_droneParams.Iz, m_droneParams.kf, m_droneParams.km };
}
void DroneTrajectory::setDroneParams(Eigen::Vector<double, NUM_DRONE_PARAMS> droneParams)
{
    m_droneParams.Ix = droneParams(Ix);
    m_droneParams.Iy = droneParams(Iy);
    m_droneParams.Iz = droneParams(Iz);
    m_droneParams.kf = droneParams(kf);
    m_droneParams.km = droneParams(km);
}

double noDist(double time){
    (void)time;
    return 0;
}

SystemState initializeState()
{
    SystemState initialState;
    initialState.plant = Eigen::Vector<double, NUM_PLANT_STATES>::Zero();
    initialState.alge = Eigen::Vector<double, NUM_ALGE_STATES>::Zero();
    return initialState;
}

double errorSquared(SimResults const & simResults)
{
    double error = 0;
    for(int i = 0; i < simResults.time.size(); i++)
    {
        error += std::pow((simResults.stateProgression.at(i).plant - interpolateState(simResults.time.at(i))).norm(), 2);
    }
    return error;
}

Eigen::Vector<double, NUM_DRONE_PARAMS> Jacobian(std::vector<dwddronep> const & trajSens, SimResults const & simResults)
{
    Eigen::Vector<double, NUM_DRONE_PARAMS> j;
    j.setZero();
    for(int i = 0; i < simResults.time.size(); i++)
    {
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_DRONE_PARAMS> xtrajSens = trajSens.at(i).dxdp;
        Eigen::Vector<double, NUM_PLANT_STATES> simPlantState = simResults.stateProgression.at(i).plant;
        Eigen::Vector<double, NUM_PLANT_STATES> actualPlantState = interpolateState(simResults.time.at(i));
        Eigen::Vector<double, NUM_PLANT_STATES> error = simPlantState - actualPlantState;
        j += xtrajSens.transpose()*error;
    }
    return j;
}

void systemId(Logger & log)
{
    std::string path1 = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/square_test/pos_pos.csv";
    std::string path2 = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/square_test/pos_att.csv";
    setActiveTrajectory(path1, path2);
    
    std::string path = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/square_test/traj.csv";
    std::vector<DataPoint> traj = loadTrajectory(path);

    setActiveTrajectory(traj);

    double finalTime = 10;
    double simTime = 1e-3;

    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {interpolateX, interpolateY, interpolateZ, interpolateYaw};
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);

    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    double error = errorSquared(simResults);
    log << "Error: " << error << std::endl;

    double stepsize = 1e-10;
    while(error > 10){
        std::vector<dwddronep> ts = droneTrajectory.trajSensDroneParam(simResults, simResults.time.size());
        Eigen::Vector<double, NUM_DRONE_PARAMS> j = Jacobian(ts, simResults);
        Eigen::Vector<double, NUM_DRONE_PARAMS> currParams = droneTrajectory.getDroneParams();
        Eigen::Vector<double, NUM_DRONE_PARAMS> newParams = currParams - stepsize*j;
        droneTrajectory.setDroneParams(newParams);

        simResults = droneTrajectory.Trajectory(initializeState());
        error = errorSquared(simResults);
        log << "currParams: " << currParams << std::endl;
        log << "newParams: " << newParams << std::endl;
        log << "j: " << j << std::endl;
        log << "Error: " << error << std::endl;
    }
}

void testTrajectory(Logger & log)
{
    std::string path = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/square_test/traj.csv";
    std::vector<DataPoint> traj = loadTrajectory(path);

    setActiveTrajectory(traj);

    double finalTime = 10;
    double simTime = 1e-3;

    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {interpolateX, interpolateY, interpolateZ, interpolateYaw};
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);

    SimResults simResults = droneTrajectory.Trajectory(initializeState(), false);
    log << "Timestamp,x,y,z" << std::endl;
    for(int i = 0; i < simResults.time.size(); i++)
    {
        Eigen::Vector<double, NUM_PLANT_STATES> plantState = simResults.stateProgression.at(i).plant;
        log << simResults.time.at(i) << "," << plantState(x)*1000 << "," << plantState(y)*1000 << "," << plantState(z)*1000 << std::endl;
    }
}

int main()
{
    Logger log("./build/log.txt");
    // systemId(log);
    testTrajectory(log);
    std::cout << ":D" << std::endl;
    return 0;
}