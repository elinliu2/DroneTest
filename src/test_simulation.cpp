#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>

double xDist(double time)
{
    if(time < 0.15){
        return 10;
    }
    return 0;
}

double yDist(double time)
{
    if(time < 0.15){
        return 10;
    }
    return 0;
}

double noDist(double time){
    (void)time;
    return 0;
}

double zeroRef(double time){
    (void)time;
    return 0;
}

double oneRef(double time){
    (void)time;
    return 0.01;
}

double oneHundredRef(double time){
    (void)time;
    return 0.101;
}

SystemState initializeState()
{
    SystemState initialState;
    initialState.plant = Eigen::Vector<double, NUM_PLANT_STATES>::Zero();
    initialState.alge = Eigen::Vector<double, NUM_ALGE_STATES>::Zero();
    for (int i = 0; i < NUM_CTRL_STATES; i++){
        PIDstate empty;
        initialState.ctrl(i) = empty;
    }
    initialState.plant(z) = 0.100;
    return initialState;
}

int main()
{
    Logger log("log.txt");
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {zeroRef, zeroRef, oneHundredRef};
    DroneTrajectory droneTrajectory(log, dist, ref);
    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    log << "INFO - simResults size: " << simResults.stateProgression.size() << std::endl;
    
    Logger splot("splot.txt");
    splotTrajectory(simResults, splot);

    Logger splotz("z.txt");
    splotPlantState(simResults, splotz, z);

    std::cout << ":D" << std::endl;
    return 0;
}