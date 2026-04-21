#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>

double windDist(double time)
{
    if(time < 1){
        return 1;
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
    return 1;
}

double negativeOneRef(double time){
    (void)time;
    return -1;
}

double oneHundredthRef(double time){
    (void)time;
    return 0.100;
}

double smoothStep(double time)
{
    double scale = 0.1;
    if (time <= 0 ){
        return 0;
    } else if (time < 1){
        return scale*(3*pow(time, 2) - 2*pow(time, 3));
    } else {
        return scale;
    }
}

SystemState initializeState()
{
    SystemState initialState;
    initialState.plant = Eigen::Vector<double, NUM_PLANT_STATES>::Zero();
    initialState.alge = Eigen::Vector<double, NUM_ALGE_STATES>::Zero();
    return initialState;
}

int main()
{
    Logger log("log.txt");
    std::array<double(*)(double), NUM_DIST_STATES> dist = {windDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    DroneTrajectory droneTrajectory(log, dist, ref);
    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    log << "INFO - simResults size: " << simResults.stateProgression.size() << std::endl;
    log << "INFO - time size: " << simResults.time.size() << std::endl;
    
    // Logger splot("splot.txt");
    // splotTrajectory(simResults, splot);

    std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    log << "INFO - trajSens size: " << ts.size() << std::endl;
    
    std::cout << ":D" << std::endl;
    return 0;
}