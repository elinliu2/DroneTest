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

void printDwdwo(Logger & log, dwdwo trajSens)
{
    // log << trajSens.dxdwo << std::endl;
    log << trajSens.dzdwo << std::endl;
    // log << trajSens.dydwo << std::endl;
}

void diffDwdwo(Logger & log, dwdwo trajSensA, dwdwo trajSensB)
{
    log << "dxdwo diff" << std::endl;
    log << trajSensA.dxdwo - trajSensB.dxdwo << std::endl;
    log << "dzdwo diff" << std::endl;
    log << trajSensA.dzdwo - trajSensB.dzdwo << std::endl;
    // log << trajSensA.dydwo - trajSensB.dydwo << std::endl;
}

std::pair<double, int> diffTrajSens(std::vector<dwdwo> const & tsA, std::vector<dwdwo> const & tsB)
{
    int numIterations = std::min(tsA.size(), tsB.size());
    double maxDiff = 0;
    int index = 0;
    for (int i = 0; i < numIterations; i++)
    {
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwoDiff = tsA[i].dxdwo-tsB[i].dxdwo;
        Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwoDiff = tsA[i].dzdwo-tsB[i].dzdwo;
        Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwoDiff = tsA[i].dydwo-tsB[i].dydwo;
        
        if (dxdwoDiff.maxCoeff() > maxDiff){
            maxDiff = dxdwoDiff.maxCoeff();
            index = i;
        }
        if (dydwoDiff.maxCoeff() > maxDiff){
            maxDiff = dydwoDiff.maxCoeff();
            index = i;
        }
        if (dzdwoDiff.maxCoeff() > maxDiff){
            maxDiff = dzdwoDiff.maxCoeff();
            index = i;
        }

        if ((-1*dxdwoDiff.maxCoeff()) > maxDiff){
            maxDiff = (-1*dxdwoDiff.maxCoeff());
            index = i;
        }
        if ((-1*dydwoDiff.maxCoeff()) > maxDiff){
            maxDiff = (-1*dydwoDiff.maxCoeff());
            index = i;
        }
        if ((-1*dzdwoDiff.maxCoeff()) > maxDiff){
            maxDiff = (-1*dzdwoDiff.maxCoeff());
            index = i;
        }

    }
    return {maxDiff, index};
}


std::pair<double, int> diffdxdwo(std::vector<dwdwo> const & tsA, std::vector<dwdwo> const & tsB)
{
    int numIterations = std::min(tsA.size(), tsB.size());
    double maxDiff = 0;
    int index = 0;
    for (int i = 0; i < numIterations; i++)
    {
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwoDiff = tsA[i].dxdwo-tsB[i].dxdwo;
     
        if (dxdwoDiff.maxCoeff() > maxDiff){
            maxDiff = dxdwoDiff.maxCoeff();
            index = i;
        }

        if ((-1*dxdwoDiff.maxCoeff()) > maxDiff){
            maxDiff = (-1*dxdwoDiff.maxCoeff());
            index = i;
        }
    }
    return {maxDiff, index};
}

std::pair<double, int> diffdzdwo(std::vector<dwdwo> const & tsA, std::vector<dwdwo> const & tsB)
{
    int numIterations = std::min(tsA.size(), tsB.size());
    double maxDiff = 0;
    int index = 0;
    for (int i = 0; i < numIterations; i++)
    {
        Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwoDiff = tsA[i].dzdwo-tsB[i].dzdwo;
        
        if (dzdwoDiff.maxCoeff() > maxDiff){
            maxDiff = dzdwoDiff.maxCoeff();
            index = i;
        }

        if ((-1*dzdwoDiff.maxCoeff()) > maxDiff){
            maxDiff = (-1*dzdwoDiff.maxCoeff());
            index = i;
        }

    }
    return {maxDiff, index};
}
int main()
{
    Logger log("log.txt");
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    DroneTrajectory droneTrajectory(log, dist, ref);
    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    log << "INFO - simResults size: " << simResults.stateProgression.size() << std::endl;
    log << "INFO - time size: " << simResults.time.size() << std::endl;
    
    // Logger splot("splot.txt");
    // splotTrajectory(simResults, splot);

    std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    log << "INFO - trajSens size: " << ts.size() << std::endl;

    std::vector<dwdwo> tsTest = droneTrajectory.trajSensTest(initializeState());
    std::pair<double, int> diff = diffTrajSens({ts.at(0)}, {tsTest.at(0)});
    log << "max diff ts 0: " << diff.first << std::endl;
    diff = diffTrajSens({ts.at(1)}, {tsTest.at(1)});
    log << "max diff ts 1: " << diff.first <<  std::endl;
    diff = diffTrajSens({ts.at(2)}, {tsTest.at(2)});
    log << "max diff ts 2: " << diff.first <<  std::endl;

    
    // diff = diffTrajSens({ts.at(411)}, {tsTest.at(411)});
    // log << "max diff ts 411: " << diff.first <<  std::endl;
    // diff = diffTrajSens({ts.at(412)}, {tsTest.at(412)});
    // log << "max diff ts 412: " << diff.first <<  std::endl;
    diff = diffTrajSens(ts, tsTest);
    log << "max diff: " << diff.first <<  std::endl;
    
    // diff = diffTrajSens({ts.at(2)}, {tsTest.at(2)});
    // log << "max diff ts 2: " << diff.first <<  std::endl;
    // diffDwdwo(log, ts.at(1), tsTest.at(1));

    // diff = diffTrajSens(ts, tsTest);
    // log << "max diff ts: " << diff.first <<  " index: " << diff.second << std::endl;    
   
    // diff = diffdxdwo({ts.at(2)}, {tsTest.at(2)});
    // log << "max diff dxdwo: " << diff.first <<  " index: " << diff.second << std::endl;  
    // for (int i = 0; i < NUM_PLANT_STATES; i++){
    //     log << "row " << i << std::endl;
    //     log << ts.at(410).dxdwo.row(i) << std::endl;
    //     log << tsTest.at(410).dxdwo.row(i) << std::endl;
    //     log << (ts.at(410).dxdwo.row(i)-tsTest.at(410).dxdwo.row(i)).maxCoeff() << std::endl;
    //     log << (-ts.at(410).dxdwo.row(i)+tsTest.at(410).dxdwo.row(i)).maxCoeff() << std::endl;  
    // } 

    // for (int i = 0; i < NUM_Y_STATES; i++){
    //     log << "row " << i << std::endl;
    //     log << ts.at(410).dydwo.row(i) << std::endl;
    //     log << tsTest.at(410).dydwo.row(i) << std::endl;
    //     log << (ts.at(410).dydwo.row(i)-tsTest.at(410).dydwo.row(i)).maxCoeff() << std::endl;
    //     log << (-ts.at(410).dydwo.row(i)+tsTest.at(410).dydwo.row(i)).maxCoeff() << std::endl;  
    // } 

    // diff = diffdzdwo({ts.at(2)}, {tsTest.at(2)});
    // log << "max diff dzdwo: " << diff.first <<  " index: " << diff.second << std::endl;    

    // for (int i = 0; i < NUM_Z_STATES; i++){
    //     log << "row " << i << std::endl;
    //     log << ts.at(410).dzdwo.row(i) << std::endl;
    //     log << tsTest.at(410).dzdwo.row(i) << std::endl;
    //     log << ( ts.at(410).dzdwo.row(i)-tsTest.at(410).dzdwo.row(i)).maxCoeff() << std::endl;
    //     log << (-ts.at(410).dzdwo.row(i)+tsTest.at(410).dzdwo.row(i)).maxCoeff() << std::endl;
    // }

    // for (int i = 0; i < 10001; i++){
    //     diff = diffTrajSens({ts.at(i)}, {tsTest.at(i)});
    //     log << diff.first << std::endl;
    // }
    
    // log << "INFO - diff" << std::endl;
    // diffDwdwo(log, ts.at(10001), tsTest.at(10001));
    
    std::cout << ":D" << std::endl;
    return 0;
}