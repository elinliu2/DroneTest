#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>

double windDist(double time)
{
    if(time < 1.003){
        return 0.7;
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

std::array<PIDParameters, NUM_PIDS> initialParameters(){
    PIDParameters posXpid = {1, 0, 0};
    PIDParameters posYpid = {1, 0, 0};
    PIDParameters posZpid = {1, 0, 0};

    PIDParameters velXpid = {1, 0, 0};
    PIDParameters velYpid = {1, 0, 0};
    PIDParameters velZpid = {1, 0, 0};

    PIDParameters rollpid  = {10, 0, 0, 20.0};
    PIDParameters pitchpid = {10, 0, 0, 20.0};
    PIDParameters yawpid   = {10, 0, 0, 360.0};

    PIDParameters rollRatepid  = {10, 0, 0, 33.3, true};
    PIDParameters pitchRatepid = {10, 0, 0, 33.3, true};
    PIDParameters yawRatepid   = {10, 0, 0, 166.7};
    return {posXpid, posYpid, posZpid, velXpid, velYpid, velZpid, rollpid, pitchpid, yawpid, rollRatepid, pitchRatepid, yawRatepid};
}

int main()
{
    // std::cout << "cwd: " << std::filesystem::current_path() << std::endl;
    Logger log("./build/log.txt");
    std::array<double(*)(double), NUM_DIST_STATES> dist = {windDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    double finalTime = 500;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(initializeState(), true);
    // std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    // G_tp gtp = droneTrajectory.calc_G_tp(ts);
    // Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = droneTrajectory.calc_dG_test(initializeState(), ts.at(gtp.tp), gtp, gtp.tp*1e-3);

    // SystemState initialState = initializeState();
    // Eigen::Vector<double, NUM_STATES> z0;
    // z0 << initialState.plant, initialState.alge;
    // Eigen::Vector<double, NUM_STATES> old_zk = droneTrajectory.closestZBar(initializeState());
    // log << "initial zbar dist " << (z0 - old_zk).norm() << std::endl; 
    // zkpk optimal = droneTrajectory.theGigaAlgo(initializeState());
    // Eigen::Vector<double, NUM_STATES> new_zk = droneTrajectory.closestZBar(initializeState());
    // log << "new zbar dist " << (z0 - new_zk).norm() << std::endl; 
    // std::chrono::time_point end = std::chrono::steady_clock::now();
    // std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // log << "Elapsed Time ERA algo: " << elapsed.count() << " us" << std::endl;
    // for(int i = 0; i < NUM_PARAMETERS; i++) { log << i << ": " << optimal.pk(i) << std::endl; }

    // simResults = droneTrajectory.Trajectory(initializeState());
    log << "stable? " << simResults.stable << std::endl;
    log << "converged? " << simResults.converged << std::endl;
    log << "INFO - simResults size: " << simResults.stateProgression.size() << std::endl;
    log << "INFO - time size: " << simResults.time.size() << std::endl;

    log <<       "x: " << simResults.stateProgression.at(simResults.time.size()-1).plant(x) <<
                " y: " << simResults.stateProgression.at(simResults.time.size()-1).plant(y) <<
                " z: " << simResults.stateProgression.at(simResults.time.size()-1).plant(z) <<
                " phi: " << simResults.stateProgression.at(simResults.time.size()-1).plant(phi) <<
                " theta: " << simResults.stateProgression.at(simResults.time.size()-1).plant(theta) <<
                " psi: " << simResults.stateProgression.at(simResults.time.size()-1).plant(psi) <<
                " xdot: " << simResults.stateProgression.at(simResults.time.size()-1).plant(xdot) <<
                " ydot: " << simResults.stateProgression.at(simResults.time.size()-1).plant(ydot) <<
                " zdot: " << simResults.stateProgression.at(simResults.time.size()-1).plant(zdot) <<
                " p: " << simResults.stateProgression.at(simResults.time.size()-1).plant(p) <<
                " q: " << simResults.stateProgression.at(simResults.time.size()-1).plant(q) <<
                " r: " << simResults.stateProgression.at(simResults.time.size()-1).plant(r) << std::endl;
    
    Logger splot("./build/splot.txt");
    splotTrajectory(simResults, splot);
    // std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);

    Logger xPlot("./build/x.txt");
    splotPlantState(simResults, xPlot, x);

    Logger yPlot("./build/y.txt");
    splotPlantState(simResults, yPlot, y);
    
    // log << "INFO - trajSens size: " << ts.size() << std::endl;
    
    // G_tp gtp = droneTrajectory.calc_G_tp(ts);

    // std::vector<dwdwo> tsTest = droneTrajectory.trajSensTest(initializeState());
    // std::pair<double, int> diff = diffTrajSens(ts, tsTest);
    // log << "max diff: " << diff.first <<  std::endl;

    // std::vector<dwdp> ts_p = droneTrajectory.trajSensParam(simResults, gtp);
    // std::vector<d2wdwo2> ts2Test = droneTrajectory.secondOrdertrajSensTest(initializeState());
    // std::vector<d2wdwodp> ts2ParamsTest = droneTrajectory.secondOrdertrajSensParamsTest(initializeState());
    // Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = droneTrajectory.calc_dG_test(initializeState(), ts.at(gtp.tp), ts_p.at(gtp.tp), gtp, simResults.time.at(gtp.tp));
    
    std::cout << ":D" << std::endl;
    return 0;
}