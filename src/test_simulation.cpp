#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>

double windDist(double time)
{
    if(time < 0.001){
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

SystemState stateCloseToRoABoundary()
{
    SystemState initialState;

    initialState.plant(0) = 9.42579;
    initialState.plant(1) = 1.04072;
    initialState.plant(2) = -0.0206441;
    initialState.plant(3) = 0.311541;
    initialState.plant(4) = -0.369181;
    initialState.plant(5) = -0.00404734;
    initialState.plant(6) = 17.9918;
    initialState.plant(7) = 1.23546;
    initialState.plant(8) = -0.00804218;
    initialState.plant(9) = 0.423235;
    initialState.plant(10) = 0.00129618;
    initialState.plant(11) = -0.00170205;
    initialState.alge(0) = -2.23041;
    initialState.alge(1) = -17.9819;
    initialState.alge(2) = -16.8511;
    initialState.alge(3) = 0.655264;
    initialState.alge(4) = -1.32267;
    initialState.alge(5) = -0.149645;
    initialState.alge(6) = 0.0114291;
    initialState.alge(7) = 0.00800204;
    initialState.alge(8) = 0.0470027;
    initialState.alge(9) = -13.9088;
    initialState.alge(10) = -6.52761;
    initialState.alge(11) = 0.408604;
    initialState.alge(12) = 3.29112;
    initialState.alge(13) = 0.0452016;
    initialState.alge(14) = 0.0802768;
    initialState.alge(15) = 38054.1;
    initialState.alge(16) = 3.17203;
    initialState.alge(17) = -24.3563;
    initialState.alge(18) = 22.4164;
    initialState.alge(19) = -2.36685;
    initialState.alge(20) = -0.100059;
    initialState.alge(21) = -0.185315;
    initialState.alge(22) = -1.08001;
    initialState.alge(23) = 0.0776442;
    initialState.alge(24) = 0.338539;
    initialState.alge(25) = 1437.5;
    initialState.alge(26) = 1446.75;
    initialState.alge(27) = -21.6226;
    initialState.alge(28) = -21.666;
    initialState.alge(29) = -0.173573;
    initialState.alge(30) = 161.226;
    initialState.alge(31) = -142.028;
    initialState.alge(32) = 0.145901;
    initialState.alge(33) = -2.41442;
    initialState.alge(34) = 2.01932;
    initialState.alge(35) = -0.0718121;
    initialState.alge(36) = -4.1235;
    initialState.alge(37) = 51.1278;
    initialState.alge(38) = 1587.72;
    initialState.alge(39) = 1583.39;
    initialState.alge(40) = 1581.73;
    initialState.alge(41) = 1577.57;
    initialState.alge(42) = 0.372691;
    initialState.alge(43) = -4.52449e-05;
    initialState.alge(44) = 7.04058e-07;
    initialState.alge(45) = 2.08106e-06;
    initialState.alge(46) = 20;
    initialState.alge(47) = -20;

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
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    double finalTime = 100;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(stateCloseToRoABoundary());
    // std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    // G_tp gtp = droneTrajectory.calc_G_tp(ts);
    // Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = droneTrajectory.calc_dG_test(initializeState(), ts.at(gtp.tp), gtp, gtp.tp*1e-3);

    // SystemState initialState = initializeState();
    // Eigen::Vector<double, NUM_STATES> z0;
    // z0 << initialState.plant, initialState.alge;
    // Eigen::Vector<double, NUM_STATES> old_zk = droneTrajectory.closestZBar(initializeState());
    // log << "initial zbar dist " << (z0 - old_zk).norm() << std::endl; 
    zkpk optimal = droneTrajectory.theGigaAlgo(stateCloseToRoABoundary());
    // Eigen::Vector<double, NUM_STATES> new_zk = droneTrajectory.closestZBar(initializeState());
    // log << "new zbar dist " << (z0 - new_zk).norm() << std::endl; 
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    log << "Elapsed Time ERA algo: " << elapsed.count() << " us" << std::endl;
    // for(int i = 0; i < NUM_PARAMETERS; i++) { log << i << ": " << optimal.pk(i) << std::endl; }

    simResults = droneTrajectory.Trajectory(initializeState());
    log << "stable? " << simResults.stable << std::endl;
    log << "converged? " << simResults.converged << std::endl;
    log << "INFO - simResults size: " << simResults.stateProgression.size() << std::endl;
    log << "INFO - time size: " << simResults.time.size() << std::endl;

    // double timeIndex = 1001;
    // log <<       "x: " << simResults.stateProgression.at(timeIndex).plant(x) <<
    //             " y: " << simResults.stateProgression.at(timeIndex).plant(y) <<
    //             " z: " << simResults.stateProgression.at(timeIndex).plant(z) <<
    //             " phi: " << simResults.stateProgression.at(timeIndex).plant(phi) <<
    //             " theta: " << simResults.stateProgression.at(timeIndex).plant(theta) <<
    //             " psi: " << simResults.stateProgression.at(timeIndex).plant(psi) <<
    //             " xdot: " << simResults.stateProgression.at(timeIndex).plant(xdot) <<
    //             " ydot: " << simResults.stateProgression.at(timeIndex).plant(ydot) <<
    //             " zdot: " << simResults.stateProgression.at(timeIndex).plant(zdot) <<
    //             " p: " << simResults.stateProgression.at(timeIndex).plant(p) <<
    //             " q: " << simResults.stateProgression.at(timeIndex).plant(q) <<
    //             " r: " << simResults.stateProgression.at(timeIndex).plant(r) << std::endl;

    // for(int i = 0; i < NUM_PLANT_STATES; i++)
    // { log << "initialState.plant(" << i << ") = " << simResults.stateProgression.at(timeIndex).plant(i) << ";" << std::endl; }

    // for(int i = 0; i < NUM_ALGE_STATES; i++)
    // { log << "initialState.alge(" << i << ") = " << simResults.stateProgression.at(timeIndex).alge(i) << ";" << std::endl; }
    
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