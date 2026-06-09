#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>

double windDist(double time)
{
    if(time < 1){
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

    initialState.plant(0) = 6.16166;
    initialState.plant(1) = 0.735878;
    initialState.plant(2) = -0.0207373;
    initialState.plant(3) = 0.13737;
    initialState.plant(4) = -0.369215;
    initialState.plant(5) = -0.00039356;
    initialState.plant(6) = 14.6402;
    initialState.plant(7) = 1.75168;
    initialState.plant(8) = 0.00533724;
    initialState.plant(9) = 1.53602;
    initialState.plant(10) = -0.0089044;
    initialState.plant(11) = -0.0421177;
    initialState.alge(0) = -0.881205;
    initialState.alge(1) = -14.5971;
    initialState.alge(2) = -10.3235;
    initialState.alge(3) = 0.638994;
    initialState.alge(4) = -2.04245;
    initialState.alge(5) = 0.524182;
    initialState.alge(6) = 0.00734216;
    initialState.alge(7) = -0.00531429;
    initialState.alge(8) = 0.0451457;
    initialState.alge(9) = -7.9455;
    initialState.alge(10) = -16.6931;
    initialState.alge(11) = 0.690496;
    initialState.alge(12) = 0.768229;
    initialState.alge(13) = 0.0361745;
    initialState.alge(14) = -0.0459063;
    initialState.alge(15) = 37537.8;
    initialState.alge(16) = 1.95924;
    initialState.alge(17) = -89.2444;
    initialState.alge(18) = 78.6534;
    initialState.alge(19) = -2.59825;
    initialState.alge(20) = 0.17759;
    initialState.alge(21) = -0.868145;
    initialState.alge(22) = -1.11389;
    initialState.alge(23) = 2.64994;
    initialState.alge(24) = -0.0511179;
    initialState.alge(25) = 5115.26;
    initialState.alge(26) = 5143.88;
    initialState.alge(27) = -41.7956;
    initialState.alge(28) = -42.2734;
    initialState.alge(29) = 0.783437;
    initialState.alge(30) = 573.228;
    initialState.alge(31) = -513.72;
    initialState.alge(32) = 0.208381;
    initialState.alge(33) = -4.71121;
    initialState.alge(34) = 2.92251;
    initialState.alge(35) = -0.319092;
    initialState.alge(36) = -22.0002;
    initialState.alge(37) = 278.117;
    initialState.alge(38) = 1583.44;
    initialState.alge(39) = 1560.19;
    initialState.alge(40) = 1561.95;
    initialState.alge(41) = 1538.94;
    initialState.alge(42) = 0.362681;
    initialState.alge(43) = -0.000161436;
    initialState.alge(44) = 2.11438e-06;
    initialState.alge(45) = 1.11666e-05;
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
    Logger log("./build/log.txt");
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    double finalTime = 1;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    
    std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    std::vector<dwdwo> tstest = droneTrajectory.trajSensTest(initializeState());
    double max = 0;
    int index = 0;
    for (int i = 0; i < simResults.time.size(); i++){
        double temp_max = (ts[i].dxdwo - tstest[i].dxdwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i; }
        temp_max = (ts[i].dydwo - tstest[i].dydwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i;  }
        temp_max = (ts[i].dzdwo - tstest[i].dzdwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i; }
    }
    log << "ts diff " << max << " " << index << std::endl;

    std::vector<dwdp> tsp = droneTrajectory.trajSensParam(simResults, simResults.time.size());
    std::vector<dwdp> tsptest = droneTrajectory.trajSensParamTest(initializeState());
    max = 0; 
    index = 0;
    for (int i = 0; i < simResults.time.size(); i++){
        double temp_max = (tsp[i].dxdp - tsptest[i].dxdp).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i; }
        temp_max = (tsp[i].dydp - tsptest[i].dydp).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i;  }
        temp_max = (tsp[i].dzdp - tsptest[i].dzdp).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i; }
    }
    log << "tsp diff " << max << " " << index << std::endl;

    std::vector<d2wdwo2> ts2 = droneTrajectory.secondOrdertrajSens(simResults, ts);
    std::vector<d2wdwo2> ts2test = droneTrajectory.secondOrdertrajSensTest(initializeState());
    max = 0; 
    index = 0;
    for (int i = 0; i < simResults.time.size(); i++){
        Eigen::Tensor<double, 0> max_tensor = (ts2[i].d2xdwo2 - ts2test[i].d2xdwo2).abs().maximum();
        double temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i; }
        max_tensor = (ts2[i].d2ydwo2 - ts2test[i].d2ydwo2).abs().maximum();
        temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i;  }
        max_tensor = (ts2[i].d2zdwo2 - ts2test[i].d2zdwo2).abs().maximum();
        temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i; }
    }
    log << "ts2 diff " << max << " " << index << std::endl;

    std::vector<d2wdwodp> ts2p = droneTrajectory.secondOrdertrajSensParams(simResults, ts, tsp);
    std::vector<d2wdwodp> ts2testp = droneTrajectory.secondOrdertrajSensParamsTest(initializeState());
    max = 0; 
    index = 0;   
    for (int i = 0; i < simResults.time.size(); i++){
        Eigen::Tensor<double, 0> max_tensor = (ts2p[i].d2xdwodp - ts2testp[i].d2xdwodp).abs().maximum();
        double temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i; }
        max_tensor = (ts2p[i].d2ydwodp - ts2testp[i].d2ydwodp).abs().maximum();
        temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i;  }
        max_tensor = (ts2p[i].d2zdwodp - ts2testp[i].d2zdwodp).abs().maximum();
        temp_max = max_tensor();
        if (temp_max > max) { max = temp_max; index = i; }
    }
    log << "ts2p diff " << max << " " << index << std::endl;
    G_tp gtp = droneTrajectory.calc_G_tp(ts);
    for(int i = 1; i < 500; i++)
    {
        gtp.tp = i;
        d2w testd2w = droneTrajectory.calc_d2w(simResults, ts, gtp, tsp, ts2, ts2p);
        Eigen::Tensor<double, 0> diff_tensor = (testd2w.dwo2.d2xdwo2 - ts2.at(i).d2xdwo2).abs().maximum();
        double diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2xdwo2 " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwo2.d2ydwo2 - ts2.at(i).d2ydwo2).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2ydwo2 " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwo2.d2zdwo2 - ts2.at(i).d2zdwo2).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2zdwo2 " << diff_double << " " << i << std::endl; }

        diff_tensor = (testd2w.dwodp.d2xdwodp - ts2p.at(i).d2xdwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2xdwodp " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwodp.d2ydwodp - ts2p.at(i).d2ydwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2xdwodp " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwodp.d2zdwodp - ts2p.at(i).d2zdwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > 0 ) { log << "d2zdwodp " << diff_double << " " << i << std::endl; }
    }
    log << "finished comparing d2w" << std::endl;
    gtp = droneTrajectory.calc_G_tp(ts);
    d2w d2w = droneTrajectory.calc_d2w(simResults, ts, gtp, tsp, ts2, ts2p);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = droneTrajectory.calc_dG(ts.at(gtp.tp), d2w,  gtp);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG_test = droneTrajectory.calc_dG_test(initializeState(), ts.at(gtp.tp), gtp, finalTime);
    
    log << "dG diff " << (dG-dG_test).cwiseAbs().maxCoeff() << std::endl;
       
    std::cout << ":D" << std::endl;
    return 0;
}