#include "DroneTrajectory.h"
#include "Splotting.h"
#include "Logger.h"
#include <iostream>
#include <cmath>
#include "test_trajectory_from_file.cpp"

double windDist(double time)
{
    if(time < 1){
        return 0.525;
    }

    // if(time < 0.1){
    //     return 0.019;
    // }

    return 0;
}

double testWindDist(double time)
{
    if(time < 0.05){
        return 0.2;
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

double ninetyRef(double time){
    (void)time;
    return 90;
}

double negativeOneRef(double time){
    (void)time;
    return -1;
}

double oneHundredthRef(double time){
    (void)time;
    return 0.01;
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

SystemState stateOnRoABoundary()
{
    SystemState initialState;

    initialState.plant(0) = 9.4078;
    initialState.plant(1) = 1.03948;
    initialState.plant(2) = -0.0206361;
    initialState.plant(3) = 0.311115;
    initialState.plant(4) = -0.369183;
    initialState.plant(5) = -0.00404598;
    initialState.plant(6) = 17.9853;
    initialState.plant(7) = 1.2388;
    initialState.plant(8) = -0.0079619;
    initialState.plant(9) = 0.425983;
    initialState.plant(10) = 0.00125436;
    initialState.plant(11) = -0.00177402;
    initialState.alge(eix) = -2.22199;
    initialState.alge(edx) = -17.9702;
    initialState.alge(desVelX) = -16.8151;
    initialState.alge(eiy) = 0.655338;
    initialState.alge(edy) = -1.32674;
    initialState.alge(desVelY) = -0.147003;
    initialState.alge(eiz) = 0.0114084;
    initialState.alge(edz) = 0.00792153;
    initialState.alge(desVelZ) = 0.0469763;

    initialState.alge(eixdot) = -13.874;
    initialState.alge(edxdot) = -16.8217;
    initialState.alge(eiydot) = 0.410062;
    initialState.alge(edydot) = 3.24285;
    initialState.alge(eizdot) = 0.0451465;
    initialState.alge(edzdot) = 0.0807363;
    initialState.alge(desThrust) = 38050.7;

    initialState.alge(eiphi) = 3.16988;
    initialState.alge(edphi) = -24.5161;
    initialState.alge(desRollRate) = 22.556;
    initialState.alge(eitheta) = -2.36801;
    initialState.alge(edtheta) = -0.0990131;
    initialState.alge(desPitchRate) = -0.188172;
    initialState.alge(desYawRef) = 0.0;
    initialState.alge(eipsi) = -1.08024;
    initialState.alge(edpsi) = 0.0827114;
    initialState.alge(desYawRate) = 0.339614;

    initialState.alge(delay_1_rollRate) = 1446.75;
    initialState.alge(delay_2_rollRate) = 1456.06;
    initialState.alge(delay_1_pitchRate) = -21.666;
    initialState.alge(delay_2_pitchRate) = -21.7089;

    initialState.alge(eip) = -0.17174;
    initialState.alge(edp) = 162.264;
    initialState.alge(desRollOutput) = -142.955;
    initialState.alge(eiq) = 0.146161;
    initialState.alge(edq) = -2.41921;
    initialState.alge(desPitchOutput) = 2.0219;
    initialState.alge(eir) = -0.0722482;
    initialState.alge(edr) = -4.1726;
    initialState.alge(desYawOutput) = 51.7445;

    initialState.alge(w1) = 1587.62;
    initialState.alge(w2) = 1583.24;
    initialState.alge(w3) = 1581.6;
    initialState.alge(w4) = 1577.38;

    initialState.alge(ft) = 0.372623;
    initialState.alge(tx) = -4.5536e-05;
    initialState.alge(ty) = 7.05955e-07;
    initialState.alge(tz) = 2.10597e-06;

    initialState.alge(desRoll) = 20;
    initialState.alge(desPitch) = -20;
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

void testTrajSens(Logger & log)
{
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    double finalTime = 10;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(stateCloseToRoABoundary());

    Logger splot("./build/splot.txt");
    splotTrajectory(simResults, splot);
    
    std::vector<dwdwo> ts = droneTrajectory.trajSens(simResults);
    std::vector<dwdwo> tstest = droneTrajectory.trajSensTest(stateCloseToRoABoundary());
    double max = 0;
    int index = 0;
    for (int i = 0; i < ts.size(); i++){
        double temp_max = (ts[i].dxdwo - tstest[i].dxdwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max;  index = i;}
        temp_max = (ts[i].dydwo - tstest[i].dydwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max;  index = i;}
        temp_max = (ts[i].dzdwo - tstest[i].dzdwo).cwiseAbs().maxCoeff();
        if (temp_max > max) { max = temp_max; index = i;}
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
    G_tp gtp = droneTrajectory.calc_G_tp(tstest);
    log << "tp: " << gtp.tp << std::endl;
    double tol = 5e-8;
    for(int i = 1; i < 50; i++)
    {
        gtp.tp = i;
        d2w testd2w = droneTrajectory.calc_d2w(simResults, ts, gtp);
        Eigen::Tensor<double, 0> diff_tensor = (testd2w.dwo2.d2xdwo2 - ts2.at(i).d2xdwo2).abs().maximum();
        double diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2xdwo2 " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwo2.d2ydwo2 - ts2.at(i).d2ydwo2).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2ydwo2 " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwo2.d2zdwo2 - ts2.at(i).d2zdwo2).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2zdwo2 " << diff_double << " " << i << std::endl; }

        diff_tensor = (testd2w.dwodp.d2xdwodp - ts2p.at(i).d2xdwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2xdwodp " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwodp.d2ydwodp - ts2p.at(i).d2ydwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2xdwodp " << diff_double << " " << i << std::endl; }
        diff_tensor = (testd2w.dwodp.d2zdwodp - ts2p.at(i).d2zdwodp).abs().maximum();
        diff_double = diff_tensor();
        if ( diff_double > tol ) { log << "d2zdwodp " << diff_double << " " << i << std::endl; }
    }
    log << "finished comparing d2w" << std::endl;
    gtp = droneTrajectory.calc_G_tp(ts);
    d2w d2w = droneTrajectory.calc_d2w(simResults, ts, gtp);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = droneTrajectory.calc_dG(ts.at(gtp.tp), d2w,  gtp);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG_test = droneTrajectory.calc_dG_test(initializeState(), ts.at(gtp.tp), gtp, finalTime);
    
    log << "dG diff " << (dG-dG_test) << std::endl;

    for(int i = 0; i < NUM_STATES+NUM_PARAMETERS; i++)
    {
        if(dG(i)!= 0){
            log << "i " << i << "dG diff %" << (dG(i)-dG_test(i))/dG(i) << std::endl;
        }
    }
}

void testERAAlgo(Logger & log)
{
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::string path = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/tall_hover/traj.csv";
    std::vector<DataPoint> traj = loadTrajectory(path);
    setActiveTrajectory(traj);
    std::array<double(*)(double), NUM_REF_STATES> ref = {interpolateX, interpolateY, interpolateZ, interpolateYaw};

    double finalTime = 12;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(initializeState());
    log << simResults.stable << std::endl;
    zkpk zkpk = droneTrajectory.theGigaAlgo(initializeState());
    log << "zkpk" << std::endl;
    log << zkpk.zk << std::endl << std::endl;
    log << zkpk.pk << std::endl;
}

void testSim(Logger & log)
{
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {zeroRef, zeroRef, oneRef, zeroRef};
    double finalTime = 100;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(stateCloseToRoABoundary(), true);

    Logger splot("./build/splot.txt");
    splotTrajectory(simResults, splot);

    Logger xPlot("./build/x.txt");
    splotPlantState(simResults, xPlot, x);

    Logger yPlot("./build/y.txt");
    splotPlantState(simResults, yPlot, y);
}

Eigen::Vector<double, NUM_PARAMETERS> get_test_param()
{
    // initial state sf
    // Eigen::Vector<double, NUM_PARAMETERS> params = {
    //             0.999881,
    //                    1,
    //                    1,
    //             0.998236,
    //                    1,
    //                    1,
    //             0.969074,
    //             0.997806,
    //                    1,
    //               0.9999,
    //             0.999999,
    //                    1,
    //             0.998259,
    //             0.999989,
    //                    1,
    //             0.874869,
    //             0.994998,
    //                    1,
    //             0.957454,
    //             0.999577,
    //                    1,
    //             0.961113,
    //             0.999624,
    //                    1,
    //             0.948504,
    //             0.999708,
    //             0.884178,
    //             0.948737,
    //             0.998083,
    //             0.973266,
    //             0.952399,
    //             0.998214,
    //             0.973419,
    //             0.107803,
    //             0.996855,
    //                    1
    // };

    // Eigen::Vector<double, NUM_PARAMETERS> params = {
    //     0.985795,
    //         1,
    //         1,
    //     1.00001,
    //         1,
    //         1,
    //         1,
    //     0.999999,
    //         1,
    //     1.00307,
    //     0.992643,
    //         1,
    //     0.999927,
    //         1,
    //         1,
    //     1.00002,
    //         1,
    //         1,
    //         1,
    //     0.999996,
    //         1,
    //     0.999244,
    //     1.00171,
    //         1,
    //     1.00007,
    //     1.00001,
    //     0.999986,
    //         1,
    //         1,
    //     0.999999,
    //     0.99987,
    //     0.999995,
    //     1.00003,
    //     1.00003,
    //         1,
    //         1
    // };

    // square traj
    Eigen::Vector<double, NUM_PARAMETERS> params = {
        0.949178,
            1,
            1,
        0.955509,
            1,
            1,
        0.988354,
        0.997718,
            1,
        0.936882,
        0.998657,
            1,
        0.944227,
        0.998879,
            1,
        0.978727,
        0.996415,
            1,
        0.908343,
        0.999385,
            1,
        0.91211,
        1.00034,
            1,
        0.958529,
        0.999797,
        0.976927,
        0.765037,
        1.00315,
        0.745839,
        0.769431,
        1.00403,
        0.766602,
        0.882404,
        0.99906,
            1,
    };
    return params;
}

void test_param(Logger & log)
{
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {zeroRef, oneRef, oneRef, ninetyRef};
    double finalTime = 10;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    // droneTrajectory.m_sf = get_test_param();
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(initializeState());

    log << "stable? " << simResults.stable << std::endl;

    Logger splot("./build/splot.txt");
    splotTrajectory(simResults, splot, "splot");

    Logger xPlot("./build/x.txt");
    splotPlantState(simResults, xPlot, x, "x");

    Logger yPlot("./build/y.txt");
    splotPlantState(simResults, yPlot, y, "y");

    Logger zPlot("./build/z.txt");
    splotPlantState(simResults, zPlot, z, "z");

    Logger psiPlot("./build/psi.txt");
    splotPlantState(simResults, psiPlot, psi, "Psi");
}

double squareRefx(double time){
    if (time <= 4){
        return 0;
    } 
    if (time <= 6){
        return 1;
    } 
    return 0;
}

double squareRefy(double time){
    if (time <= 5){
        return 0;
    } 
    if (time <= 7){
        return 1;
    } 
    return 0;
}

double squareRefz(double time){
    if (time <= 8){
        return 0.4;
    } 
    return 0;
}

double squareRefyaw(double time){
    (void) time;
    if (time <= 3){
        return 0;
    } 
    if (time <= 4){
        return 90;
    } 
    // if (time <= 5){
    //     return 180;
    // } 
    // if (time <= 6){
    //     return 270;
    // } 
    return 0;
}

void test_sysid(Logger & log)
{
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {squareRefx, squareRefy, squareRefz, squareRefyaw};
    double finalTime = 10;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    SimResults simResults = droneTrajectory.Trajectory(initializeState(), false);
    log << "stable? " << simResults.stable << std::endl;
    log << "Timestamp,x,y,z" << std::endl;
    for(int i = 0; i < simResults.time.size(); i++){
        log << simResults.time.at(i) << ", " 
            << simResults.stateProgression.at(i).plant(x) << ", "
            << simResults.stateProgression.at(i).plant(y) << ", "
            << simResults.stateProgression.at(i).plant(z) << std::endl;
    }
    Logger psiPlot("./build/psi.txt");
    splotPlantState(simResults, psiPlot, psi, "psi");
}

double hover(double time)
{
    if (time < 4)    { return interpolateZ(time); }
    return 0.4;
}

double capAngleTest(double angle)
{
  double result = angle;

  while (result > 180.0) {
    result -= 360;
  }

  while (result < -180.0) {
    result += 360;
  }

  return result;
}

void test_square_traj(Logger & log)
{
    std::string path = "../crazyflie_rl_sim/workspace/crazyflie_rl_sim/recordings/system_id_logs/square_full_actual/traj.csv";
    std::vector<DataPoint> traj = loadTrajectory(path);

    // log << "Timestamp,x,y,z,yaw" << std::endl;
    // for(int i = 0; i < traj.size(); i++){
    //    DataPoint dp = traj.at(i); 
    //    log << dp.timestamp << "," << dp.x << "," << dp.y << ","  << dp.z << "," << dp.yaw << std::endl;
    // }

    setActiveTrajectory(traj);

    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {interpolateX, interpolateY, interpolateZ, interpolateYaw};
    double finalTime = 17.5;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    // droneTrajectory.m_sf = get_test_param();

    SystemState ic = initializeState();
    // ic.plant = {-6.0/1000, 1.0/1000, 1.0/1000, -0.46661794/180*M_PI, -0.21732523/180*M_PI, 2.3516355/180*M_PI, -4.0/1000,1.0/1000,2.0/1000,7.0/1000,0.0/1000,0.0/1000};
    ic.plant = {-3.0/1000,-2.0/1000,1.0/1000,-0.3299162/180*M_PI,0.1375021/180*M_PI,0.71228117/180*M_PI,1.0/1000,-4.0/1000,1.0/1000,4.0/1000,-2.0/1000,3.0/1000};
    SimResults simResults = droneTrajectory.Trajectory(ic, false);
    // log << "stable? " << simResults.stable << std::endl;

    Logger splot("./build/splot.txt");
    splotTrajectory(simResults, splot, "testSquareTraj"); 

    Logger xPlot("./build/x.txt");
    splotPlantState(simResults, xPlot, x, "x");

    Logger yPlot("./build/y.txt");
    splotPlantState(simResults, yPlot, y, "y");

    Logger zPlot("./build/z.txt");
    splotPlantState(simResults, zPlot, z, "z");

    log << "Timestamp,x,y,z,roll,pitch,yaw,xdot,ydot,zdot,p,q,r" << std::endl;
    for(int i = 0; i < simResults.time.size(); i++){
        log << simResults.time.at(i) << ", " 
            << simResults.stateProgression.at(i).plant(x)*1000 << ","
            << simResults.stateProgression.at(i).plant(y)*1000 << ","
            << simResults.stateProgression.at(i).plant(z)*1000 << ","
            << simResults.stateProgression.at(i).plant(phi)*180/M_PI << ","
            << -simResults.stateProgression.at(i).plant(theta)*180/M_PI << ","
            << capAngleTest(simResults.stateProgression.at(i).plant(psi)*180/M_PI) << ","
            << simResults.stateProgression.at(i).plant(xdot)*1000 << ","
            << simResults.stateProgression.at(i).plant(ydot)*1000 << ","
            << simResults.stateProgression.at(i).plant(zdot)*1000 << ","
            << simResults.stateProgression.at(i).plant(p)*1000 << ","
            << simResults.stateProgression.at(i).plant(q)*1000 << ","
            << simResults.stateProgression.at(i).plant(r)*1000 << std::endl;
    }
}


void test_closestzbar(Logger & log) {
    std::array<double(*)(double), NUM_DIST_STATES> dist = {noDist, noDist, noDist, noDist, noDist, noDist};
    std::array<double(*)(double), NUM_REF_STATES> ref = {oneRef, oneRef, zeroRef, zeroRef};
    double finalTime = 500;
    double simTime = 1e-3;
    DroneTrajectory droneTrajectory(log, dist, ref, finalTime, simTime);
    std::chrono::time_point start = std::chrono::steady_clock::now();
    Eigen::Vector<double, NUM_STATES>  zbar = droneTrajectory.closestZBar(stateCloseToRoABoundary());
}

int main()
{
    Logger log("./build/log.txt");
    // testERAAlgo(log);
    // testTrajSens(log);
    // testSim(log);
    // test_vp(log);
    // test_closestzbar(log);
    // test_param(log);
    // test_sysid(log);
    // test_square_traj(log);
    std::cout << ":D" << std::endl;
    return 0;
}