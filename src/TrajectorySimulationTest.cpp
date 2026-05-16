#include "DroneTrajectory.h"

#include <thread>

std::vector<dwdwo> DroneTrajectory::trajSensTest(SystemState initialState)
{
    // std::chrono::time_point start = std::chrono::steady_clock::now();
    double delta = 1e-5;
    SimResults simResults = Trajectory(initialState);
    int numIterations = simResults.time.size();
    std::vector<dwdwo> ts(numIterations);
    
    Eigen::Vector<double, NUM_STATES> plus;
    Eigen::Vector<double, NUM_STATES> minus;
            
    for(int i = 0; i < NUM_STATES; i++) {
        SystemState testPlusState = initialState;
        SystemState testMinusState = initialState;
        if (i < NUM_PLANT_STATES){
            testPlusState.plant(i) += delta;
            testMinusState.plant(i) -= delta;
        } else if ( i < NUM_PLANT_STATES + NUM_Z_STATES) {
            testPlusState.alge(i-NUM_PLANT_STATES) += delta;
            testMinusState.alge(i-NUM_PLANT_STATES) -= delta;
        } else {
            testPlusState.alge(i-NUM_PLANT_STATES) += delta;
            testMinusState.alge(i-NUM_PLANT_STATES) -= delta;
        }
        SimResults plusSimResults = Trajectory(testPlusState);
        SimResults minusSimResults = Trajectory(testMinusState);
        
        for(int t = 0; t < numIterations; t++)
        {
            plus << plusSimResults.stateProgression.at(t).plant, plusSimResults.stateProgression.at(t).alge;
            minus << minusSimResults.stateProgression.at(t).plant, minusSimResults.stateProgression.at(t).alge;
            Eigen::Vector<double, NUM_STATES> delta_Trajsens = 1/(2*delta)*(plus-minus);
            ts.at(t).dxdwo.col(i) = delta_Trajsens.segment(0, NUM_PLANT_STATES);
            ts.at(t).dzdwo.col(i) = delta_Trajsens.segment(NUM_PLANT_STATES, NUM_Z_STATES);
            ts.at(t).dydwo.col(i) = delta_Trajsens.segment(NUM_PLANT_STATES+NUM_Z_STATES, NUM_Y_STATES);
        }
    }
    // std::chrono::time_point end = std::chrono::steady_clock::now();
    // std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // m_logger << "Elapsed Time test: " << elapsed.count() << " us" << std::endl;
    return ts;
}

std::vector<dwdp> DroneTrajectory::trajSensParamTest(SystemState initialState)
{

    double delta = 1e-5;
    SimResults simResults = Trajectory(initialState);
    int numIterations = simResults.time.size();
    std::vector<dwdp> ts_p(numIterations);
    
    Eigen::Vector<double, NUM_STATES> plus;
    Eigen::Vector<double, NUM_STATES> minus;
    
    SimResults plusSimResults;
    SimResults minusSimResults;
    std::array<PIDParameters, NUM_PIDS> og_params = m_ctrlParams;
    
    for(int i = 0; i < NUM_PARAMETERS; i++) {
        SystemState testPlusState = initialState;
        SystemState testMinusState = initialState;
         
        m_ctrlParams.at(i).kp = og_params.at(i).kp + delta;
        plusSimResults = Trajectory(initialState);
        m_ctrlParams.at(i).kp = og_params.at(i).kp - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_kp = trajSens(minusSimResults);
        m_ctrlParams.at(i).kp = og_params.at(i).kp;    

        m_ctrlParams.at(i).ki = og_params.at(i).ki + delta;
        plusSimResults = Trajectory(initialState);
        std::vector<dwdwo> plus_ki = trajSens(plusSimResults);
        m_ctrlParams.at(i).ki = og_params.at(i).ki - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_ki = trajSens(minusSimResults);
        m_ctrlParams.at(i).ki = og_params.at(i).ki;    

        m_ctrlParams.at(i).kd = og_params.at(i).kd + delta;
        plusSimResults = Trajectory(initialState);
        std::vector<dwdwo> plus_kd = trajSens(plusSimResults);
        m_ctrlParams.at(i).kd = og_params.at(i).kd - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_kd = trajSens(minusSimResults);   
        m_ctrlParams.at(i).kd = og_params.at(i).kd;     
        
    }
    return ts_p;
}

std::vector<d2wdwo2> DroneTrajectory::secondOrdertrajSensTest(SystemState initialState)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();
    double delta = 1e-5;
    SimResults simResults = Trajectory(initialState);
    int numIterations = simResults.time.size();
    std::vector<d2wdwo2> ts2(numIterations);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ;

    std::vector<std::vector<dwdwo>> plusses;
    std::vector<std::vector<dwdwo>> minuses;
            
    for(int i = 0; i < NUM_STATES; i++) {
        m_logger << i << std::endl;
        SystemState testPlusState = initialState;
        SystemState testMinusState = initialState;
        if (i < NUM_PLANT_STATES){
            testPlusState.plant(i) += delta;
            testMinusState.plant(i) -= delta;
        } else if ( i < NUM_PLANT_STATES + NUM_Z_STATES) {
            testPlusState.alge(i-NUM_PLANT_STATES) += delta;
            testMinusState.alge(i-NUM_PLANT_STATES) -= delta;
        } else {
            testPlusState.alge(i-NUM_PLANT_STATES) += delta;
            testMinusState.alge(i-NUM_PLANT_STATES) -= delta;
        }

        SimResults plusSimResults = Trajectory(testPlusState);
        SimResults minusSimResults = Trajectory(testMinusState);

        plusses.push_back(trajSens(plusSimResults));
        minuses.push_back(trajSens(minusSimResults));
    }

    for(int t = 0; t < numIterations; t++)
    {
        for(int i = 0; i < NUM_STATES; i++){
            resultX = (plusses.at(i).at(t).dxdwo - minuses.at(i).at(t).dxdwo) / (2.0 * delta);
            resultY = (plusses.at(i).at(t).dydwo - minuses.at(i).at(t).dydwo) / (2.0 * delta);
            resultZ = (plusses.at(i).at(t).dzdwo - minuses.at(i).at(t).dzdwo) / (2.0 * delta);

            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_tensor(resultX.data(), NUM_PLANT_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_tensor(resultY.data(), NUM_Y_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_tensor(resultZ.data(), NUM_Z_STATES, NUM_STATES);

            ts2.at(t).d2xdwo2.chip(i, 2) = resultX_tensor;
            ts2.at(t).d2ydwo2.chip(i, 2) = resultY_tensor;
            ts2.at(t).d2zdwo2.chip(i, 2) = resultZ_tensor;
        }
    }

    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    m_logger << "Elapsed Time test: " << elapsed.count() << " us" << std::endl;
    return ts2;
}

std::vector<d2wdwodp> DroneTrajectory::secondOrdertrajSensParamsTest(SystemState initialState)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();
    double delta = 1e-5;
    SimResults simResults = Trajectory(initialState);
    SimResults plusSimResults;
    SimResults minusSimResults;
    int numIterations = simResults.time.size();
    std::vector<d2wdwodp> ts2(numIterations);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_kp;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_kp;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_kp;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_ki;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_ki;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_ki;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_kd;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_kd;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_kd;

    std::array<PIDParameters, NUM_PIDS> og_params = m_ctrlParams;
            
    for(int i = 0; i < NUM_PIDS; i++) {
        m_ctrlParams.at(i).kp = og_params.at(i).kp + delta;
        plusSimResults = Trajectory(initialState);
        std::vector<dwdwo> plus_kp = trajSens(plusSimResults);
        m_ctrlParams.at(i).kp = og_params.at(i).kp - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_kp = trajSens(minusSimResults);
        m_ctrlParams.at(i).kp = og_params.at(i).kp;    

        m_ctrlParams.at(i).ki = og_params.at(i).ki + delta;
        plusSimResults = Trajectory(initialState);
        std::vector<dwdwo> plus_ki = trajSens(plusSimResults);
        m_ctrlParams.at(i).ki = og_params.at(i).ki - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_ki = trajSens(minusSimResults);
        m_ctrlParams.at(i).ki = og_params.at(i).ki;    

        m_ctrlParams.at(i).kd = og_params.at(i).kd + delta;
        plusSimResults = Trajectory(initialState);
        std::vector<dwdwo> plus_kd = trajSens(plusSimResults);
        m_ctrlParams.at(i).kd = og_params.at(i).kd - delta;
        minusSimResults = Trajectory(initialState);
        std::vector<dwdwo> minus_kd = trajSens(minusSimResults);
        m_ctrlParams.at(i).kd = og_params.at(i).kd;    

        for(int t = 0; t < numIterations; t++)
        {
            resultX_kp = (plus_kp.at(t).dxdwo - minus_kp.at(t).dxdwo) / (2.0 * delta);
            resultY_kp = (plus_kp.at(t).dydwo - minus_kp.at(t).dydwo) / (2.0 * delta);
            resultZ_kp = (plus_kp.at(t).dzdwo - minus_kp.at(t).dzdwo) / (2.0 * delta);

            resultX_ki = (plus_ki.at(t).dxdwo - minus_ki.at(t).dxdwo) / (2.0 * delta);
            resultY_ki = (plus_ki.at(t).dydwo - minus_ki.at(t).dydwo) / (2.0 * delta);
            resultZ_ki = (plus_ki.at(t).dzdwo - minus_ki.at(t).dzdwo) / (2.0 * delta);

            resultX_kd = (plus_kd.at(t).dxdwo - minus_kd.at(t).dxdwo) / (2.0 * delta);
            resultY_kd = (plus_kd.at(t).dydwo - minus_kd.at(t).dydwo) / (2.0 * delta);
            resultZ_kd = (plus_kd.at(t).dzdwo - minus_kd.at(t).dzdwo) / (2.0 * delta);

            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_kp_tensor(resultX_kp.data(), NUM_PLANT_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_kp_tensor(resultY_kp.data(), NUM_Y_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_kp_tensor(resultZ_kp.data(), NUM_Z_STATES, NUM_STATES);

            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_ki_tensor(resultX_ki.data(), NUM_PLANT_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_ki_tensor(resultY_ki.data(), NUM_Y_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_ki_tensor(resultZ_ki.data(), NUM_Z_STATES, NUM_STATES);


            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_kd_tensor(resultX_kd.data(), NUM_PLANT_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_kd_tensor(resultY_kd.data(), NUM_Y_STATES, NUM_STATES);
            Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_kd_tensor(resultZ_kd.data(), NUM_Z_STATES, NUM_STATES);

            ts2.at(t).d2xdwodp.chip(i*NUM_PID_STATES, 2) = resultX_kp_tensor;
            ts2.at(t).d2ydwodp.chip(i*NUM_PID_STATES, 2) = resultY_kp_tensor;
            ts2.at(t).d2zdwodp.chip(i*NUM_PID_STATES, 2) = resultZ_kp_tensor;

            ts2.at(t).d2xdwodp.chip(i*NUM_PID_STATES + 1, 2) = resultX_ki_tensor;
            ts2.at(t).d2ydwodp.chip(i*NUM_PID_STATES + 1, 2) = resultY_ki_tensor;
            ts2.at(t).d2zdwodp.chip(i*NUM_PID_STATES + 1, 2) = resultZ_ki_tensor;

            ts2.at(t).d2xdwodp.chip(i*NUM_PID_STATES + 2, 2) = resultX_kd_tensor;
            ts2.at(t).d2ydwodp.chip(i*NUM_PID_STATES + 2, 2) = resultY_kd_tensor;
            ts2.at(t).d2zdwodp.chip(i*NUM_PID_STATES + 2, 2) = resultZ_kd_tensor;
        }
    }
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    m_logger << "Elapsed Time test: " << elapsed.count() << " us" << std::endl;
    return ts2;
}

void calcStatesDG(DroneTrajectory dt, int threadIndex, int numStates, d2wdwo2 & trajSenswo2, SystemState initialState, int tp)
{
    double delta = 1e-5;
    for(int i = 0; i < numStates; i++) {
        SystemState testPlusState = initialState;
        SystemState testMinusState = initialState;
        if (i+threadIndex * numStates < NUM_PLANT_STATES){
            testPlusState.plant(i+threadIndex * numStates) += delta;
            testMinusState.plant(i+threadIndex * numStates) -= delta;
        } else {
            testPlusState.alge(i+threadIndex * numStates-NUM_PLANT_STATES) += delta;
            testMinusState.alge(i+threadIndex * numStates-NUM_PLANT_STATES) -= delta;
        }

        SimResults plusSimResults = dt.Trajectory(testPlusState, false);
        SimResults minusSimResults = dt.Trajectory(testMinusState, false);

        std::vector<dwdwo> plus_ts = dt.trajSens(plusSimResults);
        std::vector<dwdwo> minus_ts = dt.trajSens(minusSimResults);

        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX = (plus_ts.at(tp).dxdwo - minus_ts.at(tp).dxdwo) / (2.0 * delta);
        Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY = (plus_ts.at(tp).dydwo - minus_ts.at(tp).dydwo) / (2.0 * delta);
        Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ = (plus_ts.at(tp).dzdwo - minus_ts.at(tp).dzdwo) / (2.0 * delta);
        
        Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_tensor(resultX.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_tensor(resultY.data(), NUM_Y_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_tensor(resultZ.data(), NUM_Z_STATES, NUM_STATES);

        trajSenswo2.d2xdwo2.chip(i+threadIndex * numStates, 2) = resultX_tensor;
        trajSenswo2.d2ydwo2.chip(i+threadIndex * numStates, 2) = resultY_tensor;
        trajSenswo2.d2zdwo2.chip(i+threadIndex * numStates, 2) = resultZ_tensor;
    }
}

void calcController(DroneTrajectory dt, int threadIndex, d2wdwodp & trajSenswodp, std::array<PIDParameters, NUM_PIDS> og_params, SystemState initialState, int tp)
{
    SimResults plusSimResults;
    SimResults minusSimResults;
    double delta = 1e-5;
  
    dt.m_ctrlParams.at(threadIndex).kp = og_params.at(threadIndex).kp + delta;
    plusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> plus_kp = dt.trajSens(plusSimResults);
    dt.m_ctrlParams.at(threadIndex).kp = og_params.at(threadIndex).kp - delta;
    minusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> minus_kp = dt.trajSens(minusSimResults);
    dt.m_ctrlParams.at(threadIndex).kp = og_params.at(threadIndex).kp;

    dt.m_ctrlParams.at(threadIndex).ki = og_params.at(threadIndex).ki + delta;
    plusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> plus_ki = dt.trajSens(plusSimResults);
    dt.m_ctrlParams.at(threadIndex).ki = og_params.at(threadIndex).ki - delta;
    minusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> minus_ki = dt.trajSens(minusSimResults);
    dt.m_ctrlParams.at(threadIndex).ki = og_params.at(threadIndex).ki;

    dt.m_ctrlParams.at(threadIndex).kd = og_params.at(threadIndex).kd + delta;
    plusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> plus_kd = dt.trajSens(plusSimResults);
    dt.m_ctrlParams.at(threadIndex).kd = og_params.at(threadIndex).kd - delta;
    minusSimResults = dt.Trajectory(initialState, false);
    std::vector<dwdwo> minus_kd = dt.trajSens(minusSimResults);
    dt.m_ctrlParams.at(threadIndex).kd = og_params.at(threadIndex).kd;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_kp = (plus_kp.at(tp).dxdwo - minus_kp.at(tp).dxdwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_kp = (plus_kp.at(tp).dydwo - minus_kp.at(tp).dydwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_kp = (plus_kp.at(tp).dzdwo - minus_kp.at(tp).dzdwo) / (2.0 * delta);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_ki = (plus_ki.at(tp).dxdwo - minus_ki.at(tp).dxdwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_ki = (plus_ki.at(tp).dydwo - minus_ki.at(tp).dydwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_ki = (plus_ki.at(tp).dzdwo - minus_ki.at(tp).dzdwo) / (2.0 * delta);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> resultX_kd = (plus_kd.at(tp).dxdwo - minus_kd.at(tp).dxdwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> resultY_kd = (plus_kd.at(tp).dydwo - minus_kd.at(tp).dydwo) / (2.0 * delta);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> resultZ_kd = (plus_kd.at(tp).dzdwo - minus_kd.at(tp).dzdwo) / (2.0 * delta);

    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_kp_tensor(resultX_kp.data(), NUM_PLANT_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_kp_tensor(resultY_kp.data(), NUM_Y_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_kp_tensor(resultZ_kp.data(), NUM_Z_STATES, NUM_STATES);

    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_ki_tensor(resultX_ki.data(), NUM_PLANT_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_ki_tensor(resultY_ki.data(), NUM_Y_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_ki_tensor(resultZ_ki.data(), NUM_Z_STATES, NUM_STATES);

    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultX_kd_tensor(resultX_kd.data(), NUM_PLANT_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultY_kd_tensor(resultY_kd.data(), NUM_Y_STATES, NUM_STATES);
    Eigen::TensorMap<Eigen::Tensor<double, 2>> resultZ_kd_tensor(resultZ_kd.data(), NUM_Z_STATES, NUM_STATES);

    trajSenswodp.d2xdwodp.chip(threadIndex*NUM_PID_STATES, 2) = resultX_kp_tensor;
    trajSenswodp.d2ydwodp.chip(threadIndex*NUM_PID_STATES, 2) = resultY_kp_tensor;
    trajSenswodp.d2zdwodp.chip(threadIndex*NUM_PID_STATES, 2) = resultZ_kp_tensor;

    trajSenswodp.d2xdwodp.chip(threadIndex*NUM_PID_STATES + 1, 2) = resultX_ki_tensor;
    trajSenswodp.d2ydwodp.chip(threadIndex*NUM_PID_STATES + 1, 2) = resultY_ki_tensor;
    trajSenswodp.d2zdwodp.chip(threadIndex*NUM_PID_STATES + 1, 2) = resultZ_ki_tensor;

    trajSenswodp.d2xdwodp.chip(threadIndex*NUM_PID_STATES + 2, 2) = resultX_kd_tensor;
    trajSenswodp.d2ydwodp.chip(threadIndex*NUM_PID_STATES + 2, 2) = resultY_kd_tensor;
    trajSenswodp.d2zdwodp.chip(threadIndex*NUM_PID_STATES + 2, 2) = resultZ_kd_tensor;
}


Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> DroneTrajectory::calc_dG_test(SystemState initialState, dwdwo ts, G_tp gtp, double endtime)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> dwdwo_matrix;
    dwdwo_matrix << ts.dxdwo, ts.dzdwo, ts.dydwo;

    // I don't like this either but I don't have time to fix it right now
    // Really this should just be passed in
    double OG_finalTime = m_finalTime;
    m_finalTime = endtime;

    d2wdwo2 trajSenswo2;
    int num_threads = 12;
    std::vector<std::thread> threads;
    for (int threadIndex = 0; threadIndex < num_threads; threadIndex++)
    {
        threads.push_back(std::thread(calcStatesDG, *this, threadIndex, NUM_STATES/num_threads, std::ref(trajSenswo2), initialState, gtp.tp));
    }

    for (auto& t : threads) {
        t.join();
    }
    
    threads.clear();
    d2wdwodp trajSenswodp;
    for (int threadIndex = 0; threadIndex < num_threads; threadIndex++)
    {
        threads.push_back(std::thread(calcController, *this, threadIndex, std::ref(trajSenswodp), m_ctrlParams, initialState, gtp.tp));
    }

    for (auto& t : threads) {
        t.join();
    }
    
    Eigen::Tensor<double, 3> tmp1 = trajSenswo2.d2xdwo2.concatenate(trajSenswo2.d2ydwo2, 0);
    Eigen::Tensor<double, 3> secondOrderTraj = tmp1.concatenate(trajSenswo2.d2zdwo2, 0);
    Eigen::Tensor<double, 3> tmp2 = trajSenswodp.d2xdwodp.concatenate(trajSenswodp.d2ydwodp, 0);
    Eigen::Tensor<double, 3> secondOrderTrajParams = tmp2.concatenate(trajSenswodp.d2zdwodp, 0);

    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS>::Zero();
    for(int j = 0; j < NUM_STATES; j++)
    {
        for(int i = 0; i < NUM_STATES; i++)
        {
            Eigen::Vector<double, NUM_STATES> signed_dw = dwdwo_matrix.col(i).array().sign();
            Eigen::Tensor<double, 1> secondOrderTrajChip = secondOrderTraj.chip(j, 2).chip(i, 1);
            Eigen::Map<Eigen::Vector<double, NUM_STATES>> secondOrderTrajChipVec(secondOrderTrajChip.data(), secondOrderTrajChip.size());
            dG(j) += signed_dw.transpose() * secondOrderTrajChipVec;
        }
    }

    for(int j = 0; j < NUM_PARAMETERS; j++)
    {
        for(int i = 0; i < NUM_STATES; i++)
        {
            Eigen::Vector<double, NUM_STATES> signed_dw = dwdwo_matrix.col(i).array().sign();
            Eigen::Tensor<double, 1> secondOrderTrajChip = secondOrderTrajParams.chip(j, 2).chip(i, 1);
            Eigen::Map<Eigen::Vector<double, NUM_STATES>> secondOrderTrajChipVec(secondOrderTrajChip.data(), secondOrderTrajChip.size());
            dG(NUM_STATES + j) += signed_dw.transpose() * secondOrderTrajChipVec;
        }
    }
    
    dG = -1 * std::pow(gtp.G, 2) * dG;
    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    m_logger << "Elapsed Time dG: " << elapsed.count() << " us" << std::endl;
    m_finalTime = OG_finalTime;
    return dG;
}

void DroneTrajectory::dfdx_test(SystemState initialState)
{
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> exact_dfdx = dfdx(initialState);
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> delta_dfdx;

    SystemState plus; 
    SystemState minus; 
    double delta = 1e-5;
    for(int i = 0; i < NUM_PLANT_STATES; i++) 
    {
        plus = initialState;
        plus.plant(i) += delta;
        minus = initialState;
        minus.plant(i) -= delta;

        Eigen::Vector<double, NUM_PLANT_STATES> plus_f = f(plus, 0);
        Eigen::Vector<double, NUM_PLANT_STATES> minus_f = f(minus, 0);
        delta_dfdx.col(i) =  1/(2*delta)*(plus_f-minus_f);
    }
    m_logger << "dfdx max diff: " << std::max((exact_dfdx - delta_dfdx).maxCoeff(), (delta_dfdx - exact_dfdx).maxCoeff()) << std::endl;
}

void DroneTrajectory::dfdz_test(SystemState initialState)
{
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> exact_dfdz = dfdz(initialState);
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> delta_dfdz;

    SystemState plus; 
    SystemState minus; 
    double delta = 1e-5;
    for(int i = 0; i < NUM_Z_STATES; i++) 
    {
        plus = initialState;
        plus.alge(i) += delta;
        minus = initialState;
        minus.alge(i) -= delta;

        Eigen::Vector<double, NUM_PLANT_STATES> plus_f = f(plus, 0);
        Eigen::Vector<double, NUM_PLANT_STATES> minus_f = f(minus, 0);
        delta_dfdz.col(i) =  1/(2*delta)*(plus_f-minus_f);
    }
    // m_logger << "exact dfdz" << std::endl;
    // m_logger << exact_dfdz << std::endl;
    // m_logger << "delta dfdz" << std::endl;
    // m_logger << delta_dfdz << std::endl;

    m_logger << "dfdz max diff: " << std::max((exact_dfdz - delta_dfdz).maxCoeff(), (delta_dfdz - exact_dfdz).maxCoeff()) << std::endl;

}

void DroneTrajectory::dhdx_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> exact_dhdxPlus = dhdxPlus(currState, time, timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> exact_dhdxCurr = dhdxCurr(prevState, timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> delta_dhdxPlus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> delta_dhdxCurr;
    
    SystemState plus; 
    SystemState minus; 
    double delta = 1e-5;
    for(int i = 0; i < NUM_PLANT_STATES; i++) 
    {
        plus = currState;
        plus.plant(i) += delta;
        minus = currState;
        minus.plant(i) -= delta;
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = CascadedPIDController(plus.plant, prevState.plant, currState.alge, time, timestep);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = CascadedPIDController(minus.plant, prevState.plant, currState.alge, time, timestep);
        delta_dhdxPlus.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(0, NUM_Z_STATES);
    }

    // for(int i = 0; i < NUM_Z_STATES; i++)
    // {
    //     for(int j = 0; j < NUM_PLANT_STATES; j++)
    //     {
    //         if (std::abs(exact_dhdxPlus(i, j) - delta_dhdxPlus(i, j)) > 2e-5)
    //         {
    //             m_logger << "dhdx_diff: index - " << i << " " << j << " diff: " << exact_dhdxPlus(i, j) - delta_dhdxPlus(i, j) << std::endl;
    //             m_logger << "exact_dhdxPlus " << exact_dhdxPlus(i, j) << " delta_dhdxPlus " << delta_dhdxPlus(i, j) << std::endl;
    //         }
    //     }
    // }

    for(int i = 0; i < NUM_PLANT_STATES; i++) 
    {
        plus = prevState;
        plus.plant(i) += delta;
        minus = prevState;
        minus.plant(i) -= delta;
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = CascadedPIDController(currState.plant, plus.plant, currState.alge, time, timestep);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = CascadedPIDController(currState.plant, minus.plant, currState.alge, time, timestep);
        delta_dhdxCurr.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(0, NUM_Z_STATES);
    }

    // for(int i = 0; i < NUM_Z_STATES; i++)
    // {
    //     for(int j = 0; j < NUM_PLANT_STATES; j++)
    //     {
    //         if (std::abs(exact_dhdxCurr(i, j) - delta_dhdxCurr(i, j)) > 2e-5)
    //         {
    //             m_logger << "dhdx_diff: index - " << i << " " << j << " diff: " << exact_dhdxCurr(i, j) - delta_dhdxCurr(i, j) << std::endl;
    //             m_logger << "exact_dhdxPlus " << exact_dhdxCurr(i, j) << " delta_dhdxPlus " << delta_dhdxCurr(i, j) << std::endl;
    //         }
    //     }
    // }

    m_logger << "dhdxPlus max diff: " << std::max((exact_dhdxPlus - delta_dhdxPlus).maxCoeff(), (delta_dhdxPlus - exact_dhdxPlus).maxCoeff()) << std::endl;
    m_logger << "dhdxCurr max diff: " << std::max((exact_dhdxCurr - delta_dhdxCurr).maxCoeff(), (exact_dhdxCurr - delta_dhdxCurr).maxCoeff()) << std::endl;

}

void DroneTrajectory::dhdz_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> exact_dhdzPlus = dhdzPlus(currState, timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> exact_dhdzCurr = dhdzCurr();
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> delta_dhdzPlus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> delta_dhdzCurr;
    
    double delta = 1e-5;
    for(int i = 0; i < NUM_Z_STATES; i++) 
    {
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, i, delta);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, i, -delta);
        delta_dhdzPlus.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(0, NUM_Z_STATES);
    }

    // m_logger << "exact dhdzPlus" << std::endl;
    // m_logger << exact_dhdzPlus << std::endl;
    // m_logger << "delta dhdzPlus" << std::endl;
    // m_logger << delta_dhdzPlus << std::endl;

    // for(int i = 0; i < NUM_Z_STATES; i++)
    // {
    //     for(int j = 0; j < NUM_Z_STATES; j++)
    //     {
    //         if (std::abs(exact_dhdzPlus(i, j) - delta_dhdzPlus(i, j)) > 2e-5)
    //         {
    //             m_logger << "dhdz_diff: index - " << i << " " << j << " diff: " << exact_dhdzPlus(i, j) - delta_dhdzPlus(i, j) << std::endl;
    //             m_logger << "exact_dhdzPlus " << exact_dhdzPlus(i, j) << " delta_dhdzPlus " << delta_dhdzPlus(i, j) << std::endl;
    //         }
    //     }
    // }

    m_logger << "dhdzPlus max diff: " << std::max((exact_dhdzPlus - delta_dhdzPlus).maxCoeff(), (delta_dhdzPlus - exact_dhdzPlus).maxCoeff()) << std::endl;

    for(int i = 0; i < NUM_Z_STATES; i++) 
    {
        Eigen::Vector<double, NUM_ALGE_STATES> plus = currState.alge;
        Eigen::Vector<double, NUM_ALGE_STATES> minus = currState.alge;
        plus(i) += delta;
        minus(i) -= delta;
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = h(currState.plant, prevState.plant, plus, time, timestep, NUM_Z_STATES, 0);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = h(currState.plant, prevState.plant, minus, time, timestep, NUM_Z_STATES, 0);
        delta_dhdzCurr.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(0, NUM_Z_STATES);
    }

    for(int i = 0; i < NUM_Z_STATES; i++)
    {
        for(int j = 0; j < NUM_Z_STATES; j++)
        {
            if (std::abs(exact_dhdzCurr(i, j) - delta_dhdzCurr(i, j)) > 2e-5)
            {
                m_logger << "dhdz_diff: index - " << i << " " << j << " diff: " << exact_dhdzCurr(i, j) - delta_dhdzCurr(i, j) << std::endl;
                m_logger << "exact_dhdzCurr " << exact_dhdzCurr(i, j) << " delta_dhdzCurr " << delta_dhdzCurr(i, j) << std::endl;
            }
        }
    }
    m_logger << "dhdzCurr max diff: " << std::max((exact_dhdzCurr - delta_dhdzCurr).maxCoeff(), (delta_dhdzCurr - exact_dhdzCurr).maxCoeff()) << std::endl;
}

void DroneTrajectory::dgdz_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> exact_dgdz = dgdz(currState);
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> delta_dgdz;
    
    double delta = 1e-5;
    for(int i = 0; i < NUM_Z_STATES; i++) 
    {
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, i, delta);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, i, -delta);
        delta_dgdz.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(NUM_Z_STATES, NUM_Y_STATES);
    }
    for(int i = 0; i < NUM_Y_STATES; i++)
    {
        for(int j = 0; j < NUM_Z_STATES; j++)
        {
            if (std::abs(exact_dgdz(i, j) - delta_dgdz(i, j)) > 2e-5)
            {
                m_logger << "dgdz: index - " << i << " " << j << " diff: " << exact_dgdz(i, j) - delta_dgdz(i, j) << std::endl;
                m_logger << "exact_dgdz " << exact_dgdz(i, j) << " delta_dgdz " << delta_dgdz(i, j) << std::endl;
            }
        }
    }
    m_logger << "dgdz max diff: " << std::max((exact_dgdz - delta_dgdz).maxCoeff(), (delta_dgdz - exact_dgdz).maxCoeff()) << std::endl;

}

void DroneTrajectory::dhdy_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> exact_dhdy = dhdy(timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> delta_dhdy;
    
    double delta = 1e-5;
    for(int i = 0; i < NUM_Y_STATES; i++) 
    {
        Eigen::Vector<double, NUM_ALGE_STATES> plus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, NUM_Z_STATES+i, delta);
        Eigen::Vector<double, NUM_ALGE_STATES> minus_h = h(currState.plant, prevState.plant, currState.alge, time, timestep, NUM_Z_STATES+i, -delta);
        delta_dhdy.col(i) =  1/(2*delta)*(plus_h-minus_h).segment(0, NUM_Z_STATES);
    }
    for(int i = 0; i < NUM_Z_STATES; i++)
    {
        for(int j = 0; j < NUM_Y_STATES; j++)
        {
            if (std::abs(exact_dhdy(i, j) - delta_dhdy(i, j)) > 2e-5)
            {
                m_logger << "dhdy: index - " << i << " " << j << " diff: " << exact_dhdy(i, j) - delta_dhdy(i, j) << std::endl;
                m_logger << "exact_dgdz " << exact_dhdy(i, j) << " delta_dgdz " << delta_dhdy(i, j) << std::endl;
            }
        }
    }
    m_logger << "dhdy max diff: " << std::max((exact_dhdy - delta_dhdy).maxCoeff(), (delta_dhdy - exact_dhdy).maxCoeff()) << std::endl;

}