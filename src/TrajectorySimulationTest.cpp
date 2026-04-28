#include "DroneTrajectory.h"


std::vector<dwdwo> DroneTrajectory::trajSensTest(SystemState initialState)
{
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
    return ts;
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

    m_logger << "dfdx max diff: " << std::max((exact_dfdz - delta_dfdz).maxCoeff(), (delta_dfdz - exact_dfdz).maxCoeff()) << std::endl;

}

void DroneTrajectory::dhdx_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> exact_dhdxPlus = dhdxPlus(currState, time, timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> exact_dhdxCurr = dhdxCurr(currState, prevState, timestep);
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

    m_logger << "dhdxPlus max diff: " << std::max((exact_dhdxPlus - delta_dhdxPlus).maxCoeff(), (delta_dhdxPlus - exact_dhdxPlus).maxCoeff()) << std::endl;

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

    m_logger << "dhdxCurr max diff: " << std::max((exact_dhdxCurr - delta_dhdxCurr).maxCoeff(), (exact_dhdxCurr - delta_dhdxCurr).maxCoeff()) << std::endl;

}

void DroneTrajectory::dhdz_test(SystemState currState, SystemState prevState, double time, double timestep)
{
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> exact_dhdzPlus = dhdzPlus(currState, timestep);
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> exact_dhdzCurr = dhdzCurr(currState, timestep);
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

    for(int i = 0; i < NUM_Z_STATES; i++)
    {
        for(int j = 0; j < NUM_Z_STATES; j++)
        {
            if (std::abs(exact_dhdzPlus(i, j) - delta_dhdzPlus(i, j)) > 2e-5)
            {
                m_logger << "dhdz_diff: index - " << i << " " << j << " diff: " << exact_dhdzPlus(i, j) - delta_dhdzPlus(i, j) << std::endl;
                m_logger << "exact_dhdzPlus " << exact_dhdzPlus(i, j) << " delta_dhdzPlus " << delta_dhdzPlus(i, j) << std::endl;
            }
        }
    }

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