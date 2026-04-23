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
        
        for(int t = 0; t < 1; t++)
        {
            plus << plusSimResults.stateProgression.at(t).plant, plusSimResults.stateProgression.at(t).alge;
            minus << minusSimResults.stateProgression.at(t).plant, minusSimResults.stateProgression.at(t).alge;
            // m_logger << "index: " << i << std::endl;
            // m_logger << plus << std::endl;
            // m_logger << minus << std::endl;
            if (i < NUM_PLANT_STATES){
               ts.at(t).dxdwo.row(i) = 1/(2*delta)*(plus-minus).transpose();
            } else if ( i < NUM_PLANT_STATES + NUM_Z_STATES) {
                ts.at(t).dzdwo.row(i-NUM_PLANT_STATES) = 1/(2*delta)*(plus-minus).transpose();
            } else {
                ts.at(t).dydwo.row(i-NUM_PLANT_STATES-NUM_Z_STATES) = 1/(2*delta)*(plus-minus).transpose();
            }
        }
    }
    return ts;
}

