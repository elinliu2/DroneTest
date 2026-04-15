#include "DroneTrajectory.h"

std::vector<Eigen::Matrix<double, NUM_STATES, NUM_STATES>>  DroneTrajectory::trajSens(SimResults simResults)
{
    const int iterations = simResults.time.size();
    std::vector<Eigen::Matrix<double, NUM_STATES, NUM_STATES>> ts;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(NUM_PLANT_STATES, NUM_PLANT_STATES);
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> initial_ts;
    initial_ts << I , Eigen::Matrix<double, NUM_PLANT_STATES, NUM_ALGE_STATES>::Zero(), Eigen::Matrix<double, NUM_ALGE_STATES, NUM_STATES>::Zero();
    ts.push_back(initial_ts);
    
    // iterating trajectory sensitivity
    for(int i = 1; i < iterations; i++)
    {
        double timestep = simResults.time.at(i) - simResults.time.at(i-1);

        // dwdwo
        Eigen::Matrix<double, NUM_STATES, NUM_STATES> dwdwo = ts.at(i-1);
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo = dwdwo.topRows(NUM_PLANT_STATES);
        Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES>dzdwo = dwdwo.middleRows(NUM_PLANT_STATES, NUM_Z_STATES);
        Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo = dwdwo.bottomRows(NUM_Y_STATES);

        // dxdwo_plus
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus = dfdx(simResults.stateProgression.at(i));
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr = dfdx(simResults.stateProgression.at(i-1));
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_curr = dfdz(simResults.stateProgression.at(i-1));        
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus = (I-timestep/2*dfdx_plus).inverse()*(dxdwo + timestep/2*(dfdx_curr * dxdwo + 2*dfdz_curr*dzdwo));
        Eigen::Matrix<double, 2*NUM_PLANT_STATES, NUM_STATES> dxdwo_plus_curr(2*NUM_PLANT_STATES, NUM_STATES);
        dxdwo_plus_curr << dxdwo_plus, dxdwo;
        
        // dzdwo for 1 to n
        Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_PLANT_STATES> dhdx_plus_curr = dhdx(simResults.stateProgression.at(i), timestep);
        Eigen::Matrix<double, NUM_Z_STATES, 2*NUM_Z_STATES> dhdz_plus_curr = dhdz(simResults.stateProgression.at(i), timestep); 
        Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> dhdy_plus = dhdy();
        Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> dgdz_plus = dgdz(simResults.stateProgression.at(i));
        Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus = Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES>::Zero();
        dzdwo_plus += dhdx_plus_curr * dxdwo_plus_curr;
        Eigen::Matrix<double, 2*NUM_Z_STATES, NUM_STATES> dzdwo_plus_curr;
        dzdwo_plus_curr << dzdwo_plus, dzdwo;

        for(int i = 0; i < desThrust; i ++){
            dzdwo_plus_curr.row(i) += dhdz_plus_curr.row(i) * dzdwo_plus_curr;
        }
        dzdwo_plus = dzdwo_plus_curr.topRows(NUM_Z_STATES);
        
        // dydwo
        Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo_plus = dgdz_plus * dzdwo_plus;

        // dzdwo 
        dzdwo_plus_curr.topRows(NUM_Z_STATES) += dhdy_plus * dydwo_plus;
        for(int i = epphi; i < NUM_Z_STATES; i ++){
            dzdwo_plus_curr.row(i) += dhdz_plus_curr.row(i) * dzdwo_plus_curr;
        }
        dzdwo_plus = dzdwo_plus_curr.topRows(NUM_Z_STATES);
        
        Eigen::Matrix<double, NUM_STATES, NUM_STATES> dwdwo_plus;
        dwdwo_plus << dxdwo_plus, dzdwo_plus_curr.topRows(NUM_Z_STATES), dydwo;
        ts.push_back(dwdwo_plus);
    }
    return ts;
}