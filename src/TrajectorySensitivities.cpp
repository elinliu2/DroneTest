#include "DroneTrajectory.h"
#include <chrono>

std::vector<dwdwo>  DroneTrajectory::trajSens(SimResults const & simResults)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    const int iterations = simResults.time.size();
    std::vector<dwdwo> ts;
    ts.reserve(iterations);
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> initial_dxdwo;
    initial_dxdwo.setZero();
    initial_dxdwo.block(0, 0, NUM_PLANT_STATES, NUM_PLANT_STATES).setIdentity();
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> initial_dzdwo;
    initial_dzdwo.setZero();
    initial_dzdwo.block(0, NUM_PLANT_STATES, NUM_Z_STATES, NUM_Z_STATES).setIdentity();
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> initial_dydwo;
    initial_dydwo.setZero();
    initial_dydwo.block(0, NUM_PLANT_STATES+NUM_Z_STATES, NUM_Y_STATES, NUM_Y_STATES).setIdentity();
    ts.emplace_back(initial_dxdwo, initial_dzdwo, initial_dydwo);

    Eigen::SparseMatrix<double> I(NUM_PLANT_STATES, NUM_PLANT_STATES);
    I.setIdentity();

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    // Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> A;
    // Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> B;
    // Eigen::PartialPivLU<Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>> solver;
    
    // std::chrono::time_point elapsed0 = std::chrono::steady_clock::now();
    // m_logger << "Elapsed Time 0: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed0 - start).count() << " us" << std::endl;
    // m_logger << "test0" << std::endl; 

    // iterating trajectory sensitivity
    for(int i = 1; i < 2; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];
        // dwdwo
        dwdwo curr = ts[i-1];
        dxdwo = curr.dxdwo;
        dzdwo = curr.dzdwo;
        // std::chrono::time_point elapsed1 = std::chrono::steady_clock::now();
        // m_logger << "test1" << std::endl;

        // dxdwo_plus
        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::SparseMatrix<double> dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}); 
        Eigen::SparseMatrix<double> dfdz_curr = dfdz(simResults.stateProgression[i-1]);   
        // Eigen::SparseMatrix<double> dfdz_plus; 
        // Eigen::SparseMatrix<double> dfdz_curr;   
        // dfdz_plus.setZero();
        // dfdz_curr.setZero();
        
        // A = I - (timestep / 2.0) * dfdx_plus;
        // B = dxdwo + (timestep / 2.0)*(dfdx_curr*dxdwo + 2*dfdz_curr*dzdwo);
        // solver.compute(A);
        // dxdwo_plus.noalias() = solver.solve(B);
        dxdwo_plus = (I - (timestep / 2.0) * dfdx_plus).toDense().inverse() * (dxdwo + (timestep / 2.0)*(dfdx_curr*dxdwo + (dfdz_plus+dfdz_curr)*dzdwo));

        // dzdwo for 1 to n
        Eigen::SparseMatrix<double> dhdx_plus = dhdxPlus(simResults.stateProgression[i], timestep);
        dhdz_test(simResults.stateProgression[i], simResults.stateProgression[i-1], simResults.time[i], timestep);
        Eigen::SparseMatrix<double> dhdx_curr = dhdxCurr(timestep);
        Eigen::SparseMatrix<double> dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep); 
        Eigen::SparseMatrix<double> dhdz_curr = dhdzCurr(timestep); 
        Eigen::SparseMatrix<double> dhdy_plus = dhdy();
        Eigen::SparseMatrix<double> dgdz_plus = dgdz(simResults.stateProgression[i], timestep);
        dzdwo_plus = dhdx_plus * dxdwo_plus + dhdx_curr * dxdwo + dhdz_curr * dzdwo;
        // std::chrono::time_point elapsed3 = std::chrono::steady_clock::now();
        // m_logger << "test3" << std::endl;

        for(int i = 0; i < desThrust; i ++){
            Eigen::Matrix<double, 1, NUM_STATES> tmp;
            tmp.noalias() = dhdz_plus.block(i, 0, NUM_STATES, i) * dzdwo_plus.topRows(i);
            dzdwo_plus.row(i) += tmp;
        }
        // std::chrono::time_point elapsed4 = std::chrono::steady_clock::now();
        // m_logger << "test4" << std::endl;

        // dydwo
        dydwo_plus = dgdz_plus * dzdwo_plus;
        // std::chrono::time_point elapsed5 = std::chrono::steady_clock::now();
        // m_logger << "test5" << std::endl;

        // dzdwo 
        dzdwo_plus += dhdy_plus * dydwo_plus;
        for(int i = eiphi; i < NUM_Z_STATES; i ++){
            Eigen::Matrix<double, 1, NUM_STATES> tmp;
            tmp.noalias() = dhdz_plus.block(i, 0, NUM_STATES, i) * dzdwo_plus.topRows(i);
            dzdwo_plus.row(i) += tmp;
        }
        // std::chrono::time_point elapsed6 = std::chrono::steady_clock::now();
        // m_logger << "test6" << std::endl;

        ts.emplace_back(dxdwo_plus, dzdwo_plus, dydwo_plus);
        // std::chrono::time_point elapsed7 = std::chrono::steady_clock::now();
        // m_logger << "test7" << std::endl;

        // m_logger << "Elapsed Time 1: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed1 - elapsed0).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 2: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed2 - elapsed1).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 3: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed3 - elapsed2).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 4: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed4 - elapsed3).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 5: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed5 - elapsed4).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 6: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed6 - elapsed5).count() << " us" << std::endl;
        // m_logger << "Elapsed Time 7: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed7 - elapsed6).count() << " us" << std::endl;
    }

    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    m_logger << "Elapsed Time: " << elapsed.count() << " us" << std::endl;

    return ts;
}