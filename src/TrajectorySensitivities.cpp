#include "DroneTrajectory.h"
#include <chrono>

std::vector<dwdwo>  DroneTrajectory::trajSens(SimResults const & simResults)
{
    // std::chrono::time_point start = std::chrono::steady_clock::now();

    const int iterations = simResults.time.size();
    // m_logger << "iterations: " << iterations << std::endl;
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
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> A;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> B;
    Eigen::PartialPivLU<Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>> solver;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo_plus;

    // iterating trajectory sensitivity
    for(int i = 1; i < iterations; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];
        // dwdwo
        dwdwo curr = ts[i-1];
        dxdwo = curr.dxdwo;
        dzdwo = curr.dzdwo;

        // dxdwo_plus
        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::SparseMatrix<double> dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}); 
        Eigen::SparseMatrix<double> dfdz_curr = dfdz(simResults.stateProgression[i-1]);   
        
        A = I - (timestep / 2.0) * dfdx_plus;
        B = dxdwo + (timestep / 2.0)*(dfdx_curr*dxdwo + (dfdz_plus+dfdz_curr)*dzdwo);
        solver.compute(A);
        dxdwo_plus.noalias() = solver.solve(B);
        
        // dzdwo for 1 to n
        Eigen::SparseMatrix<double> dhdx_plus = dhdxPlus(simResults.stateProgression[i], simResults.time[i], timestep);
        Eigen::SparseMatrix<double> dhdx_curr = dhdxCurr(simResults.stateProgression[i-1], timestep);
        Eigen::SparseMatrix<double> dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep); 
        Eigen::SparseMatrix<double> dhdz_curr = dhdzCurr(); 
        Eigen::SparseMatrix<double> dhdy_plus = dhdy(timestep);
        Eigen::SparseMatrix<double> dgdz_plus = dgdz(simResults.stateProgression[i]);
        Eigen::SparseMatrix<double> dgdx_plus = dgdx(simResults.stateProgression[i]);

        dzdwo_plus = dhdx_plus * dxdwo_plus + dhdx_curr * dxdwo + dhdz_curr * dzdwo;

        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey ++){
            Eigen::Matrix<double, 1, NUM_STATES> tmp;
            tmp = dhdz_plus.block(zBeforey, 0, 1, zBeforey) * dzdwo_plus.topRows(zBeforey);
            dzdwo_plus.row(zBeforey) += tmp;
        }  
        
        // dydwo
        dydwo_plus = dgdz_plus * dzdwo_plus + dgdx_plus * dxdwo_plus;

        // dzdwo 
        dzdwo_plus += dhdy_plus * dydwo_plus;
        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::Matrix<double, 1, NUM_STATES> tmp;
            tmp = dhdz_plus.block(zAftery, 0, 1, zAftery) * dzdwo_plus.topRows(zAftery);
            dzdwo_plus.row(zAftery) += tmp;
        }
        
        ts.emplace_back(dxdwo_plus, dzdwo_plus, dydwo_plus);
    }

    // std::chrono::time_point end = std::chrono::steady_clock::now();
    // std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // m_logger << "Elapsed Time: " << elapsed.count() << " us" << std::endl;

    return ts;
}

std::vector<d2wdwo2>  DroneTrajectory::secondOrdertrajSens(SimResults const & simResults, std::vector<dwdwo> const & ts)
{
    // std::chrono::time_point start = std::chrono::steady_clock::now();

    const int iterations = simResults.time.size();
    m_logger << "iterations: " << iterations << std::endl;
    std::vector<d2wdwo2> ts2;
    ts2.reserve(iterations);

    Eigen::Tensor<double, 3> init_d2xdwo2;
    init_d2xdwo2.setZero();
    Eigen::Tensor<double, 3> init_d2zdwo2;
    init_d2zdwo2.setZero();
    Eigen::Tensor<double, 3> init_d2ydwo2;
    init_d2ydwo2.setZero();

    ts2.emplace_back(init_d2xdwo2, init_d2zdwo2, init_d2ydwo2);

    Eigen::Tensor<double, 3> d2xdwo2_curr;
    Eigen::Tensor<double, 3> d2zdwo2_curr;
    Eigen::Tensor<double, 3> d2xdwo2_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> I;
    I.setIdentity();


    for(int i = 1; i < iterations; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];

        // d2xdwo2
        d2xdwo2_curr = ts2[i-1].d2xdwo2;
        d2zdwo2_curr = ts2[i-1].d2zdwo2;
        dxdwo_curr = ts[i-1].dxdwo;
        dxdwo_plus = ts[i].dxdwo;

        Eigen::SparseMatrix<double> d2fdx2_plus = d2fdx2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::SparseMatrix<double> d2fdx2_curr = d2fdx2({simResults.stateProgression[i-1].plant, simResults.stateProgression[i-1].alge});
        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::SparseMatrix<double> dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}); 
        Eigen::SparseMatrix<double> dfdz_curr = dfdz(simResults.stateProgression[i-1]);   
        
        Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_STATES> a1_A = d2fdx2_plus*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 3>> a1_A_tensor(a1_A.data(), NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a1_B_tensor(dxdwo_plus.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(1, 1)};
        Eigen::Tensor<double, 3> a1 = a1_A_tensor.contract(a1_B_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3> a2 = dfdz_mult_d2zdwo2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}, d2zdwo2_curr);

        Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_STATES> a3_A = d2fdx2_curr*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 3>> a3_A_tensor(a3_A.data(), NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a3_B_tensor(dxdwo_curr.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::Tensor<double, 3> a3 = a3_A_tensor.contract(a3_B_tensor, contract_dims).eval();

        Eigen::TensorMap<Eigen::Tensor<double, 2>> a4_A_tensor(dfdx_curr.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3> a4 = a4_A_tensor.contract(d2xdwo2_plus, contract_dims).eval();

        Eigen::Tensor<double, 3> a5 = dfdz_mult_d2zdwo2(simResults.stateProgression[i-1], d2zdwo2_curr);

        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a6_A_matrix = (I - dfdx_plus).inverse();
        Eigen::Tensor<double, 3> a6_B_tensor = d2xdwo2_curr + timestep/2*(a1+a2+a3+a4+a5);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a6_A_tensor(a6_A_matrix.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        d2xdwo2_plus = a6_A_tensor.contract(a6_B_tensor, contract_dims).eval();
    }

    // std::chrono::time_point end = std::chrono::steady_clock::now();
    // std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // m_logger << "Elapsed Time: " << elapsed.count() << " us" << std::endl;

    return ts2;
}


G_tp DroneTrajectory::calc_G_tp(std::vector<dwdwo> trajSens)
{
    double G = 0;
    int index = 0;
    for(int i = 0; i < trajSens.size();i++)
    {
        double oneNormButAddUpAllElements = trajSens.at(i).dxdwo.cwiseAbs().sum() + trajSens.at(i).dydwo.cwiseAbs().sum() + trajSens.at(i).dzdwo.cwiseAbs().sum();
        if (oneNormButAddUpAllElements > G)
        {
            G = oneNormButAddUpAllElements;
            index = i;
        }
    }
    return{1/G, index};
}
