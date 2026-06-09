#include "DroneTrajectory.h"
#include <chrono>

std::vector<dwdwo> DroneTrajectory::trajSens(SimResults const & simResults) const
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

dwdwo DroneTrajectory::trajSens(SimResults const & simResults, int tp) const
{
    const int iterations = simResults.time.size();
    dwdwo ts; 
    
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> initial_dxdwo;
    initial_dxdwo.setZero();
    initial_dxdwo.block(0, 0, NUM_PLANT_STATES, NUM_PLANT_STATES).setIdentity();
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> initial_dzdwo;
    initial_dzdwo.setZero();
    initial_dzdwo.block(0, NUM_PLANT_STATES, NUM_Z_STATES, NUM_Z_STATES).setIdentity();
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> initial_dydwo;
    initial_dydwo.setZero();
    initial_dydwo.block(0, NUM_PLANT_STATES+NUM_Z_STATES, NUM_Y_STATES, NUM_Y_STATES).setIdentity();
    dwdwo curr(initial_dxdwo, initial_dzdwo, initial_dydwo);

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
        
        curr.dxdwo = dxdwo_plus;
        curr.dzdwo = dzdwo_plus;
        curr.dydwo = dydwo_plus;
    }
    return curr;
}

std::vector<dwdp> DroneTrajectory::trajSensParam(SimResults const & simResults, int iterations)
{
    std::vector<dwdp> ts;
    ts.reserve(iterations);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> initial_dxdp;
    initial_dxdp.setZero();
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> initial_dzdp;
    initial_dzdp.setZero();
    Eigen::Matrix<double, NUM_Y_STATES, NUM_PARAMETERS> initial_dydp;
    initial_dydp.setZero();
    ts.emplace_back(initial_dxdp, initial_dzdp, initial_dydp);

    Eigen::SparseMatrix<double> I(NUM_PLANT_STATES, NUM_PLANT_STATES);
    I.setIdentity();

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> A;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> B;
    Eigen::PartialPivLU<Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>> solver;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_PARAMETERS> dydp_plus;

    // iterating trajectory sensitivity
    for(int i = 1; i < iterations; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];
        // dwdwo
        dwdp curr = ts[i-1];
        dxdp = curr.dxdp;
        dzdp = curr.dzdp;

        // dxdwo_plus
        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::SparseMatrix<double> dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}); 
        Eigen::SparseMatrix<double> dfdz_curr = dfdz(simResults.stateProgression[i-1]);   
        // Note: dfdp = 0

        A = I - (timestep / 2.0) * dfdx_plus;
        B = dxdp + (timestep / 2.0)*(dfdx_curr*dxdp + (dfdz_plus+dfdz_curr)*dzdp );
        solver.compute(A);
        dxdp_plus.noalias() = solver.solve(B);
        
        // dzdwo for 1 to n
        Eigen::SparseMatrix<double> dhdx_plus = dhdxPlus(simResults.stateProgression[i], simResults.time[i], timestep);
        Eigen::SparseMatrix<double> dhdx_curr = dhdxCurr(simResults.stateProgression[i-1], timestep);
        Eigen::SparseMatrix<double> dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep); 
        Eigen::SparseMatrix<double> dhdz_curr = dhdzCurr(); 
        Eigen::SparseMatrix<double> dhdy_plus = dhdy(timestep);
        Eigen::SparseMatrix<double> dgdz_plus = dgdz(simResults.stateProgression[i]);
        Eigen::SparseMatrix<double> dgdx_plus = dgdx(simResults.stateProgression[i]);
        Eigen::SparseMatrix<double> dgdp_plus = dgdp(simResults.stateProgression[i]);
        Eigen::SparseMatrix<double> dhdp_plus = dhdp(simResults.stateProgression[i], simResults.time[i]);

        dzdp_plus = dhdx_plus * dxdp_plus + dhdx_curr * dxdp + dhdz_curr * dzdp + dhdp_plus;

        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey ++){
            Eigen::Matrix<double, 1, NUM_PARAMETERS> tmp;
            tmp = dhdz_plus.block(zBeforey, 0, 1, zBeforey) * dzdp_plus.topRows(zBeforey);
            dzdp_plus.row(zBeforey) += tmp;
        }  
        
        // dydwo
        dydp_plus = dgdz_plus * dzdp_plus + dgdx_plus * dxdp_plus + dgdp_plus;

        // dzdwo 
        dzdp_plus += dhdy_plus * dydp_plus;
        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::Matrix<double, 1, NUM_PARAMETERS> tmp;
            tmp = dhdz_plus.block(zAftery, 0, 1, zAftery) * dzdp_plus.topRows(zAftery);
            dzdp_plus.row(zAftery) += tmp;
        }
        
        ts.emplace_back(dxdp_plus, dzdp_plus, dydp_plus);
    }

    // std::chrono::time_point end = std::chrono::steady_clock::now();
    // std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // m_logger << "Elapsed Time: " << elapsed.count() << " us" << std::endl;

    return ts;
    
}

std::vector<d2wdwo2> DroneTrajectory::secondOrdertrajSens(SimResults const & simResults, std::vector<dwdwo> const & ts)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    int iterations = simResults.time.size();
    std::vector<d2wdwo2> ts2;
    ts2.reserve(iterations);

    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2xdwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    init_d2xdwo2.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2zdwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    init_d2zdwo2.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2ydwo2(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    init_d2ydwo2.setZero();

    ts2.emplace_back(init_d2xdwo2, init_d2zdwo2, init_d2ydwo2);

    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwo2_curr(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwo2_curr(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwo2_plus(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwo2_plus(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwo2_plus(NUM_Z_STATES, NUM_STATES, NUM_STATES);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_curr;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> I;
    I.setIdentity();

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edxdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edxdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eix_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eixdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVelx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVely_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_ft_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_tx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_ty_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_tz_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> c3_desRoll_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> c3_desPitch_matrix;

    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> dhdy_plus;

    Eigen::Matrix<double, NUM_Y_STATES, NUM_PLANT_STATES> dgdx_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> dgdz_plus;

    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims0 = {Eigen::IndexPair<int>(0, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(1, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims2 = {Eigen::IndexPair<int>(2, 0)};

    for(int i = 1; i < iterations; i++)
    {
        // m_logger << i << std::endl;
        double timestep = simResults.time[i] - simResults.time[i-1];

        d2xdwo2_curr = ts2[i-1].d2xdwo2;
        d2zdwo2_curr = ts2[i-1].d2zdwo2;
        dxdwo_curr = ts[i-1].dxdwo;
        dxdwo_plus = ts[i].dxdwo;
        dzdwo_curr = ts[i-1].dzdwo;
        dzdwo_plus = ts[i].dzdwo;

        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_curr_tensor(dxdwo_curr.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_plus_tensor(dxdwo_plus.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdwo_curr_tensor(dzdwo_curr.data(), NUM_Z_STATES, NUM_STATES);
        
        // d2xdwo2_plus
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_plus = d2fdx2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a1 = d2fdx2_plus.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dxdwo_plus_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_plus = d2fdxdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a2 = d2fdxdz_plus.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dzdwo_curr_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_plus = d2fdzdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a3 = d2fdzdx_plus.contract(dzdwo_curr_tensor, contract_dims).eval().contract(dxdwo_plus_tensor, contract_dims).eval();

        dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_plus_tensor(dfdz_plus.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a4 = dfdz_plus_tensor.contract(d2zdwo2_curr, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_curr = d2fdx2({simResults.stateProgression[i-1].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a5 = d2fdx2_curr.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dxdwo_curr_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_curr = d2fdxdz(simResults.stateProgression[i-1]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a6 = d2fdxdz_curr.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dzdwo_curr_tensor, contract_dims).eval();

        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdx_curr_tensor(dfdx_curr.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a7 = dfdx_curr_tensor.contract(d2xdwo2_curr, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_curr = d2fdzdx(simResults.stateProgression[i-1]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a8 = d2fdzdx_curr.contract(dzdwo_curr_tensor, contract_dims).eval().contract(dxdwo_curr_tensor, contract_dims).eval();

        dfdz_curr = dfdz(simResults.stateProgression[i-1]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_curr_tensor(dfdz_curr.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a9 = dfdz_curr_tensor.contract(d2zdwo2_curr, contract_dims);

        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a10_A_matrix = (I - (timestep/2.0) * dfdx_plus).inverse();
        Eigen::Tensor<double, 3, Eigen::ColMajor> a10_B_tensor = d2xdwo2_curr + timestep/2*(a1+a2+a3+a4+a5+a6+a7+a8+a9);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> a10_A_tensor(a10_A_matrix.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        d2xdwo2_plus = a10_A_tensor.contract(a10_B_tensor, contract_dims).eval();

        // d2zdwo2_plus 
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_curr_tensor = d2hdx2_curr(simResults.stateProgression[i-1], timestep);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b1 = d2hdx2_curr_tensor.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dxdwo_curr_tensor, contract_dims).eval();

        dhdx_curr = dhdxCurr(simResults.stateProgression[i-1], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_curr_tensor(dhdx_curr.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b2 = dhdx_curr_tensor.contract(d2xdwo2_curr, contract_dims).eval();

        // only uses plant state in derivative
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_plus_tensor = d2hdx2_plus(simResults.stateProgression[i], simResults.time[i], timestep);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b3 = d2hdx2_plus_tensor.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dxdwo_plus_tensor, contract_dims).eval();

        dhdx_plus = dhdxPlus(simResults.stateProgression[i], simResults.time[i], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_plus_tensor(dhdx_plus.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b4 = dhdx_plus_tensor.contract(d2xdwo2_plus, contract_dims).eval();

        dhdz_curr = dhdzCurr().toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdz_curr_tensor(dhdz_curr.data(), NUM_Z_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b5 = dhdz_curr_tensor.contract(d2zdwo2_curr, contract_dims).eval();

        d2zdwo2_plus = b1 + b2 + b3 + b4 + b5;

        dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep).toDense();
        Eigen::DSizes<Eigen::Index, 3> offset(0, 0, 0); 
            
        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey++){
            Eigen::MatrixXd blockMatrix = dhdz_plus.block(zBeforey, 0, 1, zBeforey);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrix.data(), zBeforey);
            Eigen::DSizes<Eigen::Index, 3> extent(zBeforey, NUM_STATES, NUM_STATES);
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp = dhdz_plus_block_tensor.contract(d2zdwo2_plus.slice(offset, extent), contract_dims0).eval();
            d2zdwo2_plus.chip(zBeforey, 0) += tmp;        
        }  

        // d2ydwo2_plus
        dgdx_plus = dgdx(simResults.stateProgression[i]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdx_plus_tensor(dgdx_plus.data(), NUM_Y_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c1 = dgdx_plus_tensor.contract(d2xdwo2_plus, contract_dims).eval();

        dgdz_plus = dgdz(simResults.stateProgression[i]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdz_plus_tensor(dgdz_plus.data(), NUM_Y_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c2 = dgdz_plus_tensor.contract(d2zdwo2_plus, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> c3(NUM_Y_STATES, NUM_STATES, NUM_STATES);
        c3.setZero();
        c3_desRoll_matrix = (d2desRolldx2(simResults.stateProgression[i])*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> c3_desRoll_tensor( c3_desRoll_matrix.data(), NUM_STATES, NUM_STATES);
        c3.chip(0, 0) = c3_desRoll_tensor;
        c3_desPitch_matrix = (d2desPitchdx2(simResults.stateProgression[i])*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> c3_desPitch_tensor( c3_desPitch_matrix.data(), NUM_STATES, NUM_STATES);
        c3.chip(1, 0) = c3_desPitch_tensor;

        d2ydwo2_plus = c1 + c2 + c3;

        // d2zdwo2_plus after y
        dhdy_plus = dhdy(timestep);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdy_plus_tensor(dhdy_plus.data(), NUM_Z_STATES, NUM_Y_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b7 = dhdy_plus_tensor.contract(d2ydwo2_plus, contract_dims).eval();
        d2zdwo2_plus+= b7;

        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::MatrixXd blockMatrix = dhdz_plus.block(zAftery, 0, 1, zAftery);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrix.data(), zAftery);
            Eigen::DSizes<Eigen::Index, 3> extent(zAftery, NUM_STATES, NUM_STATES);
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp = dhdz_plus_block_tensor.contract(d2zdwo2_plus.slice(offset, extent), contract_dims0).eval();
            d2zdwo2_plus.chip(zAftery, 0) += tmp;    
        }

        b8_ft_matrix = (d2ft_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b8_ft_tensor( b8_ft_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(ft, 0) += b8_ft_tensor;
        b8_tx_matrix = (d2tx_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b8_tx_tensor( b8_tx_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(tx, 0) += b8_tx_tensor;
        b8_ty_matrix = (d2ty_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b8_ty_tensor( b8_ty_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(ty, 0) += b8_ty_tensor;
        b8_tz_matrix = (d2tz_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b8_tz_tensor( b8_tz_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(tz, 0) += b8_tz_tensor;

        ts2.emplace_back(d2xdwo2_plus, d2zdwo2_plus, d2ydwo2_plus);

    }

    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    m_logger << "Elapsed Time ts2: " << elapsed.count() << " us" << std::endl;

    return ts2;
}

std::vector<d2wdwodp> DroneTrajectory::secondOrdertrajSensParams(SimResults const & simResults, std::vector<dwdwo> const & ts, std::vector<dwdp> const & tsp)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    int iterations = simResults.time.size();
    std::vector<d2wdwodp> ts2;
    ts2.reserve(iterations);

    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2xdwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    init_d2xdwodp.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2zdwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);
    init_d2zdwodp.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> init_d2ydwodp(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);
    init_d2ydwodp.setZero();

    ts2.emplace_back(init_d2xdwodp, init_d2zdwodp, init_d2ydwodp);

    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwodp_curr(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwodp_curr(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwodp_plus(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwodp_plus(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwodp_plus(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;
    
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_curr;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> I;
    I.setIdentity();

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edxdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edxdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eix_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiy_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eixdot_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiydot_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVelx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVely_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b8_ft_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b8_tx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b8_ty_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b8_tz_matrix;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> c3_desRoll_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> c3_desPitch_matrix;

    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> dhdy_plus;

    Eigen::Matrix<double, NUM_Y_STATES, NUM_PLANT_STATES> dgdx_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> dgdz_plus;

    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims0 = {Eigen::IndexPair<int>(0, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(1, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims2 = {Eigen::IndexPair<int>(2, 0)};

    for(int i = 1; i < iterations; i++)
    {
        // m_logger << i << std::endl;
        double timestep = simResults.time[i] - simResults.time[i-1];

        d2xdwodp_curr = ts2[i-1].d2xdwodp;
        d2zdwodp_curr = ts2[i-1].d2zdwodp;
        dxdwo_curr = ts[i-1].dxdwo;
        dxdwo_plus = ts[i].dxdwo;
        dzdwo_curr = ts[i-1].dzdwo;
        dzdwo_plus = ts[i].dzdwo;
        dydwo_plus = ts[i].dydwo;

        dxdp_curr = tsp[i-1].dxdp;
        dxdp_plus = tsp[i].dxdp;
        dzdp_curr = tsp[i-1].dzdp;
        dzdp_plus = tsp[i].dzdp;

        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_curr_tensor(dxdwo_curr.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_plus_tensor(dxdwo_plus.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdwo_curr_tensor(dzdwo_curr.data(), NUM_Z_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdwo_plus_tensor(dzdwo_plus.data(), NUM_Z_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dydwo_plus_tensor(dydwo_plus.data(), NUM_Y_STATES, NUM_STATES);

        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdp_curr_tensor(dxdp_curr.data(), NUM_PLANT_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdp_plus_tensor(dxdp_plus.data(), NUM_PLANT_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdp_curr_tensor(dzdp_curr.data(), NUM_Z_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdp_plus_tensor(dzdp_plus.data(), NUM_Z_STATES, NUM_PARAMETERS);
        
        // d2xdwo2_plus
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_plus = d2fdx2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a1 = d2fdx2_plus.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dxdp_plus_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_plus = d2fdxdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a2 = d2fdxdz_plus.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dzdp_curr_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_plus = d2fdzdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a3 = d2fdzdx_plus.contract(dzdwo_curr_tensor, contract_dims).eval().contract(dxdp_plus_tensor, contract_dims).eval();

        dfdz_plus = dfdz({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_plus_tensor(dfdz_plus.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a4 = dfdz_plus_tensor.contract(d2zdwodp_curr, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_curr = d2fdx2({simResults.stateProgression[i-1].plant, simResults.stateProgression[i-1].alge});
        Eigen::Tensor<double, 3, Eigen::ColMajor> a5 = d2fdx2_curr.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dxdp_curr_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_curr = d2fdxdz(simResults.stateProgression[i-1]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a6 = d2fdxdz_curr.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dzdp_curr_tensor, contract_dims).eval();

        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdx_curr_tensor(dfdx_curr.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a7 = dfdx_curr_tensor.contract(d2xdwodp_curr, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_curr = d2fdzdx(simResults.stateProgression[i-1]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a8 = d2fdzdx_curr.contract(dzdwo_curr_tensor, contract_dims).eval().contract(dxdp_curr_tensor, contract_dims).eval();

        dfdz_curr = dfdz(simResults.stateProgression[i-1]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_curr_tensor(dfdz_curr.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> a9 = dfdz_curr_tensor.contract(d2zdwodp_curr, contract_dims);

        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a10_A_matrix = (I - (timestep/2.0) * dfdx_plus).inverse();
        Eigen::Tensor<double, 3, Eigen::ColMajor> a10_B_tensor = d2xdwodp_curr + timestep/2*(a1+a2+a3+a4+a5+a6+a7+a8+a9);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> a10_A_tensor(a10_A_matrix.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        d2xdwodp_plus = a10_A_tensor.contract(a10_B_tensor, contract_dims).eval();

        // d2zdwodp_plus 
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_curr_tensor = d2hdx2_curr(simResults.stateProgression[i-1], timestep);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b1 = d2hdx2_curr_tensor.contract(dxdwo_curr_tensor, contract_dims).eval().contract(dxdp_curr_tensor, contract_dims).eval();

        dhdx_curr = dhdxCurr(simResults.stateProgression[i-1], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_curr_tensor(dhdx_curr.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b2 = dhdx_curr_tensor.contract(d2xdwodp_curr, contract_dims).eval();

        // only uses plant state in derivative
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_plus_tensor = d2hdx2_plus(simResults.stateProgression[i], simResults.time[i], timestep);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b3 = d2hdx2_plus_tensor.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dxdp_plus_tensor, contract_dims).eval();

        dhdx_plus = dhdxPlus(simResults.stateProgression[i], simResults.time[i], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_plus_tensor(dhdx_plus.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b4 = dhdx_plus_tensor.contract(d2xdwodp_plus, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdxplusdp = d2hdxplus_dp(simResults.stateProgression[i], simResults.time[i]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b8 = d2hdxplusdp.contract(dxdwo_plus_tensor, contract_dims);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b9 = b8.shuffle(Eigen::array<int,3>{0,2,1});

        dhdz_curr = dhdzCurr().toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdz_curr_tensor(dhdz_curr.data(), NUM_Z_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b5 = dhdz_curr_tensor.contract(d2zdwodp_curr, contract_dims).eval();

        d2zdwodp_plus = b1 + b2 + b3 + b4 + b5 + b9;

        dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep).toDense();
        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdzplusdp = d2hdzplus_dp();
        Eigen::DSizes<Eigen::Index, 3> offset3(0, 0, 0); 
        Eigen::DSizes<Eigen::Index, 2> offset2(0, 0); 
            
        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey++){
            Eigen::MatrixXd blockMatrixdhdz = dhdz_plus.block(zBeforey, 0, 1, zBeforey);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrixdhdz.data(), zBeforey);
            Eigen::DSizes<Eigen::Index, 3> extent(zBeforey, NUM_STATES, NUM_PARAMETERS);
            Eigen::DSizes<Eigen::Index, 2> extent2(zBeforey, NUM_STATES);
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp1 = dhdz_plus_block_tensor.contract(d2zdwodp_plus.slice(offset3, extent), contract_dims0).eval();

            Eigen::DSizes<Eigen::Index, 3> offsetblock(zBeforey, 0, 0);
            Eigen::DSizes<Eigen::Index, 3> extentblock(1, zBeforey, NUM_PARAMETERS);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp2 = d2hdzplusdp.slice(offsetblock, extentblock).contract(dzdwo_plus_tensor.slice(offset2, extent2), contract_dims);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp3 = tmp2.shuffle(Eigen::array<int,3>{0,2,1});
        
            d2zdwodp_plus.chip(zBeforey, 0) += tmp1 + tmp3.chip(0, 0);        
        }  

        // d2ydwodp_plus
        dgdx_plus = dgdx(simResults.stateProgression[i]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdx_plus_tensor(dgdx_plus.data(), NUM_Y_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c1 = dgdx_plus_tensor.contract(d2xdwodp_plus, contract_dims).eval();

        dgdz_plus = dgdz(simResults.stateProgression[i]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdz_plus_tensor(dgdz_plus.data(), NUM_Y_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c2 = dgdz_plus_tensor.contract(d2zdwodp_plus, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdx2_plus = d2gdx2(simResults.stateProgression[i]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c3 = d2gdx2_plus.contract(dxdwo_plus_tensor, contract_dims).eval().contract(dxdp_plus_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdxplusdp = d2gdxplus_dp(simResults.stateProgression[i]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c4 = d2gdxplusdp.contract(dxdwo_plus_tensor, contract_dims);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c6 = c4.shuffle(Eigen::array<int,3>{0,2,1});

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdzplusdp = d2gdzplus_dp(simResults.stateProgression[i]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c5 = d2gdzplusdp.contract(dzdwo_plus_tensor, contract_dims);
        Eigen::Tensor<double, 3, Eigen::ColMajor> c7 = c5.shuffle(Eigen::array<int,3>{0,2,1});

        d2ydwodp_plus = c1 + c2 + c3 + c6 + c7;

        // d2zdwodp_plus after y
        dhdy_plus = dhdy(timestep);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdy_plus_tensor(dhdy_plus.data(), NUM_Z_STATES, NUM_Y_STATES);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b7 = dhdy_plus_tensor.contract(d2ydwodp_plus, contract_dims).eval();

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdydp_plus = d2hdydp();
        Eigen::Tensor<double, 3, Eigen::ColMajor> b10 = d2hdydp_plus.contract(dydwo_plus_tensor, contract_dims).eval();
        Eigen::Tensor<double, 3, Eigen::ColMajor> b11 = b10.shuffle(Eigen::array<int,3>{0,2,1});

        d2zdwodp_plus+= b7 + b11;

        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::MatrixXd blockMatrix = dhdz_plus.block(zAftery, 0, 1, zAftery);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrix.data(), zAftery);
            Eigen::DSizes<Eigen::Index, 3> extent(zAftery, NUM_STATES, NUM_PARAMETERS);
            Eigen::DSizes<Eigen::Index, 2> extent2(zAftery, NUM_STATES);
            
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp1 = dhdz_plus_block_tensor.contract(d2zdwodp_plus.slice(offset3, extent), contract_dims0).eval();

            Eigen::DSizes<Eigen::Index, 3> offsetblock(zAftery, 0, 0);
            Eigen::DSizes<Eigen::Index, 3> extentblock(1, zAftery, NUM_PARAMETERS);

            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp2 = d2hdzplusdp.slice(offsetblock, extentblock).contract(dzdwo_plus_tensor.slice(offset2, extent2), contract_dims);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp3 = tmp2.shuffle(Eigen::array<int,3>{0,2,1});
            d2zdwodp_plus.chip(zAftery, 0) += tmp1 + tmp3.chip(0, 0);    
        }

        Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdz_plus2 = d2hdzplus2(simResults.stateProgression[i]);
        Eigen::Tensor<double, 3, Eigen::ColMajor> b12 = d2hdz_plus2.contract(dzdwo_plus_tensor, contract_dims).eval().contract(dzdp_plus_tensor, contract_dims).eval();
        d2zdwodp_plus += b12;
       
        ts2.emplace_back(d2xdwodp_plus, d2zdwodp_plus, d2ydwodp_plus);

    }

    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    m_logger << "Elapsed Time ts2p: " << elapsed.count() << " us" << std::endl;

    return ts2;
}


G_tp DroneTrajectory::calc_G_tp(std::vector<dwdwo> const& trajSens)
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
