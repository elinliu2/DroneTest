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

std::vector<dwdp> DroneTrajectory::trajSensParam(SimResults const & simResults, G_tp gtp)
{
    const int iterations = gtp.tp + 1;
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

    const int iterations = simResults.time.size();
    m_logger << "iterations: " << iterations << std::endl;
    std::vector<d2wdwo2> ts2;
    ts2.reserve(iterations);

    Eigen::Tensor<double, 3> init_d2xdwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    init_d2xdwo2.setZero();
    Eigen::Tensor<double, 3> init_d2zdwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    init_d2zdwo2.setZero();
    Eigen::Tensor<double, 3> init_d2ydwo2(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    init_d2ydwo2.setZero();

    ts2.emplace_back(init_d2xdwo2, init_d2zdwo2, init_d2ydwo2);

    Eigen::Tensor<double, 3> d2xdwo2_curr(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3> d2zdwo2_curr(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3> d2xdwo2_plus(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3> d2ydwo2_plus(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3> d2zdwo2_plus(NUM_Z_STATES, NUM_STATES, NUM_STATES);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;

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

    for(int i = 1; i < iterations; i++)
    {
        m_logger << i << std::endl;
        double timestep = simResults.time[i] - simResults.time[i-1];

        d2xdwo2_curr = ts2[i-1].d2xdwo2;
        d2zdwo2_curr = ts2[i-1].d2zdwo2;
        dxdwo_curr = ts[i-1].dxdwo;
        dxdwo_plus = ts[i].dxdwo;
        dzdwo_plus = ts[i].dzdwo;

        // d2xdwo2_plus
        Eigen::SparseMatrix<double> d2fdx2_plus = d2fdx2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        Eigen::SparseMatrix<double> d2fdx2_curr = d2fdx2({simResults.stateProgression[i-1].plant, simResults.stateProgression[i-1].alge});
        dfdx_plus = dfdx({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge});
        dfdx_curr = dfdx(simResults.stateProgression[i-1]);
        
        Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_STATES> a1_A = d2fdx2_plus*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 3>> a1_A_tensor(a1_A.data(), NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a1_B_tensor(dxdwo_plus.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(1, 0)};
        Eigen::Tensor<double, 3> a1 = a1_A_tensor.contract(a1_B_tensor, contract_dims).eval();

        Eigen::Tensor<double, 3> a2 = dfdz_mult_d2zdwo2({simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge}, d2zdwo2_curr);

        Eigen::Matrix<double, NUM_PLANT_STATES*NUM_PLANT_STATES, NUM_STATES> a3_A = d2fdx2_curr*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 3>> a3_A_tensor(a3_A.data(), NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a3_B_tensor(dxdwo_curr.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::Tensor<double, 3> a3 = a3_A_tensor.contract(a3_B_tensor, contract_dims).eval();

        Eigen::TensorMap<Eigen::Tensor<double, 2>> a4_A_tensor(dfdx_curr.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3> a4 = a4_A_tensor.contract(d2xdwo2_curr, contract_dims).eval();

        Eigen::Tensor<double, 3> a5 = dfdz_mult_d2zdwo2(simResults.stateProgression[i-1], d2zdwo2_curr);

        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a6_A_matrix = (I - dfdx_plus).inverse();
        Eigen::Tensor<double, 3> a6_B_tensor = d2xdwo2_curr + timestep/2*(a1+a2+a3+a4+a5);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> a6_A_tensor(a6_A_matrix.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        d2xdwo2_plus = a6_A_tensor.contract(a6_B_tensor, contract_dims).eval();

        // d2zdwo2_plus 
        Eigen::Tensor<double, 3> b1(NUM_Z_STATES, NUM_STATES, NUM_STATES);
        b1.setZero();
        b1_edx_matrix = (d2edx_dx_curr2(simResults.stateProgression[i-1], timestep)*dxdwo_curr).transpose()*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b1_edx_tensor( b1_edx_matrix.data(), NUM_STATES, NUM_STATES);
        b1.chip(edx, 0) = b1_edx_tensor;
        b1_edy_matrix = (d2edy_dx_curr2(simResults.stateProgression[i-1], timestep)*dxdwo_curr).transpose()*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b1_edy_tensor( b1_edy_matrix.data(), NUM_STATES, NUM_STATES);
        b1.chip(edy, 0) = b1_edy_tensor;
        b1_edxdot_matrix = (d2edxdot_dx_curr2(simResults.stateProgression[i-1], timestep)*dxdwo_curr).transpose()*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b1_edxdot_tensor( b1_edxdot_matrix.data(), NUM_STATES, NUM_STATES);
        b1.chip(edxdot, 0) = b1_edxdot_tensor;
        b1_edydot_matrix = (d2edydot_dx_curr2(simResults.stateProgression[i-1], timestep)*dxdwo_curr).transpose()*dxdwo_curr;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b1_edydot_tensor( b1_edydot_matrix.data(), NUM_STATES, NUM_STATES);
        b1.chip(edydot, 0) = b1_edydot_tensor;

        dhdx_curr = dhdxCurr(simResults.stateProgression[i-1], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b2_A_tensor(dhdx_curr.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3> b2 = b2_A_tensor.contract(d2xdwo2_curr, contract_dims).eval();

        Eigen::Tensor<double, 3> b3(NUM_Z_STATES, NUM_STATES, NUM_STATES);
        b3.setZero();
        b3_eix_matrix = (d2eix_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_eix_tensor( b3_eix_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(eix, 0) = b3_eix_tensor;
        b3_eiy_matrix = (d2eiy_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_eiy_tensor( b3_eiy_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(eiy, 0) = b3_eiy_tensor;
        b3_eixdot_matrix = (d2eixdot_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_eixdot_tensor( b3_eixdot_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(eixdot, 0) = b3_eixdot_tensor;
        b3_eiydot_matrix = (d2eiydot_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_eiydot_tensor( b3_eiydot_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(eiydot, 0) = b3_eiydot_tensor;

        b3_edx_matrix = (d2edx_dx_plus2(simResults.stateProgression[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_edx_tensor( b3_edx_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(edx, 0) = b3_edx_tensor;
        b3_edy_matrix = (d2edy_dx_plus2(simResults.stateProgression[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_edy_tensor( b3_edy_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(edy, 0) = b3_edy_tensor;
        b3_edxdot_matrix = (d2edxdot_dx_plus2(simResults.stateProgression[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_edxdot_tensor( b3_edxdot_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(edxdot, 0) = b3_edxdot_tensor;
        b3_edydot_matrix = (d2edydot_dx_plus2(simResults.stateProgression[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_edydot_tensor( b3_edydot_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(edydot, 0) = b3_edydot_tensor;

        b3_desVelx_matrix = (d2desVelx_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_desVelx_tensor( b3_desVelx_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(desVelX, 0) = b3_desVelx_tensor;
        b3_desVely_matrix = (d2desVely_dx_plus2(simResults.stateProgression[i], simResults.time[i], timestep)*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b3_desVely_tensor( b3_desVely_matrix.data(), NUM_STATES, NUM_STATES);
        b3.chip(desVelY, 0) = b3_desVely_tensor;

        dhdx_plus = dhdxPlus(simResults.stateProgression[i], simResults.time[i], timestep).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b4_A_tensor(dhdx_plus.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3> b4 = b4_A_tensor.contract(d2xdwo2_plus, contract_dims).eval();

        dhdz_curr = dhdzCurr().toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b5_A_tensor(dhdz_curr.data(), NUM_Z_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3> b5 = b5_A_tensor.contract(d2zdwo2_curr, contract_dims).eval();

        d2zdwo2_plus = b1 + b2 + b3 + b4 + b5;

        Eigen::Tensor<double, 3> b6(NUM_Z_STATES, NUM_STATES, NUM_STATES);
        dhdz_plus = dhdzPlus(simResults.stateProgression[i], timestep).toDense();
        Eigen::DSizes<Eigen::Index, 3> offset(0, 0, 0); // Start at (0, 2, 0)
            
        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey++){
            Eigen::TensorMap<Eigen::Tensor<double, 2>> b6_A_tensor(dhdz_plus.block(zBeforey, 0, 1, zBeforey).data(), 1, zBeforey);
            Eigen::DSizes<Eigen::Index, 3> extent(zBeforey, NUM_STATES, NUM_STATES);

            // Contraction yields [1, NUM_STATES, NUM_STATES] — rank 3
            Eigen::Tensor<double, 3> tmp = b6_A_tensor
                .contract(d2zdwo2_plus.slice(offset, extent), contract_dims)
                .eval();

            // Squeeze the leading 1-dimension before adding
            Eigen::array<Eigen::Index, 2> shape{NUM_STATES, NUM_STATES};
            d2zdwo2_plus.chip(zBeforey, 0) += tmp.reshape(shape);
                    
            }  

        // d2ydwo2_plus
        dgdx_plus = dgdx(simResults.stateProgression[i]).toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2>> c1_A_tensor(dgdx_plus.data(), NUM_Y_STATES, NUM_PLANT_STATES);
        Eigen::Tensor<double, 3> c1 = c1_A_tensor.contract(d2xdwo2_plus, contract_dims).eval();

        dgdz_plus = dgdz(simResults.stateProgression[i]);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> c2_A_tensor(dgdz_plus.data(), NUM_Y_STATES, NUM_Z_STATES);
        Eigen::Tensor<double, 3> c2 = c2_A_tensor.contract(d2zdwo2_plus, contract_dims).eval();

        Eigen::Tensor<double, 3> c3(NUM_Y_STATES, NUM_STATES, NUM_STATES);
        c3.setZero();
        c3_desRoll_matrix = (d2desRolldx2(simResults.stateProgression[i])*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> c3_desRoll_tensor( c3_desRoll_matrix.data(), NUM_STATES, NUM_STATES);
        c3.chip(0, 0) = c3_desRoll_tensor;
        c3_desPitch_matrix = (d2desPitchdx2(simResults.stateProgression[i])*dxdwo_plus).transpose()*dxdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> c3_desPitch_tensor( c3_desPitch_matrix.data(), NUM_STATES, NUM_STATES);
        c3.chip(1, 0) = c3_desPitch_tensor;

        d2ydwo2_plus = c1 + c2 + c3;

        // d2zdwo2_plus after y
        dhdy_plus = dhdy(timestep);
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b7_A_tensor(dhdy_plus.data(), NUM_Z_STATES, NUM_Y_STATES);
        Eigen::Tensor<double, 3> b7 = b7_A_tensor.contract(d2ydwo2_plus, contract_dims).eval();
        d2zdwo2_plus+= b7;

        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::TensorMap<Eigen::Tensor<double, 2>> b6_A_tensor(dhdz_plus.block(zAftery, 0, 1, zAftery).data(), 1, zAftery);
            Eigen::DSizes<Eigen::Index, 3> extent(zAftery, NUM_STATES, NUM_STATES); 

                Eigen::Tensor<double, 3> tmp = b6_A_tensor
            .contract(d2zdwo2_plus.slice(offset, extent), contract_dims)
            .eval();

            // Squeeze the leading 1-dimension before adding
            Eigen::array<Eigen::Index, 2> shape{NUM_STATES, NUM_STATES};
            d2zdwo2_plus.chip(zAftery, 0) += tmp.reshape(shape);
        
        }

        b8_ft_matrix = (d2ft_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b8_ft_tensor( b8_ft_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(ft, 0) += b8_ft_tensor;
        b8_tx_matrix = (d2tx_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b8_tx_tensor( b8_tx_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(tx, 0) += b8_tx_tensor;
        b8_ty_matrix = (d2ty_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b8_ty_tensor( b8_ty_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(ty, 0) += b8_ty_tensor;
        b8_tz_matrix = (d2tz_dz_plus2(simResults.stateProgression[i])*dzdwo_plus).transpose()*dzdwo_plus;
        Eigen::TensorMap<Eigen::Tensor<double, 2>> b8_tz_tensor( b8_tz_matrix.data(), NUM_STATES, NUM_STATES);
        d2zdwo2_plus.chip(tz, 0) += b8_tz_tensor;

        ts2.emplace_back(d2xdwo2_plus, d2zdwo2_plus, d2ydwo2_plus);

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
