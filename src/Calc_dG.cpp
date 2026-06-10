#include "DroneTrajectory.h"

d2w DroneTrajectory::calc_d2w(SimResults const & simResults, std::vector<dwdwo> const & ts, G_tp gtp)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp_curr; dxdp_curr.setZero();
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp_curr; dzdp_curr.setZero();
    Eigen::Matrix<double, NUM_Y_STATES, NUM_PARAMETERS> dydp_curr; dydp_curr.setZero();

    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwo2_curr(NUM_PLANT_STATES, NUM_STATES, NUM_STATES); d2xdwo2_curr.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwo2_curr(NUM_Z_STATES, NUM_STATES, NUM_STATES); d2zdwo2_curr.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwo2_curr(NUM_Y_STATES, NUM_STATES, NUM_STATES); d2ydwo2_curr.setZero();

    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwodp_curr(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS); d2xdwodp_curr.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwodp_curr(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS); d2zdwodp_curr.setZero();
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwodp_curr(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS); d2ydwodp_curr.setZero();

    // df
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx_plus;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_Z_STATES> dfdz_plus;

    // dh
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PLANT_STATES> dhdx_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Z_STATES> dhdz_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_Y_STATES> dhdy_plus;

    // dg
    Eigen::Matrix<double, NUM_Y_STATES, NUM_PLANT_STATES> dgdx_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_Z_STATES> dgdz_plus;

    // d2hdzplus2 terms
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_ft_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_tx_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_ty_matrix;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b8_tz_matrix;

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> I; I.setIdentity();
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> A_tsp;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> B_tsp;
    Eigen::PartialPivLU<Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES>> solver_tsp;

    // ts
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_curr;
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_curr;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo_plus;

    // tsp
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PARAMETERS> dxdp_plus;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_PARAMETERS> dzdp_plus;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_PARAMETERS> dydp_plus;

    // ts2
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwo2_plus(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwo2_plus(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwo2_plus(NUM_Z_STATES, NUM_STATES, NUM_STATES);

    // ts2p
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2xdwodp_plus(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2ydwodp_plus(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3,Eigen::ColMajor> d2zdwodp_plus(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims0 = {Eigen::IndexPair<int>(0, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims = {Eigen::IndexPair<int>(1, 0)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_dims2 = {Eigen::IndexPair<int>(2, 0)};

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_plus(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a1(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> a1_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a1_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_plus(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a2(NUM_PLANT_STATES, NUM_Z_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> a2_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> a2_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS); 

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_plus(NUM_PLANT_STATES, NUM_Z_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a3(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a3_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a3_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> a4_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a4_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdx2_curr(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES) ;
    Eigen::Tensor<double, 3, Eigen::ColMajor> a5(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a5_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a5_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdxdz_curr(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a6(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> a6_dwo2(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> a6_dwodp(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_Z_STATES); 

    Eigen::Tensor<double, 3, Eigen::ColMajor> a7_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a7_dwodp (NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2fdzdx_curr (NUM_PLANT_STATES, NUM_Z_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a8(NUM_PLANT_STATES, NUM_PLANT_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a8_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a8_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> a9_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a9_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);

    // note: no d2fdz2 d2fdxdp d2fdzdp
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a10_A_matrix;
    Eigen::Tensor<double, 3, Eigen::ColMajor> a10_B_tensor_dwo2(NUM_PLANT_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> a10_B_tensor_dwodp(NUM_PLANT_STATES, NUM_STATES, NUM_PARAMETERS);
    
    // d2zdwo2_plus // d2zdwodp_plus states 1 to n
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_curr_tensor(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b1(NUM_Z_STATES, NUM_PLANT_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b1_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b1_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS); 

    Eigen::Tensor<double, 3, Eigen::ColMajor> b2_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);  
    Eigen::Tensor<double, 3, Eigen::ColMajor> b2_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    // only uses plant state in derivative
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdx2_plus_tensor(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b3(NUM_Z_STATES, NUM_PLANT_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b3_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b3_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> b4_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b4_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> b5_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b5_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdxplusdp(NUM_Z_STATES, NUM_PLANT_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b6(NUM_Z_STATES, NUM_PARAMETERS, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b7(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdzplusdp(NUM_Z_STATES, NUM_Z_STATES, NUM_PARAMETERS);
    
    Eigen::Tensor<double, 3, Eigen::ColMajor> c1_dwo2(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c1_dwodp(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> c2_dwo2(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c2_dwodp(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdx2_plus(NUM_Y_STATES, NUM_PLANT_STATES, NUM_PLANT_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c3(NUM_Y_STATES, NUM_PLANT_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c3_dwo2(NUM_Y_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c3_dwodp(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdxplusdp(NUM_Y_STATES, NUM_PLANT_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c4(NUM_Y_STATES, NUM_PARAMETERS, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c5(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3, Eigen::ColMajor> d2gdzplusdp (NUM_Y_STATES, NUM_Z_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c6(NUM_Y_STATES, NUM_PARAMETERS, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> c7(NUM_Y_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> b8_dwo2(NUM_Z_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b8_dwodp(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdydp_plus(NUM_Z_STATES, NUM_Y_STATES, NUM_PARAMETERS);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b9(NUM_Z_STATES, NUM_PARAMETERS, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b10(NUM_Z_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Tensor<double, 3, Eigen::ColMajor> d2hdz_plus2(NUM_STATES, NUM_STATES, NUM_STATES);
    Eigen::Tensor<double, 3, Eigen::ColMajor> b12(NUM_STATES, NUM_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b12_dwo2(NUM_STATES, NUM_STATES, NUM_STATES); 
    Eigen::Tensor<double, 3, Eigen::ColMajor> b12_dwodp(NUM_STATES, NUM_STATES, NUM_PARAMETERS);

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edx_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edy_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edxdot_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b1_edydot_matrix_dwo2;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b1_edx_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b1_edy_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b1_edxdot_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b1_edydot_matrix_dwodp;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edx_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edy_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edxdot_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_edydot_matrix_dwo2;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eix_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiy_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eixdot_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_eiydot_matrix_dwo2;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVelx_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b3_desVely_matrix_dwo2;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_edx_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_edy_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_edxdot_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_edydot_matrix_dwodp;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_eix_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_eiy_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_eixdot_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_eiydot_matrix_dwodp;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_desVelx_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b3_desVely_matrix_dwodp;

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b12_ft_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b12_tx_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b12_ty_matrix_dwo2;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> b12_tz_matrix_dwo2;

    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b12_ft_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b12_tx_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b12_ty_matrix_dwodp;
    Eigen::Matrix<double, NUM_STATES, NUM_PARAMETERS> b12_tz_matrix_dwodp;

    for(int i = 1; i <= gtp.tp; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];
        double time = simResults.time[i];
        SystemState state_plus = simResults.stateProgression[i];
        SystemState state_xplus_zcurr = {simResults.stateProgression[i].plant, simResults.stateProgression[i-1].alge};
        SystemState state_curr = simResults.stateProgression[i-1];

        // calc dxdp
        dfdx_plus = dfdx(state_xplus_zcurr);
        dfdx_curr = dfdx(state_curr);
        Eigen::SparseMatrix<double> dfdz_plus_sparse = dfdz(state_xplus_zcurr); 
        Eigen::SparseMatrix<double> dfdz_curr_sparse = dfdz(state_curr);   
        // Note: dfdp = 0
        A_tsp = I - (timestep / 2.0) * dfdx_plus;
        B_tsp = dxdp_curr + (timestep / 2.0)*(dfdx_curr*dxdp_curr + (dfdz_plus_sparse+dfdz_curr_sparse)*dzdp_curr);
        solver_tsp.compute(A_tsp);
        dxdp_plus.noalias() = solver_tsp.solve(B_tsp);

        // dzdp for 1 to n
        Eigen::SparseMatrix<double> dhdx_plus_sparse = dhdxPlus(state_plus, time, timestep);
        Eigen::SparseMatrix<double> dhdx_curr_sparse = dhdxCurr(state_curr, timestep);
        Eigen::SparseMatrix<double> dhdz_plus_sparse = dhdzPlus(state_plus, timestep); 
        Eigen::SparseMatrix<double> dhdz_curr_sparse = dhdzCurr(); 
        Eigen::SparseMatrix<double> dhdp_plus_sparse = dhdp(state_plus, time);
        dzdp_plus = dhdx_plus_sparse * dxdp_plus + dhdx_curr_sparse * dxdp_curr + dhdz_curr_sparse * dzdp_curr + dhdp_plus_sparse;

        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey ++){
            Eigen::Matrix<double, 1, NUM_PARAMETERS> tmp;
            tmp = dhdz_plus_sparse.block(zBeforey, 0, 1, zBeforey) * dzdp_plus.topRows(zBeforey);
            dzdp_plus.row(zBeforey) += tmp;
        }  
        
        // dydp
        Eigen::SparseMatrix<double> dgdz_plus_sparse = dgdz(state_plus);
        Eigen::SparseMatrix<double> dgdx_plus_sparse = dgdx(state_plus);
        Eigen::SparseMatrix<double> dgdp_plus_sparse = dgdp(state_plus);
        dydp_plus = dgdz_plus_sparse * dzdp_plus + dgdx_plus_sparse * dxdp_plus + dgdp_plus_sparse;

        // dzdp n to m
        Eigen::SparseMatrix<double> dhdy_plus_sparse = dhdy(timestep);
        dzdp_plus += dhdy_plus_sparse * dydp_plus;
        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::Matrix<double, 1, NUM_PARAMETERS> tmp;
            tmp = dhdz_plus_sparse.block(zAftery, 0, 1, zAftery) * dzdp_plus.topRows(zAftery);
            dzdp_plus.row(zAftery) += tmp;
        }
        
        // second order sensitivities
        dxdwo_curr = ts[i-1].dxdwo;
        dxdwo_plus = ts[i].dxdwo;
        dzdwo_curr = ts[i-1].dzdwo;
        dzdwo_plus = ts[i].dzdwo;
        dydwo_plus = ts[i].dydwo;
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_curr_tensor(dxdwo_curr.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdwo_plus_tensor(dxdwo_plus.data(), NUM_PLANT_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdwo_curr_tensor(dzdwo_curr.data(), NUM_Z_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdwo_plus_tensor(dzdwo_plus.data(), NUM_Z_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dydwo_plus_tensor(dydwo_plus.data(), NUM_Y_STATES, NUM_STATES);

        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdp_curr_tensor(dxdp_curr.data(), NUM_PLANT_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dxdp_plus_tensor(dxdp_plus.data(), NUM_PLANT_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdp_curr_tensor(dzdp_curr.data(), NUM_Z_STATES, NUM_PARAMETERS);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dzdp_plus_tensor(dzdp_plus.data(), NUM_Z_STATES, NUM_PARAMETERS);

        // d2xdwo2_plus // d2xdwodp
        d2fdx2_plus = d2fdx2(state_xplus_zcurr);
        a1 = d2fdx2_plus.contract(dxdwo_plus_tensor, contract_dims).eval();
        a1_dwo2  = a1.contract(dxdwo_plus_tensor, contract_dims).eval();
        a1_dwodp = a1.contract(dxdp_plus_tensor, contract_dims).eval();
        
        d2fdxdz_plus = d2fdxdz(state_xplus_zcurr);
        a2 = d2fdxdz_plus.contract(dxdwo_plus_tensor, contract_dims).eval();
        a2_dwo2  = a2.contract(dzdwo_curr_tensor, contract_dims).eval();
        a2_dwodp = a2.contract(dzdp_curr_tensor, contract_dims).eval();

        d2fdzdx_plus = d2fdzdx(state_xplus_zcurr);
        a3 = d2fdzdx_plus.contract(dzdwo_curr_tensor, contract_dims).eval();
        a3_dwo2  = a3.contract(dxdwo_plus_tensor, contract_dims).eval();
        a3_dwodp = a3.contract(dxdp_plus_tensor, contract_dims).eval();

        dfdz_plus = dfdz_plus_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_plus_tensor(dfdz_plus.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        a4_dwo2  = dfdz_plus_tensor.contract(d2zdwo2_curr, contract_dims).eval();
        a4_dwodp = dfdz_plus_tensor.contract(d2zdwodp_curr, contract_dims).eval();

        d2fdx2_curr = d2fdx2(state_curr);
        a5 = d2fdx2_curr.contract(dxdwo_curr_tensor, contract_dims).eval();
        a5_dwo2  = a5.contract(dxdwo_curr_tensor, contract_dims).eval();
        a5_dwodp = a5.contract(dxdp_curr_tensor, contract_dims).eval();

        d2fdxdz_curr = d2fdxdz(state_curr);
        a6 = d2fdxdz_curr.contract(dxdwo_curr_tensor, contract_dims).eval();
        a6_dwo2  = a6.contract(dzdwo_curr_tensor, contract_dims).eval();
        a6_dwodp = a6.contract(dzdp_curr_tensor, contract_dims).eval();

        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdx_curr_tensor(dfdx_curr.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        a7_dwo2  = dfdx_curr_tensor.contract(d2xdwo2_curr, contract_dims).eval();
        a7_dwodp = dfdx_curr_tensor.contract(d2xdwodp_curr, contract_dims).eval();

        d2fdzdx_curr = d2fdzdx(state_curr);
        a8 = d2fdzdx_curr.contract(dzdwo_curr_tensor, contract_dims).eval();
        a8_dwo2  = a8.contract(dxdwo_curr_tensor, contract_dims).eval();
        a8_dwodp = a8.contract(dxdp_curr_tensor, contract_dims).eval();

        dfdz_curr = dfdz_curr_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> dfdz_curr_tensor(dfdz_curr.data(), NUM_PLANT_STATES, NUM_Z_STATES);
        a9_dwo2  = dfdz_curr_tensor.contract(d2zdwo2_curr, contract_dims);
        a9_dwodp = dfdz_curr_tensor.contract(d2zdwodp_curr, contract_dims);

        // note: no d2fdz2 d2fdxdp d2fdzdp
        Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> a10_A_matrix = (I - (timestep/2.0) * dfdx_plus).inverse();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> a10_A_tensor(a10_A_matrix.data(), NUM_PLANT_STATES, NUM_PLANT_STATES);
        a10_B_tensor_dwo2 = d2xdwo2_curr + timestep/2*(a1_dwo2+a2_dwo2+a3_dwo2+a4_dwo2+a5_dwo2+a6_dwo2+a7_dwo2+a8_dwo2+a9_dwo2);
        a10_B_tensor_dwodp = d2xdwodp_curr + timestep/2*(a1_dwodp+a2_dwodp+a3_dwodp+a4_dwodp+a5_dwodp+a6_dwodp+a7_dwodp+a8_dwodp+a9_dwodp);
        d2xdwo2_plus = a10_A_tensor.contract(a10_B_tensor_dwo2, contract_dims).eval();
        d2xdwodp_plus = a10_A_tensor.contract(a10_B_tensor_dwodp, contract_dims).eval();

        // d2zdwo2_plus // d2zdwodp_plus states 1 to n
        dhdx_curr = dhdx_curr_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_curr_tensor(dhdx_curr.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        b2_dwo2  = dhdx_curr_tensor.contract(d2xdwo2_curr, contract_dims).eval();
        b2_dwodp = dhdx_curr_tensor.contract(d2xdwodp_curr, contract_dims).eval();

        dhdx_plus = dhdx_plus_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdx_plus_tensor(dhdx_plus.data(), NUM_Z_STATES, NUM_PLANT_STATES);
        b4_dwo2  = dhdx_plus_tensor.contract(d2xdwo2_plus, contract_dims).eval();
        b4_dwodp = dhdx_plus_tensor.contract(d2xdwodp_plus, contract_dims).eval();

        dhdz_curr = dhdz_curr_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdz_curr_tensor(dhdz_curr.data(), NUM_Z_STATES, NUM_Z_STATES);
        b5_dwo2  = dhdz_curr_tensor.contract(d2zdwo2_curr, contract_dims).eval();
        b5_dwodp = dhdz_curr_tensor.contract(d2zdwodp_curr, contract_dims).eval();

        d2hdxplusdp = d2hdxplus_dp(state_plus, time);
        b6 = d2hdxplusdp.contract(dxdwo_plus_tensor, contract_dims);
        b7 = b6.shuffle(Eigen::array<int,3>{0,2,1});

        d2zdwo2_plus = b2_dwo2 +  b4_dwo2 + b5_dwo2;
        d2zdwodp_plus = b2_dwodp + b4_dwodp + b5_dwodp + b7;

        // d2hdx2_curr_tensor = d2hdx2_curr(state_curr, timestep);
        // b1 = d2hdx2_curr_tensor.contract(dxdwo_curr_tensor, contract_dims).eval();
        // b1_dwo2  = b1.contract(dxdwo_curr_tensor, contract_dims).eval();
        // b1_dwodp = b1.contract(dxdp_curr_tensor, contract_dims).eval();
        // d2zdwo2_plus += b1_dwo2;
        // d2zdwodp_plus += b1_dwodp;

        b1_edx_matrix_dwo2  = ((d2edx_dx_curr2(state_curr, timestep)*dxdwo_curr).transpose()*dxdwo_curr).transpose();
        b1_edx_matrix_dwodp = ((d2edx_dx_curr2(state_curr, timestep)*dxdp_curr).transpose()*dxdwo_curr).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edx_tensor_dwo2( b1_edx_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edx_tensor_dwodp( b1_edx_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edx, 0) += b1_edx_tensor_dwo2;
        d2zdwodp_plus.chip(edx, 0) += b1_edx_tensor_dwodp;

        b1_edy_matrix_dwo2  = ((d2edy_dx_curr2(state_curr, timestep)*dxdwo_curr).transpose()*dxdwo_curr).transpose();
        b1_edy_matrix_dwodp = ((d2edy_dx_curr2(state_curr, timestep)*dxdp_curr).transpose()*dxdwo_curr).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edy_tensor_dwo2( b1_edy_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edy_tensor_dwodp( b1_edy_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edy, 0) += b1_edy_tensor_dwo2;
        d2zdwodp_plus.chip(edy, 0) += b1_edy_tensor_dwodp;

        b1_edxdot_matrix_dwo2  = ((d2edxdot_dx_curr2(state_curr, timestep)*dxdwo_curr).transpose()*dxdwo_curr).transpose();
        b1_edxdot_matrix_dwodp = ((d2edxdot_dx_curr2(state_curr, timestep)*dxdp_curr).transpose()*dxdwo_curr).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edxdot_tensor_dwo2( b1_edxdot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edxdot_tensor_dwodp( b1_edxdot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edxdot, 0) += b1_edxdot_tensor_dwo2;
        d2zdwodp_plus.chip(edxdot, 0) += b1_edxdot_tensor_dwodp;

        b1_edydot_matrix_dwo2  = ((d2edydot_dx_curr2(state_curr, timestep)*dxdwo_curr).transpose()*dxdwo_curr).transpose();
        b1_edydot_matrix_dwodp = ((d2edydot_dx_curr2(state_curr, timestep)*dxdp_curr).transpose()*dxdwo_curr).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edydot_tensor_dwo2( b1_edydot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b1_edydot_tensor_dwodp( b1_edydot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edydot, 0) += b1_edydot_tensor_dwo2;
        d2zdwodp_plus.chip(edydot, 0) += b1_edydot_tensor_dwodp;

        // only uses plant state in derivative
        // d2hdx2_plus_tensor = d2hdx2_plus(state_plus, time, timestep);
        // b3 = d2hdx2_plus_tensor.contract(dxdwo_plus_tensor, contract_dims).eval();
        // b3_dwo2  = b3.contract(dxdwo_plus_tensor, contract_dims).eval();
        // b3_dwodp = b3.contract(dxdp_plus_tensor, contract_dims).eval();

        // d2zdwo2_plus += b3_dwo2;
        // d2zdwodp_plus += b3_dwodp;

        b3_edx_matrix_dwo2  = ((d2edx_dx_plus2(state_plus, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_edx_matrix_dwodp = ((d2edx_dx_plus2(state_plus, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edx_tensor_dwo2( b3_edx_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edx_tensor_dwodp( b3_edx_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edx, 0) += b3_edx_tensor_dwo2;
        d2zdwodp_plus.chip(edx, 0) += b3_edx_tensor_dwodp;

        b3_edy_matrix_dwo2  = ((d2edy_dx_plus2(state_plus, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_edy_matrix_dwodp = ((d2edy_dx_plus2(state_plus, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edy_tensor_dwo2( b3_edy_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edy_tensor_dwodp( b3_edy_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edy, 0) += b3_edy_tensor_dwo2;
        d2zdwodp_plus.chip(edy, 0) += b3_edy_tensor_dwodp;

        b3_edxdot_matrix_dwo2  = ((d2edxdot_dx_plus2(state_plus, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_edxdot_matrix_dwodp = ((d2edxdot_dx_plus2(state_plus, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edxdot_tensor_dwo2( b3_edxdot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edxdot_tensor_dwodp( b3_edxdot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edxdot, 0) += b3_edxdot_tensor_dwo2;
        d2zdwodp_plus.chip(edxdot, 0) += b3_edxdot_tensor_dwodp;

        b3_edydot_matrix_dwo2  = ((d2edydot_dx_plus2(state_plus, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_edydot_matrix_dwodp = ((d2edydot_dx_plus2(state_plus, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edydot_tensor_dwo2( b3_edydot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_edydot_tensor_dwodp( b3_edydot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(edydot, 0) += b3_edydot_tensor_dwo2;
        d2zdwodp_plus.chip(edydot, 0) += b3_edydot_tensor_dwodp;

        b3_eix_matrix_dwo2  = ((d2eix_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_eix_matrix_dwodp = ((d2eix_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eix_tensor_dwo2( b3_eix_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eix_tensor_dwodp( b3_eix_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(eix, 0) += b3_eix_tensor_dwo2;
        d2zdwodp_plus.chip(eix, 0) += b3_eix_tensor_dwodp;

        b3_eiy_matrix_dwo2  = ((d2eiy_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_eiy_matrix_dwodp = ((d2eiy_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eiy_tensor_dwo2( b3_eiy_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eiy_tensor_dwodp( b3_eiy_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(eiy, 0) += b3_eiy_tensor_dwo2;
        d2zdwodp_plus.chip(eiy, 0) += b3_eiy_tensor_dwodp;

        b3_eixdot_matrix_dwo2  = ((d2eixdot_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_eixdot_matrix_dwodp = ((d2eixdot_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eixdot_tensor_dwo2( b3_eixdot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eixdot_tensor_dwodp( b3_eixdot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(eixdot, 0) += b3_eixdot_tensor_dwo2;
        d2zdwodp_plus.chip(eixdot, 0) += b3_eixdot_tensor_dwodp;

        b3_eiydot_matrix_dwo2  = ((d2eiydot_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_eiydot_matrix_dwodp = ((d2eiydot_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eiydot_tensor_dwo2( b3_eiydot_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_eiydot_tensor_dwodp( b3_eiydot_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(eiydot, 0) += b3_eiydot_tensor_dwo2;
        d2zdwodp_plus.chip(eiydot, 0) += b3_eiydot_tensor_dwodp;

        b3_desVelx_matrix_dwo2  = ((d2desVelx_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_desVelx_matrix_dwodp = ((d2desVelx_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_desVelx_tensor_dwo2( b3_desVelx_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_desVelx_tensor_dwodp( b3_desVelx_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(desVelX, 0) += b3_desVelx_tensor_dwo2;
        d2zdwodp_plus.chip(desVelX, 0) += b3_desVelx_tensor_dwodp;

        b3_desVely_matrix_dwo2  = ((d2desVely_dx_plus2(state_plus, time, timestep)*dxdwo_plus).transpose()*dxdwo_plus).transpose();
        b3_desVely_matrix_dwodp = ((d2desVely_dx_plus2(state_plus, time, timestep)*dxdp_plus).transpose()*dxdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_desVely_tensor_dwo2( b3_desVely_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b3_desVely_tensor_dwodp( b3_desVely_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(desVelY, 0) += b3_desVely_tensor_dwo2;
        d2zdwodp_plus.chip(desVelY, 0) += b3_desVely_tensor_dwodp;

        dhdz_plus = dhdz_plus_sparse.toDense();
        d2hdzplusdp = d2hdzplus_dp();
        Eigen::DSizes<Eigen::Index, 3> offset3(0, 0, 0); 
        Eigen::DSizes<Eigen::Index, 2> offset2(0, 0); 

        for(int zBeforey = 1; zBeforey <= desThrust; zBeforey++){
            Eigen::MatrixXd blockMatrixdhdz = dhdz_plus.block(zBeforey, 0, 1, zBeforey);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrixdhdz.data(), zBeforey);

            Eigen::DSizes<Eigen::Index, 3> extent_dwo2(zBeforey, NUM_STATES, NUM_STATES);
            Eigen::DSizes<Eigen::Index, 3> extent_dwodp(zBeforey, NUM_STATES, NUM_PARAMETERS);
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp_dwo2 = dhdz_plus_block_tensor.contract(d2zdwo2_plus.slice(offset3, extent_dwo2), contract_dims0).eval();
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp_dwodp = dhdz_plus_block_tensor.contract(d2zdwodp_plus.slice(offset3, extent_dwodp), contract_dims0).eval();
            
            Eigen::DSizes<Eigen::Index, 3> offsetblock(zBeforey, 0, 0);
            Eigen::DSizes<Eigen::Index, 3> extentblock(1, zBeforey, NUM_PARAMETERS);
            Eigen::DSizes<Eigen::Index, 2> extent2(zBeforey, NUM_STATES);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp2 = d2hdzplusdp.slice(offsetblock, extentblock).contract(dzdwo_plus_tensor.slice(offset2, extent2), contract_dims);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp3 = tmp2.shuffle(Eigen::array<int,3>{0,2,1});
            
            d2zdwo2_plus.chip(zBeforey, 0) += tmp_dwo2; 
            d2zdwodp_plus.chip(zBeforey, 0) += tmp_dwodp + tmp3.chip(0, 0);        
        }  

        // d2ydwo2_plus // d2ydwodp_plus
        dgdx_plus = dgdx_plus_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdx_plus_tensor(dgdx_plus.data(), NUM_Y_STATES, NUM_PLANT_STATES);
        c1_dwo2  = dgdx_plus_tensor.contract(d2xdwo2_plus, contract_dims).eval();
        c1_dwodp = dgdx_plus_tensor.contract(d2xdwodp_plus, contract_dims).eval();

        dgdz_plus = dgdz_plus_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dgdz_plus_tensor(dgdz_plus.data(), NUM_Y_STATES, NUM_Z_STATES);
        c2_dwo2  = dgdz_plus_tensor.contract(d2zdwo2_plus, contract_dims).eval();
        c2_dwodp = dgdz_plus_tensor.contract(d2zdwodp_plus, contract_dims).eval();

        d2gdx2_plus = d2gdx2(state_plus);
        c3 = d2gdx2_plus.contract(dxdwo_plus_tensor, contract_dims).eval();
        c3_dwo2  = c3.contract(dxdwo_plus_tensor, contract_dims).eval();
        c3_dwodp = c3.contract(dxdp_plus_tensor, contract_dims).eval();

        d2gdxplusdp = d2gdxplus_dp(state_plus);
        c4 = d2gdxplusdp.contract(dxdwo_plus_tensor, contract_dims);
        c5 = c4.shuffle(Eigen::array<int,3>{0,2,1});

        d2gdzplusdp = d2gdzplus_dp(state_plus);
        c6 = d2gdzplusdp.contract(dzdwo_plus_tensor, contract_dims);
        c7 = c6.shuffle(Eigen::array<int,3>{0,2,1});

        d2ydwo2_plus = c1_dwo2 + c2_dwo2 + c3_dwo2;
        d2ydwodp_plus = c1_dwodp + c2_dwodp + c3_dwodp + c5 + c7;

        // d2zdwo2_plus // d2zdwodp_plus states after y
        dhdy_plus = dhdy_plus_sparse.toDense();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> dhdy_plus_tensor(dhdy_plus.data(), NUM_Z_STATES, NUM_Y_STATES);
        b8_dwo2  = dhdy_plus_tensor.contract(d2ydwo2_plus, contract_dims).eval();
        b8_dwodp = dhdy_plus_tensor.contract(d2ydwodp_plus, contract_dims).eval();

        d2hdydp_plus = d2hdydp();
        b9 = d2hdydp_plus.contract(dydwo_plus_tensor, contract_dims).eval();
        b10 = b9.shuffle(Eigen::array<int,3>{0,2,1});
        d2zdwo2_plus+= b8_dwo2;
        d2zdwodp_plus+= b8_dwodp + b10;

        for(int zAftery = eiphi; zAftery < NUM_Z_STATES; zAftery ++){
            Eigen::MatrixXd blockMatrixdhdz = dhdz_plus.block(zAftery, 0, 1, zAftery);
            Eigen::TensorMap<Eigen::Tensor<double, 1, Eigen::ColMajor>> dhdz_plus_block_tensor(blockMatrixdhdz.data(), zAftery);

            Eigen::DSizes<Eigen::Index, 3> extent_dwo2(zAftery, NUM_STATES, NUM_STATES);
            Eigen::DSizes<Eigen::Index, 3> extent_dwodp(zAftery, NUM_STATES, NUM_PARAMETERS);
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp_dwo2 = dhdz_plus_block_tensor.contract(d2zdwo2_plus.slice(offset3, extent_dwo2), contract_dims0).eval();
            Eigen::Tensor<double, 2, Eigen::ColMajor> tmp_dwodp = dhdz_plus_block_tensor.contract(d2zdwodp_plus.slice(offset3, extent_dwodp), contract_dims0).eval();
            
            Eigen::DSizes<Eigen::Index, 3> offsetblock(zAftery, 0, 0);
            Eigen::DSizes<Eigen::Index, 3> extentblock(1, zAftery, NUM_PARAMETERS);
            Eigen::DSizes<Eigen::Index, 2> extent2(zAftery, NUM_STATES);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp2 = d2hdzplusdp.slice(offsetblock, extentblock).contract(dzdwo_plus_tensor.slice(offset2, extent2), contract_dims);
            Eigen::Tensor<double, 3, Eigen::ColMajor> tmp3 = tmp2.shuffle(Eigen::array<int,3>{0,2,1});
            
            d2zdwo2_plus.chip(zAftery, 0) += tmp_dwo2; 
            d2zdwodp_plus.chip(zAftery, 0) += tmp_dwodp + tmp3.chip(0, 0);     
        }

        // d2hdz_plus2 = d2hdzplus2(state_plus);
        // b12 = d2hdz_plus2.contract(dzdwo_plus_tensor, contract_dims).eval();
        // b12_dwo2  = b12.contract(dzdwo_plus_tensor, contract_dims).eval();
        // b12_dwodp = b12.contract(dzdp_plus_tensor, contract_dims).eval();

        // d2zdwo2_plus += b12_dwo2;
        // d2zdwodp_plus += b12_dwodp;

        b12_ft_matrix_dwo2  = ((d2ft_dz_plus2(state_plus)*dzdwo_plus).transpose()*dzdwo_plus).transpose();
        b12_ft_matrix_dwodp = ((d2ft_dz_plus2(state_plus)*dzdp_plus).transpose()*dzdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_ft_tensor_dwo2( b12_ft_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_ft_tensor_dwodp( b12_ft_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(ft, 0) += b12_ft_tensor_dwo2;
        d2zdwodp_plus.chip(ft, 0) += b12_ft_tensor_dwodp;
        
        b12_tx_matrix_dwo2  = ((d2tx_dz_plus2(state_plus)*dzdwo_plus).transpose()*dzdwo_plus).transpose();
        b12_tx_matrix_dwodp = ((d2tx_dz_plus2(state_plus)*dzdp_plus).transpose()*dzdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_tx_tensor_dwo2( b12_tx_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_tx_tensor_dwodp( b12_tx_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(tx, 0) += b12_tx_tensor_dwo2;
        d2zdwodp_plus.chip(tx, 0) += b12_tx_tensor_dwodp;

        b12_ty_matrix_dwo2  = ((d2ty_dz_plus2(state_plus)*dzdwo_plus).transpose()*dzdwo_plus).transpose();
        b12_ty_matrix_dwodp = ((d2ty_dz_plus2(state_plus)*dzdp_plus).transpose()*dzdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_ty_tensor_dwo2( b12_ty_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_ty_tensor_dwodp( b12_ty_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(ty, 0) += b12_ty_tensor_dwo2;
        d2zdwodp_plus.chip(ty, 0) += b12_ty_tensor_dwodp;

        b12_tz_matrix_dwo2  = ((d2tz_dz_plus2(state_plus)*dzdwo_plus).transpose()*dzdwo_plus).transpose();
        b12_tz_matrix_dwodp = ((d2tz_dz_plus2(state_plus)*dzdp_plus).transpose()*dzdwo_plus).transpose();
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_tz_tensor_dwo2( b12_tz_matrix_dwo2.data(), NUM_STATES, NUM_STATES);
        Eigen::TensorMap<Eigen::Tensor<double, 2,Eigen::ColMajor>> b12_tz_tensor_dwodp( b12_tz_matrix_dwodp.data(), NUM_STATES, NUM_PARAMETERS);
        d2zdwo2_plus.chip(tz, 0) += b12_tz_tensor_dwo2;
        d2zdwodp_plus.chip(tz, 0) += b12_tz_tensor_dwodp;

        // update step
        dxdp_curr = dxdp_plus; dzdp_curr = dzdp_plus; dydp_curr = dydp_plus;
        d2xdwo2_curr = d2xdwo2_plus; d2zdwo2_curr = d2zdwo2_plus; d2ydwo2_curr = d2ydwo2_plus;
        d2xdwodp_curr = d2xdwodp_plus; d2zdwodp_curr = d2zdwodp_plus; d2ydwodp_curr = d2ydwodp_plus;
    }

    std::chrono::time_point end = std::chrono::steady_clock::now();
    std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    m_logger << "Elapsed Time d2w: " << elapsed.count() << " us" << std::endl;

    return {{d2xdwo2_curr, d2zdwo2_curr, d2ydwo2_curr}, {d2xdwodp_curr, d2zdwodp_curr, d2ydwodp_curr}};    
}

Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> DroneTrajectory::calc_dG(dwdwo ts, d2w const & d2w, G_tp gtp)
{
    std::chrono::time_point start = std::chrono::steady_clock::now();

    Eigen::Matrix<double, NUM_STATES, NUM_STATES> dwdwo_matrix;
    dwdwo_matrix << ts.dxdwo, ts.dzdwo, ts.dydwo;

    Eigen::Tensor<double, 3> tmp_dwo2  = d2w.dwo2.d2xdwo2.concatenate(d2w.dwo2.d2zdwo2, 0);
    Eigen::Tensor<double, 3> secondOrderTraj = tmp_dwo2.concatenate(d2w.dwo2.d2ydwo2, 0);

    Eigen::Tensor<double, 3> tmp_dwodp = d2w.dwodp.d2xdwodp.concatenate(d2w.dwodp.d2zdwodp, 0);
    Eigen::Tensor<double, 3> secondOrderTrajParams = tmp_dwodp.concatenate(d2w.dwodp.d2ydwodp, 0);

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
    return dG;
}