#include "DroneTrajectory.h"
#include <iostream>
#include "Logger.h"
typedef Eigen::Vector<double, 4> test_state;
typedef Eigen::Vector<double, 2> test_x_state;
enum test_state_enum{test_x1, test_x2, test_y, test_z};

struct test_timestep
{
    test_state state;
    bool stable = true;
};

struct test_simResults{
    std::vector<double> time;
    std::vector<test_state> stateProgression;
    bool stable = true;
    bool converged = false;
};

// dynamics and updated control state
test_x_state f(test_state state) 
{
    return {state(test_x2), -state(test_x1) - 0.1*state(test_x2) + state(test_z)};
}

Eigen::Matrix2d dfdx()
{
    Eigen::Matrix2d M;
    M << 0, 1,
        -1, -0.1;
    return M;
}

Eigen::Matrix<double, 2, 1> dfdz(){
    Eigen::Matrix<double, 2, 1> M;
    M << 0, 1;
    return M;
}

test_x_state H(test_state prev, test_x_state x_guess, double timestep) 
{
    test_x_state x_prev = {prev(test_x1), prev(test_x2)};
    test_x_state fprev = f(prev);
    test_state guess = prev;
    guess(test_x1) = x_guess(test_x1);
    guess(test_x2) = x_guess(test_x2);
    test_x_state fguess = f(guess);
    return x_prev - x_guess + timestep/2*(fprev + fguess);
}

Eigen::Matrix2d DH(test_x_state x_guess, double timestep) 
{
    Eigen::Matrix2d dh = -1*Eigen::Matrix2d::Identity();
    Eigen::Matrix2d dfdx_plus = dfdx();
    return dh + timestep/2*dfdx_plus;
}

test_timestep simulateTimestep(test_state prev, double time, double timestep, Logger & log) 
{
    double tol = 1e-12;
    test_state next = prev;
    test_x_state x_prev = {prev(test_x1), prev(test_x2)};
    test_x_state x_guess = x_prev + timestep*f(prev);
    int count = 0;
    int max_iterations = 100;
    bool stable = true;
    for(; count < max_iterations; count++ ){
        test_x_state h = H(prev, x_guess, timestep);
        if (!h.allFinite()) {
            stable = false;
            log.warn(std::string("simulateTimestep: non-finite residual h"));
            break;
        }
        if (h.norm() < tol) {
            break;
        }
        Eigen::MatrixX<double> dh = DH(x_guess, timestep);
        if (!dh.allFinite()) {
            stable = false;
            log.warn(std::string("simulateTimestep: non-finite Jacobian DH"));
            break;
        }

        Eigen::FullPivLU<Eigen::MatrixXd> lu(dh);
        if (!lu.isInvertible()) {
            stable = false;
            log.warn(std::string("simulateTimestep: DH is singular or ill-conditioned"));
            break;
        }

        x_guess = x_guess - lu.solve(h);
        if (!x_guess.allFinite()) {
            stable = false;
            log.warn(std::string("simulateTimestep: guess.plant became non-finite"));
            break;
        }
    }

    next(test_x1) = x_guess(test_x1);
    next(test_x2) = x_guess(test_x2);
    if(next(test_x1) >  0){
        next(test_y) = std::sin(next(test_x1));
    } else {
        next(test_y) = 0;
    }
    next(test_z) = -next(test_y);

    stable = stable && count < max_iterations;   

    if (!stable) {
        return {prev, false};
    }
    return {next, true};
}


test_simResults test_simple_trajectory(test_state initialState, double finalTime, Logger & log) 
{
    test_simResults simResults;
    simResults.stateProgression.push_back(initialState);
    double time = 0;
    double timestep = 1e-3;
    simResults.time.push_back(time);
    while(time < finalTime && simResults.stable){
        test_state prev = simResults.stateProgression.back();
        test_timestep state1 = simulateTimestep(prev, time, timestep, log);
        simResults.stateProgression.push_back(state1.state);
        simResults.stable &= state1.stable;

        time += timestep;        
        simResults.time.push_back(time);

    }
    return simResults;
}

std::vector<Eigen::Matrix4d> trajSens(test_simResults const & simResults) 
{
    // std::chrono::time_point start = std::chrono::steady_clock::now();

    const int iterations = simResults.time.size();
    std::vector<Eigen::Matrix4d> ts;
    ts.reserve(iterations);
    Eigen::Matrix2d I;
    I.setIdentity();
    Eigen::Matrix4d initial_dwdwo;
    initial_dwdwo.setIdentity();
    ts.push_back(initial_dwdwo);

    // iterating trajectory sensitivity
    for(int i = 1; i < iterations; i++)
    {
        double timestep = simResults.time[i] - simResults.time[i-1];
        // dwdwo
        Eigen::Matrix4d dwdwo = ts[i-1];
        Eigen::Matrix4d dwdwo_plus;

        // dfdx constant
        Eigen::Matrix2d dfdxMatrix = dfdx();
        Eigen::Matrix<double, 2, 1> dfdzMatrix = dfdz();

        Eigen::Matrix2d A = I - (timestep / 2.0) * dfdxMatrix;
        Eigen::Matrix<double, 2, 4> B = dwdwo.topRows(2) + (timestep / 2.0)*(dfdxMatrix*dwdwo.topRows(2) + (dfdzMatrix+dfdzMatrix)*dwdwo.row(test_z));
        Eigen::PartialPivLU<Eigen::Matrix2d> solver;
        solver.compute(A);
        dwdwo_plus.topRows(2) = solver.solve(B);

        if (simResults.stateProgression.at(i).coeff(test_x1) > 0) {
            dwdwo_plus.row(test_y) = std::cos(simResults.stateProgression.at(i).coeff(test_x1)) * dwdwo_plus.row(test_x1);
            dwdwo_plus.row(test_z) = -dwdwo_plus.row(test_y);
        } else {
            dwdwo_plus.row(test_y).setZero();
            dwdwo_plus.row(test_z) = 0.5*dwdwo.row(test_z);
        }
        ts.push_back(dwdwo_plus);
    }

    return ts;
}

std::vector<Eigen::Matrix4d> trajSensTest(test_state initialState, double finalTime, Logger & log)
{
    double delta = 1e-5;
    test_simResults simResults = test_simple_trajectory(initialState, finalTime, log);
    int numIterations = simResults.time.size();
    std::vector<Eigen::Matrix4d> ts(numIterations);
    
    test_state plus;
    test_state minus;
            
    for(int i = 0; i < 4; i++) {
        test_state testPlusState = initialState;
        test_state testMinusState = initialState;
        testPlusState(i) += delta;
        testMinusState(i) -= delta;

        test_simResults plusSimResults = test_simple_trajectory(testPlusState, finalTime, log); 
        test_simResults minusSimResults = test_simple_trajectory(testMinusState, finalTime, log);
        
        for(int t = 0; t < numIterations; t++)
        {
            plus = plusSimResults.stateProgression.at(t);
            minus = minusSimResults.stateProgression.at(t);
            test_state delta_Trajsens = 1/(2*delta)*(plus-minus);
            ts.at(t).col(i) = delta_Trajsens;
        }
    }
    return ts;
}


void splot(test_simResults simResults)
{
    
    // Use full path if PATH is not set
    const char* cmd = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent";
    FILE* gp = _popen(cmd, "w");
    if (!gp) {
        printf("Error: Could not open gnuplot!\n");
        return;
    }

    // Send gnuplot commands
    fprintf(gp, "set title 'Plant State %d' \n", test_x1);
    fprintf(gp, "set xlabel 'time'\n");
    fprintf(gp, "set ylabel 'plantSate'\n");
    fprintf(gp, "plot '-' with lines title 'trajectory'\n");

    // Send the data
    for (int i = 0; i < (int) simResults.time.size(); i++) {
        fprintf(gp, "%f %f\n", simResults.time.at(i), simResults.stateProgression.at(i)(test_x1));
    }

    fprintf(gp, "e\n");  // 'e' ends the data section
    _pclose(gp);  // close gnuplot
    
}

int main()
{
    Logger log("./build/log.txt");
    test_state initialState = {1, 0, 0, 0};
    double finalTime = 50;
    test_simResults simResults = test_simple_trajectory(initialState, finalTime, log);
    splot(simResults);
    std::vector<Eigen::Matrix4d> ts = trajSens(simResults) ;
    std::vector<Eigen::Matrix4d> tstest = trajSensTest( initialState,  finalTime,  log);
    for(int i = 0; i < simResults.time.size(); i++)
    {
        // log << "i: " << i << " test_x1 " << simResults.stateProgression.at(i)(test_x1) << " diff " << (ts.at(i)-tstest.at(i)).cwiseAbs().maxCoeff() << std::endl;
        log << "i: "  << i << " " << (ts.at(i).row(test_x1)).norm() << " " <<  (tstest.at(i).row(test_x1)).norm() << " " 
                        << (ts.at(i).row(test_x2)).norm() << " " <<  (tstest.at(i).row(test_x2)).norm() << " "
                        << (ts.at(i).row(test_y)).norm() << " " <<  (tstest.at(i).row(test_y)).norm() << " "
                        << (ts.at(i).row(test_z)).norm() << " " <<  (tstest.at(i).row(test_z)).norm() << " "<< std::endl;
    }
    std::cout << ":D" << std::endl;
    return 0;
}





