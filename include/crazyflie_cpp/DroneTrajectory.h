#ifndef DRONE_TRAJECTORY_H
#define DRONE_TRAJECTORY_H

#define _USE_MATH_DEFINES // Must be defined before including cmath or math.h
#include <cmath>
#include <vector>
#include <array>
#include <Eigen/Dense>

#include "Logger.h"

#define NUM_PLANT_STATES 12
#define NUM_CTRL_STATES 12
#define NUM_ALGE_STATES 4
#define NUM_DIST_STATES 6
#define NUM_REF_STATES 3

struct PIDParameters 
{
    double kp = 0;
    double ki = 0;
    double kd = 0;
    double integration_limit = 0;
    
};

struct PIDstate
{
    double kp_error = 0;
    double ki_error = 0;
    double kd_error = 0;
    double prev_sig = 0;
};

struct PIDCtrllers
{
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/platform/interface/platform_defaults_tag.h#L46
    PIDParameters posX = {2, 0, 0};
    PIDParameters posY = {2, 0, 0};
    PIDParameters posZ = {2, 0.5, 0};
    
    PIDParameters velX = {25, 1, 0};
    PIDParameters velY = {25, 1, 0};
    PIDParameters velZ = {25, 15, 0};

    PIDParameters attX = {6, 3, 0, 20.0};
    PIDParameters attY = {6, 3, 0, 20.0};
    PIDParameters attZ = {6, 1, 0, 260.0};

    PIDParameters attRateX = {250.0, 500, 0, 33.3};
    PIDParameters attRateY = {250.0, 500, 0, 33.3};
    PIDParameters attRateZ = {120, 16.7, 0, 166.7};
};

struct DroneParameters 
{
    double length = 0.046; // [m]
    double Ix = 16.571710e-6;
    double Iy = 16.655602e-6;
    double Iz = 29.261652e-6;
    double mass = 0.034; // [kg]
    double g = 9.81;

    // used for implementing motor controller saturation
    double maxVel = 4;
    double maxAngVel = 8;
    double maxAcc = 25;
    double maxAngAcc = 20;

    // https://arxiv.org/pdf/2512.14450
    double kf = 3.72e-8;
    double km = 7.73e-11;
    // https://github.com/IMRCLab/crazyswarm2/blob/bd54392f91d5c3aa29d5c170b11359767ab105d1/crazyflie_sim/crazyflie_sim/backend/data/dynobench/crazyflie2.yaml#L11 
    // says this ratio should be 0.006
    // and right now its 0.0021

    // x y z psi theta phi xdot ydot zdot p q r 
    double numPlantStates = 12;
    // 2 pids * 3 translational directions * 2 states per error
    double numCtrlStates = 12; 

    double pid_vel_roll_max = 20.0;
    double pid_vel_pitch_max = 20.0;
};

struct CtrlOut
{
    Eigen::Vector<double, NUM_ALGE_STATES> algeStates;
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> ctrlStates;
};

struct SystemState
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant;
    Eigen::Vector<PIDstate, NUM_CTRL_STATES> ctrl;
    Eigen::Vector<double, NUM_ALGE_STATES> alge;
    bool stable = true;
};

struct SimResults{
    std::vector<double> time;
    std::vector<SystemState> stateProgression;
    bool stable;
};

enum plantIndex{x, y, z, phi, theta, psi, xdot, ydot, zdot, p, q, r};
enum ctrlIndex {posX, posY, posZ, velX, velY, velZ, attX, attY, attZ, attRateX, attRateY, attRateZ};
enum algeIndex {ft, tx, ty, tz};
enum distIndex {Fwx, Fwy, Fwz, Twx, Twy, Twz};
enum refIndex  {refx, refy, refz};

class DroneTrajectory 
{
    Logger & m_logger;
    PIDCtrllers m_ctrlParams;
    // since attitude pid and attitude rate pid run at a 5x frequency
    // assume their dynamics are fast enough such that they behave "instantly" relative to outer loop controller
    DroneParameters m_droneParams;
    // input is time
    // direction based on distIndex enum
    std::array<double(*)(double), NUM_DIST_STATES> m_dist; 
    std::array<double(*)(double), NUM_REF_STATES> m_ref; 
    double m_simTimestep; // [s]
    double m_finalTime; // [s]
    

    CtrlOut CascadedPIDController(Eigen::Vector<double, NUM_PLANT_STATES> plantState, Eigen::Vector<PIDstate, NUM_CTRL_STATES> ctrlState, double time);
    SystemState simulateTimestep(SystemState prev, double time);
    Eigen::Vector<double, NUM_PLANT_STATES> H(SystemState prev, Eigen::Vector<double, NUM_PLANT_STATES> guess, double time);
    Eigen::MatrixX<double> DH(SystemState state);
    Eigen::MatrixX<double> dfdx(SystemState state);
    Eigen::Vector<double, NUM_PLANT_STATES> f(SystemState state, double time);

    public:
        DroneTrajectory( 
            Logger & log, 
            std::array<double(*)(double), NUM_DIST_STATES> const& dist,
            std::array<double(*)(double), NUM_REF_STATES> const& ref,
            PIDCtrllers ctrlParams = {},
            DroneParameters droneParameters = {},
            double simTimestep = 2e-3, double finalTime = 10);
        
        SimResults Trajectory(SystemState initialState); 
};

#endif