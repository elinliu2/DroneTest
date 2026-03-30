#ifndef DRONE_TRAJECTORY_H
#define DRONE_TRAJECTORY_H

#define _USE_MATH_DEFINES // Must be defined before including cmath or math.h
#include <cmath>
#include <vector>
#include <array>
#include <Eigen/Dense>

#include "Logger.h"
#include "LowPassFilter.h"

#define NUM_PLANT_STATES 12
#define NUM_CTRL_STATES 12
#define NUM_ALGE_STATES 4
#define NUM_DIST_STATES 6
#define NUM_REF_STATES 4

struct PIDParameters 
{
    double kp = 0;
    double ki = 0;
    double kd = 0;
    double integration_limit = 0;
    bool lpf_en = false;
};

struct PIDstate
{
    double kp_error = 0;
    double ki_error = 0;
    double kd_error = 0;
    double prev_sig = 0;
    LowPassFilter lpf;
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

    // // https://giuseppesilano.net/publications/rosChapter19.pdf
    // I tried kf and km gains from here but they didn't work

    // https://arxiv.org/pdf/2512.14450
    double kf = 3.72e-8;
    double km = 7.73e-11;
    
    // x y z psi theta phi xdot ydot zdot p q r 
    double numPlantStates = 12;
    // 2 pids * 3 translational directions * 2 states per error
    double numCtrlStates = 12; 

    // these are in degrees
    double pid_vel_roll_max = 20.0; 
    double pid_vel_pitch_max = 20.0;
};

struct CtrlOut
{
    Eigen::Vector<double, NUM_ALGE_STATES> algeStates;
    std::array<PIDstate, NUM_CTRL_STATES> ctrlStates;
};

struct SystemState
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant;
    std::array<PIDstate, NUM_CTRL_STATES> ctrl;
    Eigen::Vector<double, NUM_ALGE_STATES> alge;
    bool stable = true;
};

struct SimResults{
    std::vector<double> time;
    std::vector<SystemState> stateProgression;
    bool stable = true;
    bool converged = false;
};

std::array<PIDParameters, NUM_CTRL_STATES> inline defaultPIDParameters(){
    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/platform/interface/platform_defaults_tag.h#L46
    PIDParameters posXpid = {2, 0, 0};
    PIDParameters posYpid = {2, 0, 0};
    PIDParameters posZpid = {2, 0.5, 0};

    PIDParameters velXpid = {25, 1, 0};
    PIDParameters velYpid = {25, 1, 0};
    PIDParameters velZpid = {25, 15, 0};

    PIDParameters rollpid  = {6, 3, 0, 20.0};
    PIDParameters pitchpid = {6, 3, 0, 20.0};
    PIDParameters yawpid   = {6, 1, 0.35, 360.0};

    PIDParameters rollRatepid  = {250.0, 500, 2.5, 33.3, true};
    PIDParameters pitchRatepid = {250.0, 500, 2.5, 33.3, true};
    PIDParameters yawRatepid   = {120.0, 16.7, 0, 166.7};
    return {posXpid, posYpid, posZpid, velXpid, velYpid, velZpid, rollpid, pitchpid, yawpid, rollRatepid, pitchRatepid, yawRatepid};
}

enum plantIndex{x, y, z, phi, theta, psi, xdot, ydot, zdot, p, q, r};
enum ctrlIndex {posX, posY, posZ, velX, velY, velZ, roll, pitch, yaw, rollRate, pitchRate, yawRate};
enum algeIndex {ft, tx, ty, tz};
enum distIndex {Fwx, Fwy, Fwz, Twx, Twy, Twz};
enum refIndex  {refx, refy, refz, refyaw};

class DroneTrajectory 
{
    Logger & m_logger;
    std::array<PIDParameters, NUM_CTRL_STATES> m_ctrlParams;
    DroneParameters m_droneParams;
    // input is time
    // direction based on distIndex enum
    std::array<double(*)(double), NUM_DIST_STATES> m_dist; 
    std::array<double(*)(double), NUM_REF_STATES> m_ref; 
    double m_simTimestep; // [s]
    double m_finalTime; // [s]
    

    CtrlOut CascadedPIDController(Eigen::Vector<double, NUM_PLANT_STATES> plantState, std::array<PIDstate, NUM_CTRL_STATES> ctrlState, double time, double timestep);
    SystemState simulateTimestep(SystemState prev, double time, double timestep);
    Eigen::Vector<double, NUM_PLANT_STATES> H(SystemState prev, Eigen::Vector<double, NUM_PLANT_STATES> guess, double time, double timestep);
    Eigen::MatrixX<double> DH(SystemState state, double timestep);
    Eigen::MatrixX<double> dfdx(SystemState state);
    Eigen::Vector<double, NUM_PLANT_STATES> f(SystemState state, double time);
    bool isConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time);
    bool isNotConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time, std::array<PIDstate, NUM_CTRL_STATES> ctrlState);

    public:
        DroneTrajectory( 
            Logger & log, 
            std::array<double(*)(double), NUM_DIST_STATES> const& dist,
            std::array<double(*)(double), NUM_REF_STATES> const& ref,
            std::array<PIDParameters, NUM_CTRL_STATES> ctrlParams = defaultPIDParameters(),
            DroneParameters droneParameters = {},
            double simTimestep = 1e-3, double finalTime = 250);
        
        SimResults Trajectory(SystemState initialState); 
};

#endif