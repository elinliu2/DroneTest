#include <vector>
#include <Eigen/Dense>

#define NUM_PLANT_STATES 12
#define NUM_CTRL_STATES 12
#define NUM_ALGE_STATES 4
#define NUM_DIST_STATES 6
#define NUM_REF_STATES 3

struct PIDParameters 
{
    double kp = 6.0;
    double kd = 6.0;
    double ki = 6.0;
};

struct PIDstate
{
    double kp_error = 0;
    double ki_error = 0;
    double kd_error = 0;
};

struct PIDCtrllers
{
    PIDParameters posX;
    PIDParameters posY;
    PIDParameters posZ;
    
    PIDParameters velX;
    PIDParameters velY;
    PIDParameters velZ;

    PIDParameters attX;
    PIDParameters attY;
    PIDParameters attZ;

    PIDParameters attRateX;
    PIDParameters attRateY;
    PIDParameters attRateZ;
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
    PIDCtrllers m_ctrlParams;
    // since attitude pid and attitude rate pid run at a 5x frequency
    // assume their dynamics are fast enough such that they behave "instantly" relative to outer loop controller
    DroneParameters m_droneParameters;
    // input is time
    // direction based on distIndex enum
    std::array<double(*)(double), NUM_DIST_STATES> m_dist; 
    double m_simTimestep = 0.001; // [s]
    double m_finalTime = 100; // [s]

    public:
        DroneTrajectory(
            PIDCtrllers ctrlParams,
            DroneParameters droneParameters, 
            std::array<double(*)(double), 6> const& dist,
            double simTimestep, double finalTime);
        
        SimResults Trajectory(std::array<double(*)(double), NUM_REF_STATES> const& ref, SystemState const& initialState); 
};