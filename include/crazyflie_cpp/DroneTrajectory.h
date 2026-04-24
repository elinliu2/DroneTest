#ifndef DRONE_TRAJECTORY_H
#define DRONE_TRAJECTORY_H

#define _USE_MATH_DEFINES // Must be defined before including cmath or math.h
#include <cmath>
#include <vector>
#include <array>
#include <Eigen/Dense>

#include "Logger.h"

#define NUM_PLANT_STATES 12
#define NUM_PIDS 12
#define NUM_Z_STATES 64
#define NUM_Y_STATES 2
#define NUM_ALGE_STATES (NUM_Y_STATES+NUM_Z_STATES)
#define NUM_STATES (NUM_PLANT_STATES+NUM_ALGE_STATES)
#define NUM_DIST_STATES 6
#define NUM_REF_STATES 4
#define NUM_PID_STATES 3

struct PIDParameters 
{
    double kp = 0;
    double ki = 0;
    double kd = 0;
    double integration_limit = 0;
    bool lpf_en = false;
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

struct SystemState
{
    Eigen::Vector<double, NUM_PLANT_STATES> plant;
    Eigen::Vector<double, NUM_ALGE_STATES> alge;
};

struct Timestep
{
    SystemState state;
    bool stable = true;
};

struct SimResults{
    std::vector<double> time;
    std::vector<SystemState> stateProgression;
    bool stable = true;
    bool converged = false;
};

std::array<PIDParameters, NUM_PIDS> inline defaultPIDParameters(){
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

enum plantIndex {x, y, z, phi, theta, psi, xdot, ydot, zdot, p, q, r};
enum ctrlIndex  {posX, posY, posZ, velX, velY, velZ, roll, pitch, yaw, rollRate, pitchRate, yawRate};
// enum algeIndex  {ft, tx, ty, tz, setp_body_x, setp_body_y, state_body_x, state_body_y, w1, w2, w3, w4, 
//                 epx, eix, edx, desVelX, epy, eiy, edy, desVelY, epz, eiz, edz, desVelZ, 
//                 state_body_vx, epxdot, eixdot, edxdot, state_body_vy, epydot, eiydot, edydot, epzdot, eizdot, edzdot, desThrust, 
//                 epphi, eiphi, edphi, desRollRate, eptheta, eitheta, edtheta, desPitchRate, eppsi, eipsi, edpsi, desYawRate, 
//                 epp, eip, edp, desRollOutput, epq, eiq, edq, desPitchOutput, epr, eir, edr, desYawOutput, 
//                 delay_1_rollRate, delay_2_rollRate, delay_1_pitchRate, delay_2_pitchRate, desRoll, desPitch};
enum algeIndex  {eix, edx, desVelX, eiy, edy, desVelY, eiz, edz, desVelZ, 
                eixdot, edxdot, eiydot, edydot, eizdot, edzdot, desThrust, 
                eiphi, edphi, desRollRate, eitheta, edtheta, desPitchRate, eipsi, edpsi, desYawRate, 
                eip, edp, desRollOutput, eiq, edq, desPitchOutput, eir, edr, desYawOutput, 
                delay_1_rollRate, delay_2_rollRate, delay_1_pitchRate, delay_2_pitchRate, ft, tx, ty, tz, w1, w2, w3, w4, desRoll, desPitch};
enum pidIndex   {kp_error, ki_error, kd_error};
enum distIndex  {Fwx, Fwy, Fwz, Twx, Twy, Twz};
enum refIndex   {refx, refy, refz, refyaw};

struct dwdwo {
    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES> dxdwo;
    Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES> dzdwo;
    Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES> dydwo;

    dwdwo(
        const Eigen::Matrix<double, NUM_PLANT_STATES, NUM_STATES>& dx,
        const Eigen::Matrix<double, NUM_Z_STATES, NUM_STATES>& dz,
        const Eigen::Matrix<double, NUM_Y_STATES, NUM_STATES>& dy
    )
        : dxdwo(dx), dzdwo(dz), dydwo(dy)
    {}

    dwdwo() : dxdwo(), dzdwo(),  dydwo() {}
};

class DroneTrajectory 
{
    Logger & m_logger;
    std::array<PIDParameters, NUM_PIDS> m_ctrlParams;
    DroneParameters m_droneParams;
    // input is time
    // direction based on distIndex enum
    std::array<double(*)(double), NUM_DIST_STATES> m_dist; 
    std::array<double(*)(double), NUM_REF_STATES> m_ref; 
    double m_simTimestep; // [s]
    double m_finalTime; // [s]
    bool m_fixedNumIterations; 

    double lpf_a1;
    double lpf_a2;
    double lpf_b0;
    double lpf_b1;
    double lpf_b2;

    // https://github.com/bitcraze/crazyflie-firmware/blob/master/src/modules/src/controller/position_controller_pid.c#L260
    double m_thrustScale = 1000;
    double m_thrustBase = 36000;
    // need to convert motor pwm signal to angular velocity
    // want thrustBase = 36000 PWM == hover --> ft = mg
    // des_wi = sqrt(mg/(kf*4)) 
    double m_alpha = sqrt(m_droneParams.mass*m_droneParams.g/(m_droneParams.kf*4))/m_thrustBase;    

    Eigen::Vector<double, NUM_ALGE_STATES> CascadedPIDController(Eigen::Vector<double, NUM_PLANT_STATES> plantState, 
    Eigen::Vector<double, NUM_PLANT_STATES> prevPlantState,
    Eigen::Vector<double, NUM_ALGE_STATES> currAlgeStates, double time, double timestep);
    Timestep simulateTimestep(SystemState prev, double time, double timestep);
    Eigen::Vector<double, NUM_PLANT_STATES> H(SystemState prev, Eigen::Vector<double, NUM_PLANT_STATES> guess, double time, double timestep);
    Eigen::MatrixX<double> DH(SystemState state, double timestep);
    Eigen::Vector<double, NUM_PLANT_STATES> f(SystemState state, double time);
    bool isConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time);
    bool isNotConverging(SystemState state, std::array<double(*)(double), NUM_REF_STATES> const& ref, double time);

    Eigen::Matrix<double, NUM_PLANT_STATES, NUM_PLANT_STATES> dfdx(SystemState state);
    Eigen::SparseMatrix<double> dfdz( SystemState state);
    Eigen::SparseMatrix<double> dgdz(SystemState state);
    Eigen::SparseMatrix<double> dhdxPlus(SystemState state, double timestep);
    Eigen::SparseMatrix<double> dhdxCurr(double timestep);
    Eigen::SparseMatrix<double> dhdzPlus(SystemState state, double timestep);
    Eigen::SparseMatrix<double> dhdzCurr(double timestep);
    Eigen::SparseMatrix<double> dhdy();

    Eigen::Vector<double, NUM_ALGE_STATES> ImTooTiredToCareThatItsUglyAndDumb_h( 
    Eigen::Vector<double, NUM_PLANT_STATES> plantState,
    Eigen::Vector<double, NUM_PLANT_STATES> prevPlantState,
    Eigen::Vector<double, NUM_ALGE_STATES> currAlgeStates,
    double time, double timestep, int index);

    public:
        DroneTrajectory( 
            Logger & log, 
            std::array<double(*)(double), NUM_DIST_STATES> const& dist,
            std::array<double(*)(double), NUM_REF_STATES> const& ref,
            std::array<PIDParameters, NUM_PIDS> ctrlParams = defaultPIDParameters(),
            DroneParameters droneParameters = {},
            double simTimestep = 1e-3, double finalTime = 10, double sampleRate = 500, double cutoffFreq = 30, bool fixedNumIterations = true);
        
        SimResults Trajectory(SystemState initialState); 
        std::vector<dwdwo> trajSens(SimResults const & simResults);
        std::vector<dwdwo> trajSensTest(SystemState initialState);
        void dfdx_test(SystemState initialState);
        void dfdz_test(SystemState initialState);
        void dhdx_test(SystemState currState, SystemState prevState, double time, double timestep);
        void dhdz_test(SystemState currState, SystemState prevState, double time, double timestep);
};

#endif