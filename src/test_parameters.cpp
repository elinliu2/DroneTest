#include "Crazyflie.h"
#include <iostream>
#include <cstdint>

#define NUM_PARAMETERS 36
class TestLogger : public Logger
{
    public: 
    void info(const std::string& msg) override {
        std::cout << "INFO: " << msg << std::endl;
    }
    void warning(const std::string& msg) override {
        std::cout << "WARN: " << msg << std::endl;
    }
    void error(const std::string& msg) override {
        std::cout << "ERROR: " << msg << std::endl;
    }
};

int main()
{
    std::ofstream file("./build/parameterlog.txt");
    std::string link_uri = "radio://0/80/2M/E7E7E7E7E7";
    TestLogger testLogger; 
    Crazyflie crazyflie(link_uri, testLogger);
    crazyflie.requestParamToc(true);
    // 164,6,0,pid_attitude,roll_kp
    // float gain = crazyflie.getParam<float>(164);
    // std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;

    std::vector<float> nominalPIDgains = {  2, 1e-10, 1e-10,
                                            2, 1e-10, 1e-10,
                                            2, 0.5, 1e-10,

                                            25, 1, 1e-10,
                                            25, 1, 1e-10,
                                            25, 15, 1e-10,

                                            6, 3, 1e-10, 1e-10,
                                            6, 3, 1e-10, 1e-10,
                                            6, 1, 0.35, 1e-10,

                                            250.0, 500, 2.5,
                                            250.0, 500, 2.5,
                                            120.0, 16.7, 1e-10 };
                                        
    std::vector<float> scaling_factors = {  0.999881, 1, 1,
                                            0.998236, 1, 1,
                                            0.969074, 0.997806, 1,

                                            0.9999, 0.999999, 1,
                                            0.998259, 0.999989, 1,
                                            0.874869, 0.994998, 1,

                                            0.957454, 0.999577, 1,
                                            0.961113, 0.999624, 1,
                                            0.948504, 0.999708, 0.884178,

                                            0.948737, 0.998083, 0.973266,
                                            0.952399, 0.998214, 0.973419,
                                            0.107803, 0.996855, 1};  
                                      
    std::vector<std::pair<std::string, std::string>> paramKeys = {
        {"posCtlPid","xKp"},
        {"posCtlPid","xKi"},
        {"posCtlPid","xKd"},

        {"posCtlPid","yKp"},
        {"posCtlPid","yKi"},
        {"posCtlPid","yKd"},

        {"posCtlPid","zKp"},
        {"posCtlPid","zKi"},
        {"posCtlPid","zKd"},

        {"velCtlPid","vxKp"},
        {"velCtlPid","vxKi"},
        {"velCtlPid","vxKd"},

        {"velCtlPid","vyKp"},
        {"velCtlPid","vyKi"},
        {"velCtlPid","vyKd"},

        {"velCtlPid","vzKp"},
        {"velCtlPid","vzKi"},
        {"velCtlPid","vzKd"},

        {"pid_attitude","roll_kp"},
        {"pid_attitude","roll_ki"},
        {"pid_attitude","roll_kd"},

        {"pid_attitude","pitch_kp"},
        {"pid_attitude","pitch_ki"},
        {"pid_attitude","pitch_kd"},

        {"pid_attitude","yaw_kp"},
        {"pid_attitude","yaw_ki"},
        {"pid_attitude","yaw_kd"},

        {"pid_rate","roll_kp"},
        {"pid_rate","roll_ki"},
        {"pid_rate","roll_kd"},

        {"pid_rate","pitch_kp"},
        {"pid_rate","pitch_ki"},
        {"pid_rate","pitch_kd"},

        {"pid_rate","yaw_kp"},
        {"pid_rate","yaw_ki"},
        {"pid_rate","yaw_kd"},
    };
    
    std::vector<int> paramIndex(NUM_PARAMETERS);
    for(int i = 0; i < NUM_PARAMETERS; i++ ){
         const Crazyflie::ParamTocEntry* entry = crazyflie.getParamTocEntry(paramKeys.at(i).first, paramKeys.at(i).second);

        if (entry) {
            paramIndex.at(i) = entry->id;
        } else {
            std::cout << "Parameter " << i << " not found in TOC" << std::endl;
        }
    }
   

    for(int i = 0; i < NUM_PARAMETERS; i++){
        float pid_gain = nominalPIDgains.at(i) * scaling_factors.at(i);
        file << "i: " << i << " pid_gain " << pid_gain << std::endl;
        
        if (pid_gain > 0.1)
        {
            crazyflie.setParam(paramIndex.at(i), pid_gain);
        }
    }

    std::cout << "set params" << std::endl;

    // gain = crazyflie.getParam<float>(164);
    // std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;
    // float roll_ki = crazyflie.getParam<float>(165);
    // std::cout << "pid attitude roll_ki gain is: " << roll_ki << std::endl;

    // int8_t attFiltEn = crazyflie.getParam<int8_t>(177);
    // std::cout << "attFiltEn is: " << std::to_string(attFiltEn) << std::endl;
    // float attFiltCut = crazyflie.getParam<float>(178);
    // std::cout << "attFiltCut is: " << attFiltCut << std::endl;
    // int8_t rateFiltEn = crazyflie.getParam<int8_t>(191);
    // std::cout << "rateFiltEn is: " << std::to_string(rateFiltEn) << std::endl;

    // float omxFiltCut = crazyflie.getParam<float>(192);
    // std::cout << "omxFiltCut is: " << omxFiltCut << std::endl;
    // float omyFiltCut = crazyflie.getParam<float>(193);
    // std::cout << "omyFiltCut is: " << omyFiltCut << std::endl;
    
    return 0;
}