#include "Crazyflie.h"
#include <iostream>
#include <cstdint>

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
    std::string link_uri = "radio://0/80/2M/E7E7E7E7E7";
    TestLogger testLogger; 
    Crazyflie crazyflie(link_uri, testLogger);
    crazyflie.requestParamToc(true);
    // 164,6,0,pid_attitude,roll_kp
    // float gain = crazyflie.getParam<float>(164);
    // std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;
    // crazyflie.setParam(164, 3.5f);
    // gain = crazyflie.getParam<float>(164);
    // std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;
    // float roll_ki = crazyflie.getParam<float>(165);
    // std::cout << "pid attitude roll_ki gain is: " << roll_ki << std::endl;
    int8_t attFiltEn = crazyflie.getParam<int8_t>(177);
    std::cout << "attFiltEn is: " << std::to_string(attFiltEn) << std::endl;
    float attFiltCut = crazyflie.getParam<float>(178);
    std::cout << "attFiltCut is: " << attFiltCut << std::endl;
    int8_t rateFiltEn = crazyflie.getParam<int8_t>(191);
    std::cout << "rateFiltEn is: " << std::to_string(rateFiltEn) << std::endl;

    float omxFiltCut = crazyflie.getParam<float>(192);
    std::cout << "omxFiltCut is: " << omxFiltCut << std::endl;
    float omyFiltCut = crazyflie.getParam<float>(193);
    std::cout << "omyFiltCut is: " << omyFiltCut << std::endl;
    
    return 0;
}