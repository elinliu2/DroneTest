#include <crazyflie_cpp/Crazyflie.h>
#include <iostream>

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
    float gain = crazyflie.getParam<float>(164);
    std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;
    crazyflie.setParam(164, 3.5f);
    gain = crazyflie.getParam<float>(164);
    std::cout << "pid attitude roll_kp gain is: " << gain << std::endl;
    
    return 0;
}