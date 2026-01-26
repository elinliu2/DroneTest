#include "Crazyflie.h"
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
    std::string link_uri = "E7E7E7E7E7";
    TestLogger testLogger; 
    Crazyflie crazyflie(link_uri, testLogger);
    crazyflie.requestParamToc(true);
    // figure out id of one of the pid gains
    crazyflie.setParam(12, 1.0);
    crazyflie.requestParamToc(true);

    return 0;
}