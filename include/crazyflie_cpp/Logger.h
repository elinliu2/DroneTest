#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string>
#include <ostream>

class Logger {
private:
    std::ofstream file;

public:
    Logger(const std::string& filename) : file(filename){};

    // for strings int etc
    template<typename T>
    Logger& operator<<(const T& value) {
        file << value;
        return *this;
    }

    // for std::endl;
    Logger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        file << manip;
        return *this;
    }

    void log(std::string const& message) {
        file << message << std::endl;
    }

    void info(std::string const& message) {
        file << "INFO - " <<  message << std::endl;
    }

    void warn(std::string const& message) {
        file << "WARN - " << message << std::endl;
    }

    void error(std::string const& message) {
        file << "ERROR - " << message << std::endl;
    }


};
#endif