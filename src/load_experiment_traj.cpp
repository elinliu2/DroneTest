#include "DroneTrajectory.h"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// ---------------------------------------------
// Data structures
// ---------------------------------------------

struct PositionData {
    double timestamp;
    double x;
    double y;
    double z;
    double xdot;
    double ydot;
    double zdot;
    double p;
    double q;
    double r;
};

struct AttitudeData {
    double timestamp;
    double phi;
    double theta;
    double psi;
};


struct PlantTrajectory {
    std::vector<PositionData> position;
    std::vector<AttitudeData> attitude;
};


PlantTrajectory g_activeTraj_experiment;


// ---------------------------------------------
// Load position CSV
// timestamp,x,y,z,xdot,ydot,zdot
// ---------------------------------------------

std::vector<PositionData> loadPositionTrajectory(
    const std::string& filepath)
{
    double firstTimestamp = 0;
    std::vector<PositionData> data;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open position CSV");
    }

    std::string line;
    bool header = true;

    while (std::getline(file, line)) {

        if (header) {
            header = false;
            continue;
        }

        std::stringstream ss(line);
        std::string field;
        std::vector<double> values;

        while (std::getline(ss, field, ',')) {
            values.push_back(std::stod(field));
        }

        if (values.size() != 10)
            continue;


        if (firstTimestamp == 0)
            firstTimestamp = values[0];


        PositionData d;

        d.timestamp = (values[0] - firstTimestamp) / 1e9;

        // assume mm -> m conversion
        d.x = values[1] / 1000.0;
        d.y = values[2] / 1000.0;
        d.z = values[3] / 1000.0;

        d.xdot = values[4] / 1000.0;
        d.ydot = values[5] / 1000.0;
        d.zdot = values[6] / 1000.0;

        d.p    = values[7] / 1000.0;
        d.q    = values[8] / 1000.0;
        d.r    = values[9] / 1000.0;

        data.push_back(d);
    }

    return data;
}


// ---------------------------------------------
// Load attitude CSV
// timestamp,roll,pitch,yaw,
// ---------------------------------------------

std::vector<AttitudeData> loadAttitudeTrajectory(
    const std::string& filepath)
{
    double firstTimestamp = 0;
    std::vector<AttitudeData> data;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open attitude CSV");
    }


    std::string line;
    bool header = true;

    while (std::getline(file, line)) {

        if (header) {
            header = false;
            continue;
        }


        std::stringstream ss(line);
        std::string field;
        std::vector<double> values;


        while (std::getline(ss, field, ',')) {
            values.push_back(std::stod(field));
        }


        if (values.size() != 4)
            continue;


        if (firstTimestamp == 0)
            firstTimestamp = values[0];


        AttitudeData d;

        d.timestamp = (values[0] - firstTimestamp) / 1e9;

        d.phi = values[1];
        d.theta = values[2];
        d.psi = values[3];

        data.push_back(d);
    }

    return data;
}

template <typename T, typename F>
double interpolateField(
    const std::vector<T>& data,
    double t,
    F getter)
{
    if(data.empty())
        throw std::runtime_error("Empty trajectory");


    if(t <= data.front().timestamp)
        return getter(data.front());


    if(t >= data.back().timestamp)
        return getter(data.back());


    auto it = std::lower_bound(
        data.begin(),
        data.end(),
        t,
        [](const T& d, double time)
        {
            return d.timestamp < time;
        });


    const T& after = *it;
    const T& before = *(it-1);


    double alpha =
        (t-before.timestamp) /
        (after.timestamp-before.timestamp);


    return getter(before)
        + alpha*(getter(after)-getter(before));
}

Eigen::Vector<double, NUM_PLANT_STATES>
interpolateState(double t)
{
    Eigen::Vector<double, NUM_PLANT_STATES> state;


    state(x) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.x;});


    state(y) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.y;});


    state(z) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.z;});


    state(xdot) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.xdot;});


    state(ydot) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.ydot;});


    state(zdot) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.zdot;});


    state(phi) =
        interpolateField(
            g_activeTraj_experiment.attitude,
            t,
            [](auto& d){return d.phi;});


    state(theta) =
        interpolateField(
            g_activeTraj_experiment.attitude,
            t,
            [](auto& d){return d.theta;});


    state(psi) =
        interpolateField(
            g_activeTraj_experiment.attitude,
            t,
            [](auto& d){return d.psi;});


    state(p) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.p;});


    state(q) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.q;});


    state(r) =
        interpolateField(
            g_activeTraj_experiment.position,
            t,
            [](auto& d){return d.r;});


    return state;
}

void setActiveTrajectory(
    const std::string& positionFile,
    const std::string& attitudeFile)
{
    g_activeTraj_experiment.position =
        loadPositionTrajectory(positionFile);

    g_activeTraj_experiment.attitude =
        loadAttitudeTrajectory(attitudeFile);
}