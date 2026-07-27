#include "DroneTrajectory.h"
#include "Logger.h"
#include <iostream>
#include <cmath>

#include <sstream>
#include <algorithm>
#include <stdexcept>

// ---------------------------------------------
// Data structure for a single CSV row
// ---------------------------------------------
struct DataPoint {
    double timestamp;
    double x;
    double y;
    double z;
    double yaw;
};

// Which field to interpolate
enum class Field { X, Y, Z, YAW };

std::vector<DataPoint> g_activeTraj;

void setActiveTrajectory(std::vector<DataPoint> traj) {
    g_activeTraj = std::move(traj);
}

static double getField(const DataPoint& dp, Field field) {
    switch (field) {
        case Field::X:   return dp.x;
        case Field::Y:   return dp.y;
        case Field::Z:   return dp.z;
        case Field::YAW: return dp.yaw;
    }
    return 0.0;
}

static bool parseLine(const std::string& line, DataPoint& dp, double & firstTimestep) {
    std::stringstream ss(line);
    std::string field;
    std::vector<double> values;

    try {
        while (std::getline(ss, field, ',')) {
            // trim whitespace
            size_t start = field.find_first_not_of(" \t\r\n");
            size_t end = field.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) {
                values.push_back(0.0);
                continue;
            }
            field = field.substr(start, end - start + 1);
            values.push_back(std::stod(field));
        }
    } catch (const std::exception&) {
        return false;
    }

    if (values.size() != 5) return false;

    if (firstTimestep == 0){
        firstTimestep = values[0];
    }
    dp.timestamp = (values[0]-firstTimestep)/1e9;
    dp.x         = values[1]/1000.0;
    dp.y         = values[2]/1000.0;
    dp.z         = values[3]/1000.0;
    dp.yaw       = values[4];
    return true;
}

// ---------------------------------------------
// Load CSV file into memory. Expects header row: timestamp,x,y,z,yaw
// Returns the data sorted by timestamp. Empty vector on failure.
// ---------------------------------------------
std::vector<DataPoint> loadTrajectory(const std::string& filepath) {
    std::vector<DataPoint> data;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filepath << "\n";
        return data;
    }

    std::string line;
    bool first_line = true;
    double firstTimestep = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // Skip header (assumes first line is a header if it can't be parsed as a number)
        if (first_line) {
            first_line = false;
            continue;
        }

        DataPoint dp;
        if (parseLine(line, dp, firstTimestep)) {
            data.push_back(dp);
        } else {
            std::cerr << "Warning: skipping malformed line: " << line << "\n";
        }
    }

    return data;
}

// ---------------------------------------------
// Interpolate a given field at query_time.
// data must already be sorted by timestamp (loadTrajectory guarantees this).
// Throws std::out_of_range if data is empty.
// Clamps to endpoints if query_time is outside the recorded range.
// ---------------------------------------------
double interpolate(const std::vector<DataPoint>& data, double query_time, Field field) {
    if (data.empty()) {
        throw std::out_of_range("Trajectory has no data loaded.");
    }
    if (data.size() == 1) {
        return getField(data[0], field);
    }

    // Clamp outside range
    if (query_time <= data.front().timestamp) {
        return getField(data.front(), field);
    }
    if (query_time >= data.back().timestamp) {
        return getField(data.back(), field);
    }

    // Binary search for first element with timestamp >= query_time
    auto it = std::lower_bound(
        data.begin(), data.end(), query_time,
        [](const DataPoint& dp, double t) { return dp.timestamp < t; });

    // it points to the point at/after query_time; the previous point brackets it below
    const DataPoint& after = *it;
    const DataPoint& before = *(it - 1);

    double t0 = before.timestamp;
    double t1 = after.timestamp;
    double v0 = getField(before, field);
    double v1 = getField(after, field);

    if (t1 == t0) return v0; // duplicate timestamps guard

    double alpha = (query_time - t0) / (t1 - t0);
    return v0 + alpha * (v1 - v0);
}

// Convenience wrappers
double interpolateX(double t)   { return interpolate(g_activeTraj, t, Field::X); }
double interpolateY(double t)   { return interpolate(g_activeTraj, t, Field::Y); }
double interpolateZ(double t)   { return interpolate(g_activeTraj, t, Field::Z); }
double interpolateYaw(double t) { return interpolate(g_activeTraj, t, Field::YAW); }

