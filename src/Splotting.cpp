#include "Splotting.h"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

void splotTrajectory(SimResults simResults, Logger & log)
{
    
    // Use full path if PATH is not set
    const char* cmd = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent";
    FILE* gp = _popen(cmd, "w");
    if (!gp) {
        printf("Error: Could not open gnuplot!\n");
        return;
    }

    // Send gnuplot commands
    fprintf(gp, "set title 'Drone Trajectory'\n");
    fprintf(gp, "set xlabel 'x[m]'\n");
    fprintf(gp, "set ylabel 'y[m]'\n");
    fprintf(gp, "set zlabel 'z[m]'\n");
    fprintf(gp, "splot '-' with lines title 'trajectory'\n");

    // Send the data
    for (SystemState state : simResults.stateProgression) {
        log << "INFO - PLANT STATE:" << state.plant(x) << " "
             << state.plant(y) << " "
             << state.plant(z) << "\n";
        fprintf(gp, "%f %f %f\n", state.plant(x), state.plant(y), state.plant(z));
    }

    fprintf(gp, "e\n");  // 'e' ends the data section
    _pclose(gp);  // close gnuplot
    
}

void splotPlantState(SimResults simResults, Logger & log, int plantIndex)
{
    
    // Use full path if PATH is not set
    const char* cmd = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent";
    FILE* gp = _popen(cmd, "w");
    if (!gp) {
        printf("Error: Could not open gnuplot!\n");
        return;
    }

    // Send gnuplot commands
    fprintf(gp, "set title 'Plant State %d' \n", plantIndex);
    fprintf(gp, "set xlabel 'time'\n");
    fprintf(gp, "set ylabel 'plantSate'\n");
    fprintf(gp, "plot '-' with lines title 'trajectory'\n");

    // Send the data
    for (int i = 0; i < (int) simResults.time.size(); i++) {
        log << "INFO - PLANT STATE:" << simResults.time.at(i) << " "
             << simResults.stateProgression.at(i).plant(plantIndex) << "\n";
        fprintf(gp, "%f %f\n", simResults.time.at(i), simResults.stateProgression.at(i).plant(plantIndex));
    }

    fprintf(gp, "e\n");  // 'e' ends the data section
    _pclose(gp);  // close gnuplot
    
}



void testPlot()
{
    // Use full path if PATH is not set
    const char* cmd = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent";
    FILE* gp = _popen(cmd, "w");
    if (!gp) {
        printf("Error: Could not open gnuplot!\n");
        return;
    }

    // Send gnuplot commands
    fprintf(gp, "set title 'Sine Wave'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'sin(x)'\n");
    fprintf(gp, "plot '-' with lines title 'sin(x)'\n");

    // Send the data
    for (double x = 0; x <= 10; x += 0.1) {
        fprintf(gp, "%f %f\n", x, sin(x));
    }

    fprintf(gp, "e\n");  // 'e' ends the data section

    _pclose(gp);  // close gnuplot
}