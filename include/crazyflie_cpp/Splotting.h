#include "DroneTrajectory.h"

#include "Logger.h"

void splotTrajectory(SimResults simResults, Logger & log, std::string const& plotTitle = "");
void splotPlantState(SimResults simResults, Logger & log, int plantIndex, std::string const& plotTitle = "");
void testPlot();