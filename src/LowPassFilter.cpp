#include "LowPassFilter.h"

#define _USE_MATH_DEFINES // Must be defined before including cmath or math.h
#include <cmath>

// https://github.com/bitcraze/crazyflie-firmware/blob/master/src/utils/src/filter.c#L71
LowPassFilter::LowPassFilter(double sampleRate, double cutoffFreq)
{
    double fr = sampleRate/cutoffFreq;
    double ohm = tan(M_PI/fr);
    double c = 1.0+2.0*cos(M_PI/4.0)*ohm+ohm*ohm;
    b0 = ohm*ohm/c;
    b1 = 2.0*b0;
    b2 = b0;
    a1 = 2.0*(ohm*ohm-1.0f)/c;
    a2 = (1.0-2.0*cos(M_PI/4.0)*ohm+ohm*ohm)/c;
    delay_element_1 = 0.0;
    delay_element_2 = 0.0;
}

double LowPassFilter::lpf2pApply(double sample)
{
    float delay_element_0 = sample - delay_element_1 * a1 - delay_element_2 * a2;
    if (!std::isfinite(delay_element_0)) {
        // don't allow bad values to propigate via the filter
        delay_element_0 = sample;
    }

    float output = delay_element_0 * b0 + delay_element_1 * b1 + delay_element_2 * b2;

    delay_element_2 = delay_element_1;
    delay_element_1 = delay_element_0;
    return output;
}
