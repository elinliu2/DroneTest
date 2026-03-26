class LowPassFilter
{
    double a1;
    double a2;
    double b0;
    double b1;
    double b2;
    double delay_element_1;
    double delay_element_2;

    public: 
        LowPassFilter(double sampleRate = 1, double cutoffFreq = 1);
        double lpf2pApply(double sample);
};