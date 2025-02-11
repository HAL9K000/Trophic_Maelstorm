#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

// Store wind-speed and flipping duration of wind direction in an inline global variable.
inline std::pair<double, double> wind_speed = {0, 0};
// Stores wind-speed and direction in x and y components.
inline int wind_duration = 0;
//wind_params[1] = wind_duration (in hrs) 
// NOTE: [0 indicates no-flipping of wind direction].