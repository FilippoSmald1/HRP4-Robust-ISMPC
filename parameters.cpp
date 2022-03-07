#pragma once

#include <math.h>

// Definable parameters
// ********************

// Times
const double mpcTimeStep = 0.01; //0.05;
const double controlTimeStep = 0.01;
const double singleSupportDuration = 0.6; //0.3
const double doubleSupportDuration = 0.2; // 0.2
const double predictionTime = 0.7;  //0.7
const double previewTime = 2.5;

// Disturbance
const double midrange_x = 0.0;
const double midrange_y = 0.0;
const double dist_range_x = 0.5;
const double dist_range_y = 0.5;

// Walk parameters
const double stepHeight = 0.025; //0.03
const double comTargetHeight = 0.71;
const double kSwingFoot = 0.05; 

// Constraints
const double thetaMax = 0.30;
const double footConstraintSquareWidth = 0.08;
const double wkx = 0.5;
const double ell = 0.2;
const double wky = 0.16;
const double deltaYIn = 0.14;
const double deltaYOut = 0.3;

// Cost function weights for horizontal QP
const double qZd = 1.0;
const double qZ = 0.0;
const double qF = 100000.0;
// Cost function weights for vertical QP
const double q_force = 1.0;
const double q_position = 1000000000000.0;

// Kinematic control
const double IKerrorGain = 1.0; 

// Used in the code
// ****************
const double mass_hrp4 = 50.0; 
const double g = 9.81; 
const double eta = sqrt(g/comTargetHeight);
const int N = round(predictionTime/mpcTimeStep);
const int S = round(singleSupportDuration/mpcTimeStep);
const int F = round(doubleSupportDuration/mpcTimeStep);
const int Prev = round((previewTime-predictionTime)/mpcTimeStep);
const int M = 2; 


