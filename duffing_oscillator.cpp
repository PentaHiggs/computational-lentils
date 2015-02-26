/* Basic program to explore the chaotic properties of the Duffing Oscillator.
 * Uses cern's ROOT library/software in order to draw figures/graphs.
 * Program CANNOT compile on it's own without modification.  However, all code used is 
 * valid C++ code, the file is just not an independent C program.
 */

// C library includes
#include<functional>
#include<cmath>
#include<iostream>
#include<vector>

// ROOT includes
#include"TRoot.h"
#include"TF1.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TMultiGraph.h"

// My integration routine
#include"runge_kutta1.cpp"


