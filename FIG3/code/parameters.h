/******************************************************************************/
/* Parameters                                                                 */   
/* All model parameters are defined in this file                              */
/******************************************************************************/

#pragma once
#include <math.h>
#include <vector>
#define KIND float
using namespace std;

/* Physical and Mathematical constants */
const float pi = acos(-1.0);   // Pi
const float twoPi = 2.0*pi;    // 2 pi
const float gammaFriction = 3.92;//*activeVelocity; // Repulsive force between agents due to overlapping

/*** Setup ***/

/* Time */
const float timeStep           = 0.05;                              // Time step
const float termalizationTime  = 100000;                            // Termalization time
const int   termalizationSteps = int( termalizationTime/timeStep ); // Termalization steps
const float correlationTime    = 10;                                // Correlation time
const int   correlationSteps   = int( correlationTime/timeStep );   // Correlation steps

/* Board and Patches size */

const int escala     = 2; // Esta variable sirve para reescalear el sistema completo en dimension
const int proportion = 1; // Proporcion entre los lados del sistema

const float Ls = 150 * escala; // System side
const int   M  = 15  * escala; // Patch number per side

const float Lp  = Ls/M;            // Patch side
const float Lx  = Ls * proportion; // System side x
const float Ly  = Ls;              // System side y
const int   Mx  = M  * proportion; // Patch number in x
const int   My  = M;               // Patch number in y
const int   M2  = Mx*My;           // Total patch number

/* Agents size */
const float radius   = 1.0;        // Agent radius
const float diameter = 2.0*radius; // Agent diameter
const float rCut     = diameter;  // Cut-off radius for repulsive forces

/* Linked list grid */
const int lxGrid = floor(Lx/rCut); // Grids by side x
const int lyGrid = floor(Ly/rCut); // Grids by side y
const int nGrid  = lxGrid*lyGrid;  // Total grid number

const vector<vector<int>> ijGrid = {{0, 0}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

/* Attributes */
const bool birthsAndDeaths              = true;  // Activate births and deaths
const bool saturation                   = false;  // Activate intake saturation
const bool foodRotation                 = false;  // Activate chemotaxis