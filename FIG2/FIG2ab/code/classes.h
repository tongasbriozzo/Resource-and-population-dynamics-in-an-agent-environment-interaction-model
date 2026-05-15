/******************************************************************************/
/* Clases                                                                     */   
/* All model clases, Agents and Patches, are defined in this file             */
/******************************************************************************/

#include "parameters.h"    //Módulo con parametros para la simulación.
#include <bits/stdc++.h>
#include <random>
using namespace std;

/******************************************************************************/
/* Random Numbers                                                             */
/******************************************************************************/

/* Homogeneous distribution in the interval 0 1 */
mt19937::result_type seed = time(0);
mt19937 gen(seed);                             //Standard mersenne_twister_engine seeded time(0)
uniform_real_distribution<float> dis(0., 1.); // dis(gen), número aleatorio real homogeneamente distribuido entre 0 y 1

/* Gaussian distribution of mean 0 and variance 1 */
float media = 0.0;
float desviacion_estandar = 1.0;
std::normal_distribution<float> disGauss(media, desviacion_estandar); // Objeto que genera números aleatorios distribuidos normalmente

/******************************************************************************/
/* Clase Agente                                                               */
/******************************************************************************/
class Agent
{
public:
    Agent();
    float x, y;                         // x position, y position
    float a, Dt, g;                     // angle, angular diffusion, angular alignment
    float v0, vx, vy;                   // active velocity, x velocity, y velocity
    float e, es, eb, er;                // inner energy, starvation energy, birth energy, reproduction energy
    float Is, Ic, m, k;                 // intake slope, intake capacity, metabolic rate, kinetic rate
    int   GetPatchIndex();              // indice del parche donde se encuentra el agente
    int   GetGridIndex();               // indice de la grilla donde se encuentra el agente
    void  RandomTumble();               // aplicar ruido gaussiano al angulo a del agente
    void  PeriodicBoundaryConditions(); // aplicar condiciones periodicas de contorno a la posicion x y del agente
};

/* Inicializacion del Agente */
Agent::Agent()
{
    /* position */
    x = dis(gen)*Lx;           
    y = dis(gen)*Ly;            

    /* angles */
    a  = (dis(gen) - 0.5)*twoPi; 
    Dt = 0;
    g  = 0;

    /* valocity */
    v0 = 1.0;
    vx = v0 * cos(a);
    vy = v0 * sin(a);

    /* energy */
    e  = 0.5;
    es = 0.1;
    eb = 0.4;
    er = 0.9;

    /* metabolism */
    Is = 0.0150;
    Ic = 0.0150;
    m  = 0.0050;
    k =  0.0005;
}

/* Indice del parche donde se encuentra el agente */
int Agent::GetPatchIndex()
{
    int index;
    index = floor(x/Lp) + Mx*floor(y/Lp);
    if ( ! ( 0<=index && index<M2 ) ) // WARNING
    {
        std::cout << std::endl << "WARNING! Patch index: " << index            << std::endl;
        std::cout <<              "Patches number:       " << M2               << std::endl;
        std::cout <<              "x , y :               " << x  << ", " << y  << std::endl;
        std::cout <<              "Lx, Ly:               " << Lx << ", " << Ly << std::endl;
        std::cout <<              "Lp:                   " << Lp               << std::endl;
    }
    return index;
}

/* Indice de la grilla donde se encuentra el agente */
int Agent::GetGridIndex()
{
    int index;
    index = floor(x/Lx * lxGrid) + lxGrid*floor(y/Ly * lyGrid);
    if ( ! ( 0<=index && index<nGrid ) ) // WARNING
    {
        std::cout << std::endl << "WARNING! Grid index: " << index                    << std::endl;
        std::cout <<              "Grids number:        " << nGrid                    << std::endl;
        std::cout <<              "x     , y     :      " << x      << ", " << y      << std::endl;
        std::cout <<              "Lx    , Ly    :      " << Lx     << ", " << Ly     << std::endl;
        std::cout <<              "lxGrid, lyGrid:      " << lxGrid << ", " << lyGrid << std::endl;
    }
    return index;
}

/* Aplicar difusion angular */
void Agent::RandomTumble()
{
    a += sqrt(2.0 * Dt * timeStep) * disGauss(gen);
    while (a >   pi) a -= twoPi;
    while (a <= -pi) a += twoPi;
}

/* Aplicar condiciones periodicas de contorno */
void Agent::PeriodicBoundaryConditions()
{
    if ( ! ( -Lx<=x && x<=2*Lx ) ) std::cout << std::endl << "WARNING! Position x: " << x << std::endl; // WARNING
    if ( ! ( -Ly<=y && y<=2*Ly ) ) std::cout << std::endl << "WARNING! Position y: " << y << std::endl; // WARNING
    while (x < 0) x = x+Lx;
    while (Lx<=x) x = x-Lx;
    while (y < 0) y = y+Ly;
    while (Ly<=y) y = y-Ly;
}

/******************************************************************************/
/* Clase Patch                                                                */
/******************************************************************************/
class Patch
{
public:
    Patch();
    float f, r, c; // food, growth rate, patch capacity
};

/* Inicializacion del Parche */
Patch::Patch()
{
    f = 0.5;
    r = 0.015;
    c = 1.0;
}