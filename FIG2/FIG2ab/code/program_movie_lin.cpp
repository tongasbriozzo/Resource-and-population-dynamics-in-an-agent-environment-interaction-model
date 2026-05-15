/******************************************************************************/
/* Programa para Peliculas                                                    */   
/* Este programa reporta las posiciones y energias internas de agentes y parc */
/* hes para realizar una pelicula de la evolucion temporal                    */
/******************************************************************************/

#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // Para usar std::setw y std::setfill
#include <vector>
#include <omp.h>
#include "functions.h"

int main()
{  
    /* INICIALIZATION*/

    /* Parameters */
    float foodRotationRate = 0.1;  // Fuerza de chemotaxis, Parametro de control
    float angularDiffusion = 0.01; // Coeficiente de difusion angular
    float activeVelocity   = 0.5;  // Velocidad activa
            
    float starvationEnergy   = 0.1;
    float birthEnergy        = 0.4;
    float reproductiveEnergy = 0.9;

    float intakeSlope    = 0.0150;
    float intakeCapacity = 0.0150;
    float metabolicRate  = 0.0050;
    float kineticRate    = 0.0005;

    float growthRate    = 0.015;
    float patchCapacity = 1.000;

    /* Constants */            
    float coefDif = activeVelocity*activeVelocity / (2.0 * angularDiffusion); // Coeficiente de difusion

    float nStarvation   = growthRate/intakeSlope * (1.0 - (metabolicRate*starvationEnergy  )/(intakeSlope*patchCapacity));
    float nReproduction = growthRate/intakeSlope * (1.0 - (metabolicRate*reproductiveEnergy)/(intakeSlope*patchCapacity));

    int n0    = 64;// int( M2*(nStarvation + nReproduction)*0.5 ); // Agents Initial Number
    if ( ! (0<n0) ) n0 = 1;
    int nRoot = ceil(sqrt(n0)); // Para acomodar a los agentes segun una red uadrada (ECITA ERROR NUMERICO POR SUPERPOSICION)
            
    float f0 = patchCapacity*(1.0 - intakeSlope/growthRate*n0/M2); // Equilibrium energy for the patch in a homogeneous set
    float e0 = intakeSlope/metabolicRate * f0;      
    if ( ! (               0<f0 && f0<patchCapacity     ) ) f0 = 0.5*patchCapacity;
    if ( ! (starvationEnergy<e0 && e0<reproductiveEnergy) ) e0 = 0.5*(starvationEnergy+reproductiveEnergy);
   
    /* Report */
    std::cout << "Difussion Coefficient" << " " << coefDif << std::endl;
    std::cout << "activeVelocity" << " " << activeVelocity << " " << "angularDiffusion" << " " << angularDiffusion << " " << "foodRotationRate" << " " << foodRotationRate << std::endl;

    /* Inicializamos los agentes */
    vector<Agent> oldSystem;
    vector<Agent> newSystem;
    oldSystem.resize(n0);
    newSystem.resize(n0);
    for (int i=0; i<n0; i++)
    {
        newSystem[i].v0 = activeVelocity;
        newSystem[i].Dt = angularDiffusion;
        newSystem[i].g  = foodRotationRate;
        newSystem[i].es = starvationEnergy;
        newSystem[i].eb = birthEnergy;
        newSystem[i].er = reproductiveEnergy;
        newSystem[i].Is = intakeSlope;
        newSystem[i].Ic = intakeCapacity;
        newSystem[i].m  = metabolicRate;
        newSystem[i].k  = kineticRate;

        newSystem[i].x = (Lx/nRoot)*(float(i%nRoot)+0.5);
        newSystem[i].y = (Ly/nRoot)*(floor(i/nRoot)+0.5);
        if (Lx/nRoot > radius) newSystem[i].x += (Lx/nRoot-radius)*0.5*(dis(gen) - 0.5);
        if (Ly/nRoot > radius) newSystem[i].y += (Lx/nRoot-radius)*0.5*(dis(gen) - 0.5);
        newSystem[i].PeriodicBoundaryConditions();

        float ei = e0 + (e0-starvationEnergy)*(reproductiveEnergy-e0)*disGauss(gen);
        if ( ! (starvationEnergy<ei && ei<reproductiveEnergy) ) ei = e0;
        newSystem[i].e = ei;
    }
    oldSystem = newSystem;
    /* NOTA: Para evolucion estos valores no tienen por que ser identicos para todos los agentes */

    /* Inicializamos los parches */
    vector<Patch> oldBoard;
    vector<Patch> newBoard;
    oldBoard.resize(M2);
    newBoard.resize(M2);
    for (int i=0; i<M2; i++)
    {
        newBoard[i].r = growthRate;
        newBoard[i].c = patchCapacity;

        float fi = f0 + (f0-0)*(patchCapacity-f0)*disGauss(gen);
        if ( ! (0<fi && fi<patchCapacity) ) fi = f0;
        newBoard[i].f = fi;
    }
    oldBoard = newBoard;

    /* Inicializamos la grilla */
    vector<vector<set<size_t>>> grid;
    grid.resize(lxGrid);
    for (int i=0; i<lxGrid; i++)
    {
        grid[i].resize(lyGrid);
    }

    /* Inicializamos el vector de energias medias */
    vector<float> averageInnerEnergy;
    averageInnerEnergy.resize(2);

    /* RECORDS */

    int N_rec = 1440;
    if (newSystem.size()==0) N_rec = 1;

    ofstream movie_A("m_A.txt");
    ofstream movie_P("m_P.txt");
    ofstream movie_G("m_G.txt");

    for (int t=0; t<N_rec; t++)
    {
        ThermalizeSystem(correlationSteps*0.5, newSystem, oldSystem, newBoard, oldBoard, grid, false);

        AverageInnerEnergy(newSystem, newBoard, averageInnerEnergy);
        /* INSTANT */
        float N_i = newSystem.size();
        float A_i = averageInnerEnergy[0];
  	    float E_i = averageInnerEnergy[1];
	   	float P_i = PolarOrder(newSystem);
        float O_i = NematicOrder(newSystem);
        float S_i = Entropy(newSystem);
        
        movie_G << t*correlationTime << " " << N_i << " " << A_i << " " << E_i 
                                     << " " << P_i << " " << O_i << " " << S_i << std::endl;

        movie_A << std::endl;
        movie_P << std::endl;

        for (int iAgent=0; iAgent<newSystem.size(); iAgent++)
        {
            movie_A << iAgent << " " << newSystem[iAgent].x  << " " << newSystem[iAgent].y
                              << " " << newSystem[iAgent].vx << " " << newSystem[iAgent].vy 
                              << " " << newSystem[iAgent].e  << std::endl;
        }

        for (int iPatch=0; iPatch<M2; iPatch++)
        {
            int xPatch, yPatch;
            PatchIndexs(iPatch, xPatch, yPatch);
            movie_P << iPatch << " " << xPatch << " " << yPatch << " " << newBoard[iPatch].f/newBoard[iPatch].c << std::endl;
        }
        

    }

    movie_A.close();
    movie_P.close();
    movie_G.close();

    std::cout << "\nEl programa finalizo con exito!\n";

  return 0;
}
