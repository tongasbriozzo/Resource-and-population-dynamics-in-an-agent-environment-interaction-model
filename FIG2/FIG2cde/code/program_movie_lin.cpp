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
    float foodRotationRate = 0.00;
    float angularDiffusion = 0.10; 
    float activeVelocity   = 0.50;
    float coefDif          = activeVelocity*activeVelocity / (2.0 * angularDiffusion);
    int   nRoot            = ceil(sqrt(N));

    /* Report */
    std::cout << "Difussion Coefficient" << " " << coefDif << std::endl;
    std::cout << "activeVelocity" << " " << activeVelocity << " " << "angularDiffusion" << " " << angularDiffusion << std::endl;

    /* Inicializamos los agentes */
    vector<Agent> oldSystem;
    vector<Agent> newSystem;
    oldSystem.resize(N);
    newSystem.resize(N);
    for (int i=0; i<N; i++)
    {
        newSystem[i].SetActiveVelocity(activeVelocity);
        newSystem[i].SetAngularDiffusion(angularDiffusion);
        newSystem[i].SetFoodRotationRate(foodRotationRate);
        newSystem[i].x = radius + (Lx/nRoot) * float(i%nRoot) ;
        newSystem[i].y = radius + (Ly/nRoot) * floor(i/nRoot) ;
    }
    oldSystem = newSystem;

    /* Inicializamos los parches */
    vector<Patch> oldBoard;
    vector<Patch> newBoard;
    oldBoard.resize(M2);
    newBoard.resize(M2);
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

    int N_rec = 1200;
    if (newSystem.size()==0) N_rec = 1;

    ofstream movie_A("m_A.txt");
    ofstream movie_P("m_P.txt");
    ofstream movie_G("m_G.txt");

    for (int t=0; t<N_rec; t++)
    {
        ThermalizeSystem(1000, newSystem, oldSystem, newBoard, oldBoard, grid, false);

        AverageInnerEnergy(newSystem, newBoard, averageInnerEnergy);
        /* INSTANT */
        float N_i = newSystem.size();
        float A_i = averageInnerEnergy[0];
  	    float E_i = averageInnerEnergy[1];
	   	float P_i = PolarOrder(newSystem);
        float O_i = NematicOrder(newSystem);
        float S_i = Entropy(newSystem);
        movie_G << t*timeStep*1000 << " " << N_i << " " << A_i << " " << E_i 
                                   << " " << P_i << " " << O_i << " " << S_i << std::endl;

        movie_A << std::endl;
        movie_P << std::endl;

        for (int iAgent=0; iAgent<newSystem.size(); iAgent++)
        {
            movie_A << iAgent << " " << newSystem[iAgent].x  << " " << newSystem[iAgent].y
                              << " " << newSystem[iAgent].vx << " " << newSystem[iAgent].vy 
                              << " " << newSystem[iAgent].GetInnerEnergy() << std::endl;
        }

        for (int iPatch=0; iPatch<M2; iPatch++)
        {
            int xPatch, yPatch;
            PatchIndexs(iPatch, xPatch, yPatch);
            movie_P << iPatch << " " << xPatch << " " << yPatch << " " << newBoard[iPatch].GetInnerEnergy()/newBoard[iPatch].GetCapacity() << std::endl;
        }
        

    }

    movie_A.close();
    movie_P.close();
    movie_G.close();

    std::cout << "\nEl programa finalizo con exito!\n";

  return 0;
}
