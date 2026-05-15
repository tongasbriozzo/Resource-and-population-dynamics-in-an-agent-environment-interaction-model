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
    /* INITIALIZATION */

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

    /* TERMALIZATION */
    ThermalizeSystem(0, newSystem, oldSystem, newBoard, oldBoard, grid, true);

    /* File */
    std::string nombreArchivo = std::string("m_G_") + 
                                   (N < 10 ? "00" : (N < 100 ? "0" : "")) + 
                                   std::to_string(N) + ".txt";
    ofstream movie_G(nombreArchivo);

    /* RECORDS */

    int N_rec = 102;//60000;
    if (newSystem.size()==0) N_rec = 1;

    int tTot = 0;
    int tTer = 0;
    int tAnt = 0;

    for (int t=0; t<N_rec; t++)
    {
        ThermalizeSystem(20 * tTer, newSystem, oldSystem, newBoard, oldBoard, grid, true);

        AverageInnerEnergy(newSystem, newBoard, averageInnerEnergy);
        /* INSTANT */
        float N_i = newSystem.size();
        float A_i = averageInnerEnergy[0];
  	    float E_i = averageInnerEnergy[1];
	   	float P_i = PolarOrder(newSystem);
        float O_i = NematicOrder(newSystem);
        float S_i = Entropy(newSystem);
        movie_G << tTot << " " << N_i << " " << A_i << " " << E_i 
                        << " " << P_i << " " << O_i << " " << S_i << std::endl;

        tTot = ceil(pow(10, 0.05 * t));
        tTer = tTot - tAnt;
        if (tTer < 1)
        {
            tTot = tAnt + 1;
            tTer = 1;
        }
        tAnt = tTot;
    }

    movie_G.close();

    std::cout << "\nEl programa finalizo con exito!\n";

  return 0;
}
