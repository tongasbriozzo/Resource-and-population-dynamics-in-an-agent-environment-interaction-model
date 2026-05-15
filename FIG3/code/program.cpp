/******************************************************************************/
/* Programa Principal                                                         */   
/* Este programa reporta los parametros de orden del sistema una vez alcanzad */
/* o el equilibrio. Se Inicializa y deja termalizar al sistema. Existe un par */
/* ametro de control, en este caso el foodRotationRate, que se varia para obt */
/* tener la dependecia de los parametros de orden con este. El programa esta  */
/* paralelizado                                                               */
/******************************************************************************/

#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <vector>
#include <omp.h>
#include "functions.h"

int main()
{
    /* Numero de simulaciones a realizar (cantidad de valores para el parametro de control) */
    int nD = 5;
    int nN = 45;

    for (int iD=0; iD<nD; iD++)
    {

        /* PARALELIZATION */
        omp_set_num_threads(21);      // Seteamos el numero de nucleos a usar
        #pragma omp parallel for     // Paralelizamos
        for (int iN=0; iN<nN; iN++)
        {
            /* INITIALIZATION */

            /* Parameters */
            float coefDif          = pow(10.0, iD + 1.0);         // Coeficiente de difusion
            float angularDiffusion = 0.01*pow(2.0*coefDif, -1.0/3.0); // Coeficiente de difusion angular
            float activeVelocity   = 0.10*pow(2.0*coefDif,  1.0/3.0); // Velocidad activa

            float foodRotationRate = 0.00;                            // Fuerza de chemotaxis, Parametro de control
            
            float starvationEnergy   = 0.1;
            float birthEnergy        = 0.4;
            float reproductiveEnergy = 0.9;

            float intakeSlope    = 0.0100;
            float intakeCapacity = 0.0100;
            float metabolicRate  = 0.0050;
            float kineticRate    = 0.0005;

            float growthRate    = 0.015;
            float patchCapacity = 1.000;

            /* Constants */            

            float nStarvation   = growthRate/intakeSlope * (1.0 - (metabolicRate*starvationEnergy  )/(intakeSlope*patchCapacity));
            float nReproduction = growthRate/intakeSlope * (1.0 - (metabolicRate*reproductiveEnergy)/(intakeSlope*patchCapacity));

            int n0 = int( float(M2) * 2.25*nStarvation * float(iN)/float(nN) ); // Agents Initial Number

            if ( ! (0<n0) ) n0 = 1;
            int nRoot = ceil(sqrt(n0)); // Para acomodar a los agentes segun una red uadrada (EVITA ERROR NUMERICO POR SUPERPOSICION)
            
            float f0 = patchCapacity*(1.0 - intakeSlope/growthRate*n0/M2); // Equilibrium energy for the patch in a homogeneous set
            float e0 = intakeSlope/metabolicRate * f0;      
            if ( ! (               0<f0 && f0<patchCapacity     ) ) f0 = 0.5*patchCapacity;
            if ( ! (starvationEnergy<e0 && e0<reproductiveEnergy) ) e0 = 0.5*(starvationEnergy+reproductiveEnergy);
   
            /* Report */
            std::cout << "iN: " << iN << ", " << std::endl;
            std::cout << "Difussion Coefficient" << " " << coefDif << std::endl;
            std::cout << "activeVelocity" << " " << activeVelocity << " " << "angularDiffusion" << " " << angularDiffusion << " " << "kineticRate" << " " << kineticRate << std::endl;

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
            /* NOTA: Se puede variar la capacidad de los parches en funcion de la posicion por ejemplo si se lo desea */

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

            int   nRec = 60; // Numero de registros
            float tOld = 0.0;
            float tNew = 0.0;

            std::ostringstream oss;
            oss << std::setw(2) << std::setfill('0') << iN;

            ofstream records("population_" + std::to_string(iD) + "_" + oss.str() + ".txt");

            /* Realizamos las mediciones */
            for (int iRec=0; iRec<nRec; iRec++)
            {
                tOld = tNew;
                if (iRec<=10) tNew = iRec;
                if (iRec> 10) tNew = floor(pow(10, 0.1*iRec));

                // We let pass correlationSteps time steps to avoid correlation between measurements
                ThermalizeSystem(int((tNew-tOld)/timeStep), newSystem, oldSystem, newBoard, oldBoard, grid, false);

                // Obtenemos las energias medias
                AverageInnerEnergy(newSystem, newBoard, averageInnerEnergy);

                /* INSTANT */
                float A_i = averageInnerEnergy[0];
  	            float E_i = averageInnerEnergy[1];
                float S_i = Entropy(newSystem);
                float N_i = float(newSystem.size());
	   	        float P_i = PolarOrder(newSystem);
                float Q_i = NematicOrder(newSystem);
        
                records << tNew << " " << N_i << " " << A_i << " " << E_i 
                                << " " << P_i << " " << Q_i << " " << S_i << std::endl;
        
            }

            // close the document
            records.close();
        }

    }

    // motivational message IMPORTANT!!!
    std::cout << "\nProgram finished\n";

    return 0;
}