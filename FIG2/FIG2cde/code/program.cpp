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
    int nv0 = 21;
    int nDt = 21;

    /* PARALELIZATION */
    //  int numThreads = omp_get_num_procs(); // Numero de nucles
    //  if (numThreads>3) numThreads -= 2;    // No usamos todos los nucleos para no matar la CPU
    omp_set_num_threads(21);      // Seteamos el numero de nucleos a usar
    #pragma omp parallel for     // Paralelizamos
    for (int iv0=0; iv0<nv0; iv0++)
    {
        for (int iDt=0; iDt<nDt; iDt++)
        {
            /* INITIALIZATION */

            /* Parameters */
            float foodRotationRate = 0.00;                          // Fuerza de chemotaxis, Parametro de control
            float angularDiffusion = pow(10, float(iDt)*0.2 - 2.0); // Coeficiente de difusion angular
            float activeVelocity   = pow(10, float(iv0)*0.2 - 3.0); // Velocidad activa
            
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
            float coefDif = activeVelocity*activeVelocity / (2.0 * angularDiffusion); // Coeficiente de difusion

            float nStarvation   = growthRate/intakeSlope * (1.0 - (metabolicRate*starvationEnergy  )/(intakeSlope*patchCapacity));
            float nReproduction = growthRate/intakeSlope * (1.0 - (metabolicRate*reproductiveEnergy)/(intakeSlope*patchCapacity));

            int n0    = int( M2*(nStarvation + nReproduction)*0.5 ); // Agents Initial Number
            if ( ! (0<n0) ) n0 = 1;
            int nRoot = ceil(sqrt(n0)); // Para acomodar a los agentes segun una red uadrada (EVITA ERROR NUMERICO POR SUPERPOSICION)
            
            float f0 = patchCapacity*(1.0 - intakeSlope/growthRate*n0/M2); // Equilibrium energy for the patch in a homogeneous set
            float e0 = intakeSlope/metabolicRate * f0;      
            if ( ! (               0<f0 && f0<patchCapacity     ) ) f0 = 0.5*patchCapacity;
            if ( ! (starvationEnergy<e0 && e0<reproductiveEnergy) ) e0 = 0.5*(starvationEnergy+reproductiveEnergy);
   
            /* Report */
            std::cout << "iv0: " << iv0 << ", " << "iDt: " << iDt << std::endl;
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

                //newSystem[i].x = (Lx/nRoot)*(float(i%nRoot)+0.5);
                //newSystem[i].y = (Ly/nRoot)*(floor(i/nRoot)+0.5);
                //if (Lx/nRoot > radius) newSystem[i].x += (Lx/nRoot-radius)*0.5*(dis(gen) - 0.5);
                //if (Ly/nRoot > radius) newSystem[i].y += (Lx/nRoot-radius)*0.5*(dis(gen) - 0.5);
                //newSystem[i].PeriodicBoundaryConditions();

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

            /* TERMALIZATION */

            /* Termalizamos el sistema hasta alcanzar el equilibrio. Reportamos parametros de interes al finalizar */
            ThermalizeSystem(termalizationSteps, newSystem, oldSystem, newBoard, oldBoard, grid, true);

            /* RECORDS */    
            float A1m = 0.0, A2m = 0.0, A_i = 0.0; // Energia media de los agentes
            float E1m = 0.0, E2m = 0.0, E_i = 0.0; // Energia media de los parches
            float S1m = 0.0, S2m = 0.0, S_i = 0.0; // Entropia
            float N1m = 0.0, N2m = 0.0, N_i = 0.0; // Numero de agentes
            float P1m = 0.0, P2m = 0.0, P_i = 0.0; // Polarizacion
            float Q1m = 0.0, Q2m = 0.0, Q_i = 0.0; // Orden nematico

            int nRec = int(10000); // Numero de registros

            /* Realizamos las mediciones */
            for (int t=0; t<nRec; t++)
            {

                // We let pass correlationSteps time steps to avoid correlation between measurements
                ThermalizeSystem(correlationSteps, newSystem, oldSystem, newBoard, oldBoard, grid, false);

                // Obtenemos las energias medias
                AverageInnerEnergy(newSystem, newBoard, averageInnerEnergy);

                /* INSTANT */
                A_i = averageInnerEnergy[0];
  	            E_i = averageInnerEnergy[1];
                S_i = Entropy(newSystem);
                N_i = float(newSystem.size());
	   	        P_i = PolarOrder(newSystem);
                Q_i = NematicOrder(newSystem);

                /* MEAN */
                A1m += A_i;
                E1m += E_i;
                S1m += S_i;
                N1m += N_i;
                P1m += P_i;
                Q1m += Q_i;

                /* VARIANCE */
                A2m += (A_i*A_i);
                E2m += (E_i*E_i);
                S2m += (S_i*S_i);
                N2m += (N_i*N_i);
                P2m += (P_i*P_i);
                Q2m += (Q_i*Q_i);
        
            }

            /* NORMALIZATION */
            A1m /= float(nRec);
            A2m /= float(nRec);
            E1m /= float(nRec);
            E2m /= float(nRec);
            S1m /= float(nRec);
            S2m /= float(nRec);
            N1m /= float(nRec);
            N2m /= float(nRec);
            P1m /= float(nRec);
            P2m /= float(nRec);
            Q1m /= float(nRec);
            Q2m /= float(nRec);

            /* TRANSITORY REPORT */

            // name of the document
            std::string iSimState = std::string("iSim_") 
                                  + (iv0 < 10 ? "00" : (iv0 < 100 ? "0" : "")) + std::to_string(iv0)
                                  + "_"
                                  + (iDt < 10 ? "00" : (iDt < 100 ? "0" : "")) + std::to_string(iDt)
                                  + ".txt";                                                            
            
            // open the document
            ofstream ISimState(iSimState);     

            // save the data        
            ISimState << coefDif << ", "
                      << activeVelocity << ", "
                      << angularDiffusion << ", "
                      << foodRotationRate << ", "
                      << A1m << ", " << A2m << ", "
                      << E1m << ", " << E2m << ", "
                      << S1m << ", " << S2m << ", "
                      << N1m << ", " << N2m << ", "
                      << P1m << ", " << P2m << ", "
                      << Q1m << ", " << Q2m << ", "
                      << std::endl;

            // close the document
            ISimState.close();

        }
    }

    /* FINAL REPORT */  

    // name of the document
    ofstream FinalState("FinalState.txt");  

    // for each of the transitory documents...
    for (int iv0=0; iv0<nv0; iv0++)
    {
        for (int iDt=0; iDt<nDt; iDt++)
        {
            // name of the transitory report
            std::string iSimState = std::string("iSim_") 
                                  + (iv0 < 10 ? "00" : (iv0 < 100 ? "0" : "")) + std::to_string(iv0)
                                  + "_"
                                  + (iDt < 10 ? "00" : (iDt < 100 ? "0" : "")) + std::to_string(iDt)
                                  + ".txt";
            // open the transitory report
            std::ifstream ISimState(iSimState); 

            // reed transitory report and rewrite in final report
            std::string line;
            while (std::getline(ISimState, line)) 
            {
                FinalState << line << std::endl;
            }

            // close and remove transitory report
            ISimState.close();     
            remove(iSimState.c_str());
        }
    }

    // close final report
    FinalState.close();

    // motivational message IMPORTANT!!!
    std::cout << "\nProgram finished without errors\n";

    return 0;
}