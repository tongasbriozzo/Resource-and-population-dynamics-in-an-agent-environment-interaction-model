/******************************************************************************/
/* Functions                                                                  */   
/* All model functions are defined in this file                               */
/******************************************************************************/

#include "classes.h"
using namespace std;

/******************************************************************************/
/* Distances                                                                  */
/******************************************************************************/

/* Distancia en x entre dos agentes */
float DistanceX(Agent &A, Agent &B)
{
    float dx = A.x - B.x;
    dx -= Lx*round(dx/Lx);
    if ( ! (-Lx<=dx && dx<=Lx ) ) std::cout << std::endl << "WARNING! Distance x: " << dx << std::endl; // WARNING
    return dx;
}

/* Distancia en y entre dos agentes */
float DistanceY(Agent &A, Agent &B)
{
    float dy = A.y - B.y;
    dy -= Ly*round(dy/Ly);
    if ( ! (-Ly<=dy && dy<=Ly ) ) std::cout << std::endl << "WARNING! Distance y: " << dy << std::endl; // WARNING
    return dy;
}

/* Distancia total entre dos agentes */
float Distance(Agent &A, Agent &B)
{
    float dx = DistanceX(A, B);
    float dy = DistanceY(A, B);
    float dd = sqrt((dx*dx)+(dy*dy));
    return dd;
}

/******************************************************************************/
/* Evolution of Inner Energy                                                  */
/******************************************************************************/

/* Un agente se alimenta de un parche */
void Feed(Agent &newAgent, Agent &oldAgent, Patch &newPatch, Patch &oldPatch)
{
    float foodIntake = 0;

    if (!saturation) foodIntake = timeStep*oldAgent.Is*oldPatch.f; // Sin saturacion, el intake es lineal con I_s
    if ( saturation) foodIntake = timeStep*oldAgent.Ic*tanh(oldPatch.f*oldAgent.Is/oldAgent.Ic); // Con saturacion, el intake satura a I_c

    newAgent.e += foodIntake;
    newPatch.f -= foodIntake;
}

/* Evolucion de la energia de un agente */
void UpdateInnerEnergy(Agent &newAgent, Agent &oldAgent)
{  
    newAgent.e -= timeStep*(oldAgent.m*oldAgent.e + oldAgent.k*pow(oldAgent.v0,2));;
}

/* Evolucion de la energia de un parche */
void UpdateInnerEnergy(Patch &newPatch, Patch &oldPatch)
{
    newPatch.f += timeStep*oldPatch.r*oldPatch.f*(1.0 - oldPatch.f/oldPatch.c);
}

/* Evolucion de la energia de todo el sistema agentes + parches */
void UpdateInnerEnergy(vector<Agent> &newSystem, vector<Agent> &oldSystem, vector<Patch> &newBoard, vector<Patch> &oldBoard)
{
    for (int i=0; i<M2; i++)
    {
        UpdateInnerEnergy(newBoard[i], oldBoard[i]);
    }

    for (int i=0; i<newSystem.size(); i++)
    {
        int j = oldSystem[i].GetPatchIndex();
        Feed(newSystem[i], oldSystem[i], newBoard[j], oldBoard[j]);
        UpdateInnerEnergy(newSystem[i], oldSystem[i]);
    }
}

/******************************************************************************/
/* Patch Index                                                                */
/******************************************************************************/

/* Dados los indices x y del parche, regresa su indice global */
int PatchIndex(int ix, int iy)
{
    ix -= Mx*floor(float(ix)/float(Mx));
    iy -= My*floor(float(iy)/float(My));
    int iPatch = ix + iy*Mx;
    if ( ! ( 0<=iPatch && iPatch<M2 ) ) std::cout << std::endl << "WARNING! Patch Index: " << iPatch << std::endl; // WARNING
    return iPatch;
}

/* Dado el indice global del parche, regresa sus indices x y */
void PatchIndexs(int iP, int &ix, int &iy)
{
    if ( ! ( 0<=iP && iP<M2 ) ) std::cout << std::endl << "WARNING! Patch Indexs: " << iP << std::endl; // WARNING
    ix = iP%Mx;
    iy = int(floor(float(iP)/float(Mx)))%My;
}

/******************************************************************************/
/* Linked List Grid                                                           */
/******************************************************************************/

/* Dadas las coordenadas x y de un agente, regresa el indice de la grilla donde se encuentra */
int GridIndex(KIND x, KIND y)
{
    int iGrid = floor(x/Lx * lxGrid) + floor(y/Ly * lyGrid)*lxGrid;
    if ( ! ( 0<=iGrid && iGrid<nGrid ) ) std::cout << std::endl << "WARNING! Grid Index: " << iGrid << std::endl; // WARNING
    return iGrid;
}

/* Dado un sistema de agentes, ubica sus indices en la grilla */
void UpdateGrid(vector<Agent> &system, vector<vector<set<size_t>>> &grid)
{
	for (int ig=0; ig<lxGrid; ig++) for (int jg=0; jg<lyGrid; jg++)
    {
        grid[ig][jg].clear();
    }
  	for(int i=0; i<system.size(); i++)
	{
        int ig = floor(system[i].x/Lx * lxGrid);
        if ( ! ( 0<=ig && ig<lxGrid ) ) std::cout << std::endl << "WARNING! Grid Index I: " << ig << std::endl; // WARNING
        int jg = floor(system[i].y/Ly * lyGrid);
        if ( ! ( 0<=jg && jg<lyGrid ) ) std::cout << std::endl << "WARNING! Grid Index J: " << jg << std::endl; // WARNING
        grid[ig][jg].insert(i);
    }
}

/******************************************************************************/
/* Food Field                                                                 */
/******************************************************************************/

/* argtan function */
float arg_tan(float x, float y)
{
    if ( ! ( (x<=0||0<=x) && (y<=0||0<=y) ) ) std::cout << std::endl << "WARNING! ArgTan xy: " << x << ", " << y << std::endl; // WARNING
    float a = 0;
    if (x<0 && y< 0) a = - pi + atan(y/x);
    if (x<0 && y>=0) a =   pi + atan(y/x);
    if (x>0 && y< 0) a =        atan(y/x);
    if (x>0 && y>=0) a =        atan(y/x);
    if (x==0)
    {
        if (y<0) a = - pi/2;
        if (y>0) a =   pi/2;
    }    
    if ( ! ( -pi<=a && a<=pi ) ) std::cout << std::endl << "WARNING! ArgTan a: " << a << std::endl; // WARNING
    return a;
}

/* Dado un parche, regresa el angulo del gradiente de alimentos */
float FoodGradientAngle(int iPatch, vector<Patch> &board)
{
    if ( ! ( 0<=iPatch && iPatch<M2 ) ) std::cout << std::endl << "WARNING! FGA Patch Index: " << iPatch << std::endl; // WARNING
    int xPatch = 0;
    int yPatch = 0;
    PatchIndexs(iPatch, xPatch, yPatch);

    float FoodField[2];
    FoodField[0] = 0;
    FoodField[1] = 0;

    /* Seria interesante agregar al agente una variavle interna, lengthOfSight, que determine hasta que distancia puede sensar el ambiente, para evolucion */
    for (int iy=-2; iy<3; iy++)
    {
        for (int ix=-2; ix<3; ix++)
        {
            if (ix==0 && iy==0) continue;
            int jPatchX = xPatch+ix;
            int jPatchY = yPatch+iy;
            int jPatch  = PatchIndex(jPatchX, jPatchY);
            float f = board[jPatch].f;
            if ( ! ( 0<=f && f<=board[jPatch].c ) ) std::cout << std::endl << "WARNING! FGA Agent Energy: " << f << std::endl; // WARNING

            FoodField[0] += f * ix / ( pow(ix,2) + pow(iy,2) );
            FoodField[1] += f * iy / ( pow(ix,2) + pow(iy,2) );
        }
    }
    return arg_tan(FoodField[0], FoodField[1]);
}

/* Dado un sistema de agentes y parches, los agentes performan chemotaxis */
void FoodRotation(vector<Agent> &system, vector<Patch> &board)
{
  for (int iAgent=0; iAgent<system.size(); iAgent++)
  {
    int iPatch = system[iAgent].GetPatchIndex();

    float maxFoodAngle = FoodGradientAngle(iPatch, board);
    float     oldAngle = system[iAgent].a;
    float     rotAngle = maxFoodAngle - oldAngle;

    float newAngle = oldAngle + timeStep * system[iAgent].g * sin(rotAngle);

    /* WARNING */
    if ( ! ( -10*pi<=newAngle && newAngle<=10*pi ))
    {
        std::cout << std::endl << "WARNING! FR New Angle: " << iPatch << ", "  << maxFoodAngle << ", "  << oldAngle << ", "  << rotAngle << ", " << newAngle << std::endl;
        newAngle -= twoPi*round(newAngle/twoPi);
    }    

    while (newAngle >   pi) newAngle -= twoPi;
    while (newAngle <= -pi) newAngle += twoPi;

    system[iAgent].a = newAngle;
  }
}

/* Dado un sistema de agentes y parches, los agentes rebotan al pasar a un parche con menos alimentos */
void FoodReflecion(vector<Agent> &newSystem, vector<Agent> &oldSystem, vector<Patch> &board)
{
  for (int iAgent=0; iAgent<newSystem.size(); iAgent++)
  {
    int oldPI = oldSystem[iAgent].GetPatchIndex();
    int newPI = newSystem[iAgent].GetPatchIndex();

    if (oldPI==newPI) continue;

    float oldF = board[oldPI].f;
    float newF = board[newPI].f;

    if (newF < oldF)
    {
      float a = newSystem[iAgent].a + pi;
      a -= twoPi*round(a/twoPi);
      while (a >   pi) a -= twoPi;
      while (a <= -pi) a += twoPi;
      newSystem[iAgent].a = a;
    }
  }
}

/******************************************************************************/
/*Spatial Evolution of the System*/
/******************************************************************************/

/* Emplea linked list para obtener las fuerzas y presiones a las que esta sometido cada agente */
void UpdateSystemField(vector<Agent> &system, vector<vector<float>> &field,
                       vector<vector<set<size_t>>> &grid)
{
    for (int iGrid=0; iGrid<lxGrid; iGrid++) for (int jGrid=0; jGrid<lyGrid; jGrid++) // elegimos una grilla
    {        
        for (auto& iAgent: grid[iGrid][jGrid])                                         // de esa grilla, elegimos un agente
        {   
            for (auto& ij : ijGrid)                   // elegimos otra grilla entre las vecinas 
            {   
                int iGrid2 = iGrid + ij[0];
                int jGrid2 = jGrid + ij[1];  
                while (iGrid2 <       0) iGrid2 += lxGrid;
                while (lxGrid <= iGrid2) iGrid2 -= lxGrid;
                while (jGrid2 <       0) jGrid2 += lyGrid;
                while (lyGrid <= jGrid2) jGrid2 -= lyGrid;              
                for (auto& jAgent: grid[iGrid2][jGrid2])                               // elegimos un agente vecino
                {                    
                    if (iGrid2==iGrid && jGrid2==jGrid && iAgent>=jAgent) continue;
                    float dd = Distance(system[iAgent], system[jAgent]);
                    if (dd!=0 && dd < rCut)
                    {
                        if (dd < radius*0.0001) std::cout << "WARNING! dd = " << dd << std::endl;

                        float d2 = pow(dd,-2);
                            
                        float dx = DistanceX(system[iAgent], system[jAgent])/dd;
                        float dy = DistanceY(system[iAgent], system[jAgent])/dd;

                        field[iAgent][0] += d2*dx;
			            field[iAgent][1] += d2*dy;
                        field[jAgent][0] -= d2*dx;
			            field[jAgent][1] -= d2*dy;                            
                    }                    
                }
            }
        }
    }
}

/* Dado un sistema de agentes, actualiza sus posiciones */
void UpdateSystemPositions(vector<Agent> &newSystem, vector<Agent> &oldSystem,
                           vector<vector<set<size_t>>> grid)
{
    vector<vector<float>> field;
    field.resize(newSystem.size());
    for (int i=0; i<newSystem.size(); i++)
    {
        field[i].resize(2);
        field[i][0] = 0;
        field[i][1] = 0;
    }

    UpdateSystemField(oldSystem, field, grid);

    for (int i=0; i<newSystem.size(); i++)
    {
        float v = newSystem[i].v0;
        float a = newSystem[i].a;
        float vx = v*cos(a) + field[i][0]*v*gammaFriction;
        float vy = v*sin(a) + field[i][1]*v*gammaFriction;
        newSystem[i].vx = vx;
        newSystem[i].vy = vy;

        float vv = sqrt(vx*vx+vy*vy);
        if ( ! (vv<radius/timeStep) )
        {
            vx = vx/vv * radius/timeStep;
            vy = vy/vv * radius/timeStep;
        }
        newSystem[i].x = oldSystem[i].x + vx*timeStep;
        newSystem[i].y = oldSystem[i].y + vy*timeStep;
        newSystem[i].PeriodicBoundaryConditions();
        newSystem[i].RandomTumble();
    }
}

/******************************************************************************/
/*Births and Deaths*/
/******************************************************************************/

/* Dado un sistema de agentes, aplicamos nacimientos y muertes */
void BirthsAndDeaths(vector<Agent> &newSystem, vector<Agent> &oldSystem)
{        
    oldSystem.resize(0);
	int nCounter = 0;
	for (int i=0; i<newSystem.size(); i++)
	{
		float e  = newSystem[i].e;
        float eb = newSystem[i].eb;
        float es = newSystem[i].es;
        float er = newSystem[i].er;
		if (es<e && e<er)
		{
			oldSystem.push_back(newSystem[i]);
			nCounter += 1;
		}
		if (er<=e)
		{
			newSystem[i].e = eb;
		    Agent newAgent = newSystem[i];
			float angle = dis(gen)*twoPi;
			newAgent.x     += 0.5*cos(angle);
			newAgent.y     += 0.5*sin(angle);
			newSystem[i].x -= 0.5*cos(angle);
			newSystem[i].y -= 0.5*sin(angle);
			newAgent.PeriodicBoundaryConditions();
			newSystem[i].PeriodicBoundaryConditions();
			oldSystem.push_back(newSystem[i]);
			oldSystem.push_back(newAgent);
			nCounter += 2;
		}
	}
	newSystem.resize(0);
	for (int i=0; i<nCounter; i++)
	{
		Agent newAgent = oldSystem[i];
		newSystem.push_back(newAgent);
	}       
}

/******************************************************************************/
/* Evolution of the System                                                    */
/******************************************************************************/

/* Dado un sistema de agentes y parches, evolucionamos en un time step */
void UpdateSystem(vector<Agent> &newSystem, vector<Agent> &oldSystem, vector<Patch> &newBoard, vector<Patch> &oldBoard,
                  vector<vector<set<size_t>>> &grid)
{
    //std::cout << std::endl << "UpdateSystem: Flag 1" << std::endl;
    UpdateGrid           (           oldSystem,                     grid);
    //std::cout << std::endl << "UpdateSystem: Flag 2" << std::endl;
    if (foodRotation)
    {
    FoodRotation         (newSystem,                      oldBoard      );
    }
    //std::cout << std::endl << "UpdateSystem: Flag 3" << std::endl;
    UpdateSystemPositions(newSystem, oldSystem,                     grid);
    //std::cout << std::endl << "UpdateSystem: Flag 4" << std::endl;             
    UpdateInnerEnergy    (newSystem, oldSystem, newBoard, oldBoard      );    
    //std::cout << std::endl << "UpdateSystem: Flag 5" << std::endl;
    //std::cout << std::endl << "UpdateSystem: Flag 6" << std::endl;
    if (birthsAndDeaths)
    {
    BirthsAndDeaths      (newSystem, oldSystem);
    }
    if (!birthsAndDeaths)
    {
    oldSystem = newSystem;
    }
    //std::cout << std::endl << "UpdateSystem: Flag 7" << std::endl;
    oldBoard = newBoard;
    //std::cout << std::endl << "UpdateSystem: Flag 8" << std::endl;
    UpdateGrid           (newSystem,                                grid);
    //std::cout << std::endl << "UpdateSystem: Flag 9" << std::endl;
}

/******************************************************************************/
/* AQUI TERMINAN LAS FUNCIONES PRINCIPALES DEL MODELO                         */
/******************************************************************************/

/******************************************************************************/
/* Funciones Auxiliares                                                       */
/******************************************************************************/

/* Dado un sistema de agentes y parches, devuelve la energia media de estos */
void AverageInnerEnergy(vector<Agent> &system, vector<Patch> &board, vector<float> &averageInnerEnergy)
{
    averageInnerEnergy[0] = 0.0;
    averageInnerEnergy[1] = 0.0;

    int nAgents = system.size();
    if (0<nAgents)
    {
        for (int i=0; i<nAgents; i++)
        {
            averageInnerEnergy[0] += system[i].e;
        }
        averageInnerEnergy[0] /= float(nAgents);
    }

    for (int i=0; i<M2; i++)
    {
        averageInnerEnergy[1] += board[i].f;
    }    
    averageInnerEnergy[1] /= float(M2);
}

/******************************************************************************/
/* System Spatial Distribution                                                */
/******************************************************************************/

/* Dado un sistema de agentes, regresa su velocidad cuadratica media */
float VelocityRMS(vector<Agent> &system)
{
    float vrms = 0;
    for (int i=0; i<system.size(); i++)
    {
        float vx = system[i].vx;
        float vy = system[i].vy;
        vrms += vx*vx + vy*vy;
    }
    vrms = vrms/system.size();
    return vrms;
}

/* Dado un sistema de agentes, devuelve el numero medio de primeros vecinos */
float NumberFirstNeighbors(vector<Agent> &system, vector<vector<set<size_t>>> &grid)
{
    float NFN = 0;
    for (int iGrid=0; iGrid<lxGrid; iGrid++) for (int jGrid=0; jGrid<lyGrid; jGrid++)
    {
        for (auto iAgent: grid[iGrid][jGrid])
        {
            for (int iG=-1; iG<2; iG++) for (int jG=-1; jG<2; jG++)
            {
                int iGrid2 = iGrid + iG;
                int jGrid2 = jGrid + jG;
                iGrid2 -= lxGrid*floor(float(iGrid2)/float(lxGrid));
                jGrid2 -= lyGrid*floor(float(jGrid2)/float(lyGrid));
                for (auto jAgent: grid[iGrid2][jGrid2])
                {
                    if (iAgent!=jAgent)
                    {
                        float dd = Distance(system[iAgent], system[jAgent]);
                        if (dd!=0 && dd < rCut)
                        {
                            NFN += 1;
                        }
                    }
                }
            }
        }
    }
    NFN = NFN/system.size();
    return NFN;
}

/* Dado un sistema de agentes, devuelve el orden polar medio */
float PolarOrder(vector<Agent> system)
{
    float px = 0.0;
    float py = 0.0;

    int nAgents = system.size();
    if (0<nAgents)
    {
        for (int i=0; i<nAgents; i++)
        {
            float vx = system[i].vx;
            float vy = system[i].vy;
            float vt = sqrt( vx*vx + vy*vy );
            if (vt != 0.0)
            {
                px += vx/vt;
                py += vy/vt;
            }
        }

        px = px/float(nAgents);
        py = py/float(nAgents);
    }

    float p = sqrt( px*px + py*py ) ;

    return p;
}

/* Dado un sistema de agentes, devuelve el orden nematico medio */
float NematicOrder(vector<Agent> system)
{
    float px = 0.0;
    float py = 0.0;
    
    int nAgents = system.size();
    if (0<nAgents)
    {
        for (int i=0; i<nAgents; i++)
        {
            float vx = system[i].vx;
            float vy = system[i].vy;
            float tt = arg_tan(vx, vy);
            px += cos(2.0*tt);
            py += sin(2.0*tt);
        }
        px = px/float(nAgents);
        py = py/float(nAgents);
    }

    float p = sqrt( px*px + py*py ) ;

    return p;
}

/* Dado un sistema de agentes, devuelve el momento de orden nn medio */
float MomentumOrder(int nn, vector<Agent> system)
{
    float px = 0;
    float py = 0;
    
    int nAgents = system.size();
    if (0<nAgents)
    {
        for (int i=0; i<nAgents; i++)
        {
            float vx = system[i].vx;
            float vy = system[i].vy;
            float tt = arg_tan(vx, vy) * nn;
            px += cos(tt);
            py += sin(tt);
        }
        px = px/nAgents;
        py = py/nAgents;
    }

    float p = sqrt( px*px + py*py ) ;

    return p;
}

/* Dado un sistema de agentes, devuelve la entropia media */
float Entropy(vector<Agent> &system)
{
    int nAgents = system.size();
    float s = 0.0;

    if (0 < nAgents)
    {
        vector<int> sBoard;
        sBoard.resize(M2);
        for (int i=0; i<M2; i++)
        {
            sBoard[i] = 0;
        }

        for (int i=0; i<nAgents; i++)
        {
            int ind = system[i].GetPatchIndex();
            if (ind<0 || ind>=M2)
            {
                std::cout << "ENTROPY ERROR: " << ind << system[i].x << system[i].y << std::endl;
            }
            if (ind>=0 && ind<M2)
            {
            sBoard[ind] += 1;
            }
        }

        for (int i=0; i<M2; i++)
        {
            if (sBoard[i]!=0)
            {
              float ni = float(sBoard[i]) / float(nAgents);
              float si = ni * log(ni);
              s += si;
            }
        }

        s = -s/log(float(M2));
    }    

    return s;
}

/* Dado un sistema de agentes viejo y uno nuevo, devuelve el desplazamiento cuadratico medio entre estos */
/* WARNING: no usar si hay nacimientos y muertes, ya que el numero de agentes no coincide */
float MeanSquaredDisplacement(vector<Agent> &newSystem, vector<Agent> &oldSystem)
{
    int nA = newSystem.size();
    float msd = 0;

    for (int i=0; i<nA; i++)
    {
        float dd = Distance(newSystem[i], oldSystem[i]);
        msd += dd*dd;
    }

    msd = msd/nA;

    return msd;
}

/* Dado un sistema de agentes, devuelve la posicion del centro de masa */
vector<float> MassCenter(vector<Agent> &system)
{
    float x = 0;
    float y = 0;
    int N_agents = system.size();

    for (int i=0; i<N_agents; i++)
    {
        x += system[i].x;
        y += system[i].y;
    }

    x = x/N_agents - Lx/2;
    y = y/N_agents - Ly/2;

    return {x,y};
}

/* Dado un sistema de agentes, devuelve la posicion de mayor densidad */
vector<int> MaxDenPos(vector<Agent> &system, int n_d)
{
    int n_dx = n_d * proportion;
    int n_dy = n_d;
    vector<int> X(n_dx,0);
    vector<int> Y(n_dy,0);
    int N_agents = system.size();

    for (int i=0; i<N_agents; i++)
    {
        int xInd = floor( system[i].x * float(n_dx) / Lx );
        int yInd = floor( system[i].y * float(n_dy) / Ly );
        X[xInd] += 1;
        Y[yInd] += 1;
    }

    int xMax = 0;
    int yMax = 0;
    int xIM  = 0;
    int yIM  = 0;
    for (int i=0; i<n_dx; i++)
    {
        if (X[i] > xMax)
        {
            xMax = X[i];
            xIM  = i;
        }
    }    
    for (int i=0; i<n_dy; i++)
    {
        if (Y[i] > yMax)
        {
            yMax = Y[i];
            yIM  = i;
        }
    }
    return {xIM,yIM};
}

/******************************************************************************/
/* Thermalize system                                                          */
/******************************************************************************/

/* Dado un sistema de agentes y parches, termaliza un time step y reporta los parametros deseados */
void ThermalizeSystem(int timeTherm, vector<Agent> &newSystem, vector<Agent> &oldSystem, vector<Patch> &newBoard, vector<Patch> &oldBoard,
                  vector<vector<set<size_t>>> &grid, bool report=true)
{
    //if ( ! (0<=newSystem.size() && newSystem.size()<=10*NStarvation )) std::cout << std::endl << "WARNING! TS new system size: " << newSystem.size() << std::endl;
    //if ( ! (0<=oldSystem.size() && oldSystem.size()<=10*NStarvation )) std::cout << std::endl << "WARNING! TS old system size: " << oldSystem.size() << std::endl;
    
    if ( ! (newBoard.size() == M2 )) std::cout << std::endl << "WARNING! TS new board size: " << newBoard.size() << std::endl;
    if ( ! (oldBoard.size() == M2 )) std::cout << std::endl << "WARNING! TS old board size: " << oldBoard.size() << std::endl;
    
    for (int time=0; time<timeTherm; time++)
    {
        if (newSystem.size()==0) 
        {
            for (int iPatch=0; iPatch<M2; iPatch++)
            {
                newBoard[iPatch].f = newBoard[iPatch].c;
                oldBoard[iPatch].f = oldBoard[iPatch].c;
            }
            break;
        } 
        UpdateSystem(newSystem, oldSystem, newBoard, oldBoard, grid);
    }

    if (report)
    {
        vector<float> aIE;        
        aIE.resize(2);
        AverageInnerEnergy(newSystem, newBoard, aIE);
        float nematicOrder = NematicOrder(newSystem);
        float polarOrder   = PolarOrder(newSystem);
        float entropy      = Entropy(newSystem);
        std::cout << "Thermalization:" << " " << newSystem.size() << " " << aIE[0] << " " << aIE[1] << " " << nematicOrder << " " << polarOrder << " " << entropy << std::endl;
    }
}

/* Dado un sistema de agentes, agrega o elimina aleatoriamente un numero de agentes */
void ShuffleSystem(float amplitud, vector<Agent> &newSystem, vector<Agent> &oldSystem)
{
    int nAgents = newSystem.size();
    int nNoise  = ceil( nAgents * amplitud * dis(gen) );
    if (nNoise<=0) nNoise = 1;

    float newEnergy = 1000;
    float r = dis(gen);
    if (r < 0.5) newEnergy = -1000;

    vector<int> indices(nAgents);
    for (int i = 0; i < nAgents; ++i) indices[i] = i;
    shuffle(indices.begin(), indices.end(), gen);

    for (int i = 0; i < nNoise; ++i) 
    {
        int indice = indices[i];
        newSystem[indice].e = newEnergy;
        oldSystem[indice].e = newEnergy;
    }
}

/* Dado un sistema de agentes, agrega o elimina aleatoriamente un numero de agentes */
void RandomBnD(float amplitud, vector<Agent> &newSystem, vector<Agent> &oldSystem)
{
    int nAgents = newSystem.size();
    int nNoise  = ceil( nAgents * amplitud * dis(gen) );
    if (nNoise<=0) nNoise = 1;

    for (int i = 0; i < nNoise; i++)
    {
        int iAgent = floor( nAgents * dis(gen) );
        if ( ! (0<=iAgent && iAgent<nAgents) ) continue;

        float newEnergy = 1000;
        float r = dis(gen);
        if (r < 0.5) newEnergy = -1000;
        
        newSystem[iAgent].e = newEnergy;
        oldSystem[iAgent].e = newEnergy;
    } 
}

/* Si el sistema de agentes esta muerto, lo resucita */
void ResurrectSystem(int amplitud, vector<Agent> &newSystem, vector<Agent> &oldSystem)
{
    if (newSystem.size() == 0)
    {
        int nAgents = ceil( amplitud * dis(gen) );
        if (nAgents<=0) nAgents=1;
        newSystem.resize(nAgents);
        oldSystem.resize(nAgents);
        oldSystem = newSystem;
    }
}

/* Randomizamos el sistema */
void RandomizeSystem(int ampRes, float ampRan, vector<Agent> &newSystem, vector<Agent> &oldSystem)
{
    int nAgents = newSystem.size();
    if (nAgents==0) ResurrectSystem(ampRes, newSystem, oldSystem);
    if (nAgents >0) RandomBnD(      ampRan, newSystem, oldSystem);
}