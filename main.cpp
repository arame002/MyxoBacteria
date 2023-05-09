

#include "TissueBacteria.hpp"
//#include "SignalCommon.hpp"
#include <chrono>
#include <string.h>
#include <stdio.h>
using namespace std;





// Inverse function is on inside the InverseTime
// Function: WriteInFile, RunMode ,


void RandomForce() ;
void PositionUpdating (double ) ;                               // (time step)
void ParaView2 () ;
void InitialPili () ;
void SurfaceCoverage () ;
void VisitsPerGrid () ;
int PowerLawExponent () ;
void Update_SurfaceCoverage (ofstream ProteinLevelFile ,ofstream FrequencyOfVisit ) ;
void InitializeMatrix () ;
void initializeSlurmConfig(int argc, char* argv[]) ;

//-----------------------------------------------------------------------------------------------------
//domain and medium properties

GlobalConfigVars globalConfigVars ;
TissueBacteria tissueBacteria ;

//const double dx= 0.5 ;
//const double dy= 0.5 ;
//const int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
//const int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
//const int nz = 1 ;
//double X[nx] ;
//double Y[ny] ;
//double Z[nz]= {0.0} ;
//int visit[nx][ny] ;
vector<double> frequency ;
int numOfClasses = 1 ;
int fNoVisit = 0 ;

//-----------------------------------------------------------------------------------------------------
//slime properties
//int surfaceCoverage[nx][ny] ;
//double coveragePercentage = 0 ;


//-----------------------------------------------------------------------------------------------------
//pili properties
double piliMaxLength = 2.5 ;                       //pili maximum length
double vProPili = 0.5 ;                                 // constant protrude velocity
double vRetPili0 = 0.5 ;                                // constant retraction rate when the pili is not connected
double kPullPili = 2000 ;
double fStall = 180.0 ;
double averageLengthFree = 0.0 ;
double nAttachedPili = 0.0 ;

//-----------------------------------------------------------------------------------------------------
//simulation parameters

//int index1 = 0 ;                                     // needed for ParaView

//-----------------------------------------------------------------------------------------------------

int main (int argc, char* argv[])
{

   auto start = std::chrono::high_resolution_clock::now() ;
   
   
   initializeSlurmConfig(argc, argv);
   tissueBacteria.UpdateTissue_FromConfigFile() ;
   Fungi fungi;
  
   srand(time (0)) ;
   cout<<"Search area for slime is "<<tissueBacteria.searchAreaForSlime<<endl ;
    
    
    //---------------------------- Defining output files and directories ------------------------------
    ofstream ProteinLevelFile ;
    ProteinLevelFile.open(tissueBacteria.statsFolder + "ProteinLevelFile.txt") ;
    ofstream FrequencyOfVisit ;
    FrequencyOfVisit.open(tissueBacteria.statsFolder +"FrequncyOfVisit.txt") ;
    ofstream strSwitchP (tissueBacteria.statsFolder +"SwitchProbablities.txt" ) ;
    ofstream trajectories ( tissueBacteria.statsFolder +"trajectories.txt") ;
    
    cout<<"program is running"<<endl  ;

    //Print Header to SwitchProbabilities Data File 
    for (int i=0; i<nbacteria; i++)
    {
    ofstream strSwitchP2;
    strSwitchP2.open(tissueBacteria.statsFolder + "WriteSwitchProbabilities_Bacteria" + to_string(i) + ".txt", ios::app);
    {
       // strSwitchP << bacteria[i].motilityMetabolism.switchProbability << '\t'  ;
        strSwitchP2 << setw(10) << "maxRunDurat" << '\t' 
                    << setw(10) << "receptAct" << '\t'
                    << setw(10) << "legand" << '\t'
                    << setw(10) << "methylation" << '\t'
                    << setw(10) << "changeRate" << '\t'
                    << setw(10) << "LegandEnergy" << '\t'
                    << setw(10) << "MethylEnergy" << '\t'
                    << setw(10) << "switchProb" << '\t'
                    << setw(10) << "switchMode" << '\t';
                    
    }
    strSwitchP2<< endl ;
    }

   //--------------------------- Time controlling parameters ------------------------------------------
    double nt= tissueBacteria.runTime/tissueBacteria.dt + tissueBacteria.initialTime/tissueBacteria.initialStep ;                   //total number of steps
    nt =static_cast<int>(nt) ;
    double initialNt = tissueBacteria.initialTime/tissueBacteria.initialStep ;                     // run time for initialization
    initialNt =static_cast<int>(initialNt) ;
    //    int inverseInitialStep = static_cast<int>(initialNt/initialTime) ;  // used for visualization(ParaView)
    int inverseDt =static_cast<int>((nt-initialNt)/(tissueBacteria.runTime)) ;              // used for visualization(ParaView)
    // each frame is 0.1 second
    inverseDt = inverseDt/10 ;
   
   //--------------------------- Initializations -----------------------------------------------------
   tissueBacteria.Bacteria_Initialization() ;
   tissueBacteria.ljNodesPosition() ;
   tissueBacteria.InitialProtein() ;
   tissueBacteria.LogLinear_RNG() ;
   tissueBacteria.Update_BacteriaMaxDuration() ;
   tissueBacteria.InitialReversalTime() ;
   InitialPili () ;
   InitializeMatrix() ;
   tissueBacteria.Initialze_AllRandomForce() ;
   
   //--------------------------- Simulation setup -----------------------------------------------------
   fungi = driver(fungi) ;
   //fungi.UpdateFungiFolderNames( tissueBacteria.machineID ) ;
   vector<HyphaeSegment> hyphaeSegments_main = fungi.hyphaeSegments ;
   
   //Bacteria would try to follow hyphae as a highway
   tissueBacteria.sourceAlongHyphae = true ;
   tissueBacteria.SlimeTraceHyphae2(fungi) ;
   
   vector<vector<double> > pointSources ;
   
   // no gradient
   //pointSources.resize(2) ;
   
   // production at the tips
   pointSources = fungi.tips ;
   
   tissueBacteria.sourceProduction.resize(pointSources.at(0).size()) ;
   fill(tissueBacteria.sourceProduction.begin(), tissueBacteria.sourceProduction.end(), fungi.production) ;
   
   //fungi.WriteSourceLoc() ;
   //source for simulation 1
   //pointSources = tissueBacteria.GridSources() ;
   //pointSources = tissueBacteria.sourceChemo ;
   fungi.WriteSourceLoc( pointSources) ;
   
   tissueBacteria.gridInMain = tissueBacteria.Cal_Diffusion2D(0.0, tissueBacteria.domainx, 0.0, tissueBacteria.domainy ,tissueBacteria.tGrids.numberGridsX , tissueBacteria.tGrids.numberGridsY ,pointSources, tissueBacteria.sourceProduction ) ;
   tissueBacteria.Pass_PointSources_To_Bacteria(pointSources) ;
   
   //tissueBacteria.WriteNumberReverse() ;
   
   //--------------------------- Main loop -----------------------------------------------------------
    for (int l=0; l< (nt+1); l++)
    {
        if (l < initialNt)
        {
            
            tissueBacteria.AllNodes() ;
            tissueBacteria.Spring () ;
            tissueBacteria.Bending() ;
            //tissueBacteria.u_lj() ;
            //      RandomForce() ;
            /*
             if (l%inverseInitialStep==0)
             {
             ParaView() ;
             ParaView2() ;
             }
             */
            tissueBacteria.Update_MotilityMetabolism_Only2(.01) ;
            tissueBacteria.WriteSwitchProbabilitiesByBacteria();
            PositionUpdating(tissueBacteria.initialStep) ;
            
        }
        
        
        else
        {
            tissueBacteria.ReversalTime() ;
            tissueBacteria.AllNodes() ;
            tissueBacteria.Spring () ;
            tissueBacteria.Bending() ;
            tissueBacteria.u_lj() ;
            //  RandomForce() ;
            tissueBacteria.Motor () ;
            tissueBacteria.TurnOrientation() ;
            //  SlimeTrace() ;
            //  Connection() ;
            //  ProteinExchange() ;
            //  PiliForce() ;
            
            
            if (l%1000==0 && l%inverseDt!=0)
            {
               //Update_SurfaceCoverage(ProteinLevelFile, FrequencyOfVisit) ; //This would slow down the code bc of high number of grids
               //ParaView2() ;
               //tissueBacteria.ParaView ()  ;
               //tissueBacteria.WriteTrajectoryFile() ;
               cout<<(l-initialNt)/inverseDt<<endl ;
               //tissueBacteria.WriteBacteria_AllStats() ;
               //tissueBacteria.WriteSwitchProbabilitiesByBacteria();
                //   cout << averageLengthFree<<'\t'<<nAttachedPili<<endl ;
                //    cout << coveragePercentage <<endl ;
               
               // reducing time step in the simulation will result in changing the reversal frequencies. Lower number as inputs
               //tissueBacteria.UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
               tissueBacteria.Update_MotilityMetabolism_Only(.01) ;
            }
            
            else if  (l%inverseDt==0)
            {
               //Update_SurfaceCoverage(ProteinLevelFile, FrequencyOfVisit) ; //This would slow down the code bc of high number of grids
               ParaView2() ;
               tissueBacteria.ParaView ()  ;
               tissueBacteria.WriteTrajectoryFile() ;
               cout<<(l-initialNt)/inverseDt<<endl ;
               tissueBacteria.WriteBacteria_AllStats() ;
               tissueBacteria.WriteSwitchProbabilitiesByBacteria();
               
                //   cout << averageLengthFree<<'\t'<<nAttachedPili<<endl ;
                //    cout << coveragePercentage <<endl ;
               
               // reducing time step in the simulation will result in changing the reversal frequencies. Lower number as inputs
               tissueBacteria.UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
               tissueBacteria.Update_MotilityMetabolism(.01) ;
               
            }
            
            PositionUpdating(tissueBacteria.dt) ;
        }
        
    }
    
   tissueBacteria.WriteNumberReverse() ;
   std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
   std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
   cout << "Time to simulate "<<nbacteria<<" bacteria is : " << duration.count() << " seconds" << endl;
    cout<<"No Bug"<<endl ;
    return 0 ;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------


void PositionUpdating (double t)
{   double etaInverse = tissueBacteria.diffusion / (tissueBacteria.kblz * tissueBacteria.temp) ;
    double totalForceX = 0 ;
    double totalForceY = 0 ;
    for (int i=0; i<nbacteria; i++)
    {
       tissueBacteria.bacteria[i].oldLoc.at(0) = tissueBacteria.bacteria[i].nodes[(nnode-1)/2].x ;
       tissueBacteria.bacteria[i].oldLoc.at(1) = tissueBacteria.bacteria[i].nodes[(nnode-1)/2].y ;
       int tmpXIndex = static_cast<int>( fmod( round(fmod( tissueBacteria.bacteria[i].nodes[(nnode-1)/2].x + tissueBacteria.domainx, tissueBacteria.domainx ) / tissueBacteria.dx ), tissueBacteria.nx) ) ;
       int tmpYIndex = static_cast<int>( fmod( round(fmod( tissueBacteria.bacteria[i].nodes[(nnode-1)/2].y + tissueBacteria.domainy, tissueBacteria.domainy ) / tissueBacteria.dy ), tissueBacteria.ny )) ;
       if (tmpXIndex< 0 || tmpXIndex > tissueBacteria.nx -1 || tmpYIndex< 0 || tmpYIndex > tissueBacteria.ny -1 )
       {
          cout<<"positionUpdating, index out of range "<<tmpXIndex<<'\t'<<tmpYIndex<<endl ;
       }
       // oldChem is now updated in Cal_Chemical2
       //tissueBacteria.bacteria[i].oldChem = tissueBacteria.gridInMain.at(tmpYIndex).at(tmpXIndex) ;
        for(int j=0 ; j<nnode; j++)
        {
            tissueBacteria.Diffusion(tissueBacteria.bacteria[i].nodes[j].x , tissueBacteria.bacteria[i].nodes[j].y) ;
            etaInverse = tissueBacteria.diffusion / (tissueBacteria.kblz * tissueBacteria.temp) ;
            totalForceX = tissueBacteria.bacteria[i].nodes[j].fSpringx ;       // it is not += because the old value is related to another node
            totalForceX += tissueBacteria.bacteria[i].nodes[j].fBendingx ;
                 totalForceX += tissueBacteria.bacteria[i].nodes[j].fljx ;
            totalForceX += tissueBacteria.bacteria[i].nodes[j].fMotorx ;
           tissueBacteria.bacteria[i].nodes[j].x += t * totalForceX * etaInverse ;
           tissueBacteria.bacteria[i].nodes[j].x += tissueBacteria.bacteria[i].nodes[j].xdev ;
            
            totalForceY = tissueBacteria.bacteria[i].nodes[j].fSpringy ;       // it is not += because the old value is related to another node
            totalForceY += tissueBacteria.bacteria[i].nodes[j].fBendingy ;
                 totalForceY += tissueBacteria.bacteria[i].nodes[j].fljy ;
            totalForceY += tissueBacteria.bacteria[i].nodes[j].fMotory ;
           tissueBacteria.bacteria[i].nodes[j].y += t * totalForceY * etaInverse ;
           tissueBacteria.bacteria[i].nodes[j].y += tissueBacteria.bacteria[i].nodes[j].ydev ;
            
            
            
            // F pili
            
            if (j==0)
            {
                for (int m=0 ; m<nPili ; m++)
                {
                   tissueBacteria.bacteria[i].nodes[j].x += t * (tissueBacteria.bacteria[i].pili[m].fx) * etaInverse ;
                   tissueBacteria.bacteria[i].nodes[j].y += t * (tissueBacteria.bacteria[i].pili[m].fy) * etaInverse ;
                }
            }
           
           if (j== (nnode -1)/2 )
           {
              tissueBacteria.bacteria[i].locVelocity = etaInverse * sqrt(totalForceX * totalForceX + totalForceY * totalForceY) ;
              tissueBacteria.bacteria[i].locFriction = 1.0 / etaInverse ;
           }
            
        }
    }
   tissueBacteria.ljNodesPosition() ;
}



//-----------------------------------------------------------------------------------------------------

void ParaView2 ()
{
    int index = tissueBacteria.index1 ;
    if (index > 2)
    {
      return ;
    }
    string vtkFileName2 = tissueBacteria.folderName + "Grid"+ to_string(index)+ ".vtk" ;
    ofstream SignalOut;
    SignalOut.open(vtkFileName2.c_str());
    SignalOut << "# vtk DataFile Version 2.0" << endl;
    SignalOut << "Result for paraview 2d code" << endl;
    SignalOut << "ASCII" << endl;
    SignalOut << "DATASET RECTILINEAR_GRID" << endl;
    SignalOut << "DIMENSIONS" << " " << tissueBacteria.nx  << " " << " " << tissueBacteria.ny << " " << tissueBacteria.nz  << endl;
    
    SignalOut << "X_COORDINATES " << tissueBacteria.nx << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int i = 0; i < tissueBacteria.nx ; i++) {
        SignalOut << tissueBacteria.X[i] << endl;
    }
    
    SignalOut << "Y_COORDINATES " << tissueBacteria.ny << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int j = 0; j < tissueBacteria.ny; j++) {
        SignalOut << tissueBacteria.Y[j] << endl;
    }
    
    SignalOut << "Z_COORDINATES " << tissueBacteria.nz << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int k = 0; k < tissueBacteria.nz ; k++) {
        SignalOut << 0 << endl;
    }
    
    SignalOut << "POINT_DATA " << (tissueBacteria.nx )*(tissueBacteria.ny )*(tissueBacteria.nz ) << endl;
    SignalOut << "SCALARS liquid float 1" << endl;
    SignalOut << "LOOKUP_TABLE default" << endl;
    
    for (int k = 0; k < tissueBacteria.nz ; k++) {
        for (int j = 0; j < tissueBacteria.ny; j++) {
            for (int i = 0; i < tissueBacteria.nx; i++) {
                SignalOut << tissueBacteria.slime[i][j] << endl;
            }
        }
    }
    
}


//-----------------------------------------------------------------------------------------------------
void InitialPili ()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for ( int j=0 ; j<nPili ; j++)
        {
           tissueBacteria.bacteria[i].pili[j].lFree = rand() / (RAND_MAX + 1.0) * piliMaxLength ;
           tissueBacteria.bacteria[i].pili[j].attachment = false ;
           tissueBacteria.bacteria[i].pili[j].retraction = false ;
           tissueBacteria.bacteria[i].pili[j].fx = 0.0 ;
           tissueBacteria.bacteria[i].pili[j].fy = 0.0 ;
           tissueBacteria.bacteria[i].pili[j].F = 0.0 ;
            
            if(rand() / (RAND_MAX + 1.0) < 0.5 )
            {
               tissueBacteria.bacteria[i].pili[j].retraction = true ;
            }
            else
            {
               tissueBacteria.bacteria[i].pili[j].retraction = false ;
            }
            if (rand() / (RAND_MAX + 1.0) < 0.5 )
            {
               tissueBacteria.bacteria[i].pili[j].attachment = true ;
               tissueBacteria.bacteria[i].pili[j].retraction = true ;
            }
            else
            {
               tissueBacteria.bacteria[i].pili[j].attachment = false ;
            }
            
        }
    }
}

//-----------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void SurfaceCoverage ()
{
    double percentage = 0 ;
    
    for (int m=0; m< tissueBacteria.nx; m++)
    {
        for (int n=0 ; n < tissueBacteria.ny ; n++)
        {
            if (tissueBacteria.slime[m][n] > tissueBacteria.s0 && tissueBacteria.surfaceCoverage[m][n] == 0 )
            {
               tissueBacteria.surfaceCoverage[m][n] = 1 ;
                percentage += 1 ;
            }
        }
    }
    percentage = percentage / (tissueBacteria.nx * tissueBacteria.ny * 1.0) ;
   tissueBacteria.coveragePercentage += percentage ;
    
}

//------------------------------------------------------------------------------------------------------------
void VisitsPerGrid ()
{
    //SlimeTrace is needed for this function, it calculate #visit in each grid
    int m = 0 ;
    int n = 0 ;
    //  int newVisit[nx][ny] = {0} ;
    
    
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<2*nnode-1 ; j++)
        {
            // we add domain to x and y in order to make sure m and n would not be negative integers
            m = (static_cast<int> (round ( fmod (tissueBacteria.bacteria[i].allnodes[j].x + tissueBacteria.domainx , tissueBacteria.domainx) / tissueBacteria.dx ) ) ) % tissueBacteria.nx  ;
            n = (static_cast<int> (round ( fmod (tissueBacteria.bacteria[i].allnodes[j].y + tissueBacteria.domainy , tissueBacteria.domainy) / tissueBacteria.dy ) ) ) % tissueBacteria.ny  ;
           tissueBacteria.visit[m][n] += 1 ;
            //  newVisit[m][n] = 1 ;
            
        }
    }
    
}

//------------------------------------------------------------------------------------------------------------
int PowerLawExponent ()
{
    double discrete = 1.5 ;
    double log10discrete = log10(discrete) ;
    fNoVisit = 0 ;
    int minValue = tissueBacteria.visit[0][0] ;
    int maxValue = tissueBacteria.visit[0][0] ;
    int alfaMin = 0 ;
    int alfaMax = 0 ;
    
    // alfa and a must be integer, if you want to change them to a double variable you need to make some changes in the code
    int alfa = 0 ;
    for (int n = 0; n< tissueBacteria.ny; n++)
    {
        for (int m = 0 ; m < tissueBacteria.nx; m++)
        {
            if(tissueBacteria.visit[m][n] > maxValue ) maxValue = tissueBacteria.visit[m][n];
            if(tissueBacteria.visit[m][n] < minValue ) minValue = tissueBacteria.visit[m][n];
        }
    }
    alfaMin = max(floor (log10(minValue)/log10discrete) , 0.0 ) ;
    alfaMax = max(floor (log10(maxValue)/log10discrete) , 0.0 );
    numOfClasses = (alfaMax - alfaMin +1 )  ;
    frequency.assign(numOfClasses, 0) ;
    
    for (int n = 0; n < tissueBacteria.ny; n++)
    {
        for (int m = 0 ; m < tissueBacteria.nx; m++)
        {   if(tissueBacteria.visit[m][n] == 0)
        {
            fNoVisit += 1 ;
            continue ;
        }
            alfa = floor(log10(tissueBacteria.visit[m][n])/log10discrete) ;
            frequency.at(alfa ) += 1 ;
        }
    }
    double normalization = 1 ;
    for (int i = 0 ; i<numOfClasses ; i++)
    {
        normalization = floor(pow(discrete, i+1)) - ceil(pow(discrete, i)) + 1 ;    // it counts integers between this two numbers
        
        if(normalization == 0 )    normalization = 1 ;  // for that i, frequesncy is zero. normalized or unnormalized .
        frequency.at(i) = frequency.at(i) / normalization ;
    }
    
    return alfaMin ;
}

void InitializeMatrix ()
{
   for (int m=0; m < tissueBacteria.nx; m++)
   {
       for (int n=0; n < tissueBacteria.ny; n++)
       {
          tissueBacteria.slime[m][n] = tissueBacteria.s0  ;
          tissueBacteria.surfaceCoverage[m][n] = 0 ;
          tissueBacteria.visit[m][n] = 0 ;
       }
   }
   
   for (int i=0; i < tissueBacteria.nx; i++)
   {
      tissueBacteria.X[i]= i * tissueBacteria.dx ;
   }
   for (int j=0; j < tissueBacteria.ny; j++)
   {
      tissueBacteria.Y[j]= j * tissueBacteria.dy ;
   }
}

void initializeSlurmConfig(int argc, char* argv[]) {
   
   // read configuration.
   ConfigParser parser;
   std::string configFileNameDefault = "./resources/bacteria_M.cfg";
   globalConfigVars = parser.parseConfigFile(configFileNameDefault);
   std::string configFileNameBaseL = "./resources/test";
   std::string configFileNameBaseR = ".cfg";

   // Unknown number of input arguments.
   if (argc != 1 && argc != 3) {
      std::cout << "ERROR: Incorrect input argument count.\n"
            << "Expect either no command line argument or three arguments"
            << std::endl;
      exit(0);
   }
   // one input argument. It has to be "-slurm".
   else if (argc == 3) {
      if (strcmp(argv[1], "-slurm") != 0) {
         std::cout
               << "ERROR: one argument received from commandline but it's not recognized.\n"
               << "Currently, the argument value must be -slurm"
               << std::endl;
         exit(0);
      } else {
         std::string configFileNameM(argv[2]);
         std::string configFileNameCombined = configFileNameBaseL
               + configFileNameM + configFileNameBaseR;
         parser.updateConfigFile(globalConfigVars, configFileNameCombined);
      }
   }
   // no input argument. Take default.
}

void Update_SurfaceCoverage (ofstream ProteinLevelFile ,ofstream FrequencyOfVisit )
{
   string NumberOfVisits = tissueBacteria.statsFolder + "NumberOfVisits"+ to_string(tissueBacteria.index1)+ ".txt" ;
   ofstream NumberOfVisit;
   NumberOfVisit.open(NumberOfVisits.c_str());
   SurfaceCoverage() ;
   VisitsPerGrid() ;
   int alfaMin = PowerLawExponent() ;
   
   for (int m=0; m< tissueBacteria.nx; m++)
   {
   for (int n=0 ; n < tissueBacteria.ny ; n++)
   {
   NumberOfVisit<< tissueBacteria.visit[m][n]<<'\t' ;
   }
   NumberOfVisit<<endl ;
   }
   NumberOfVisit<<endl<<endl ;
   NumberOfVisit.close() ;
   
   for (int i=0; i<nbacteria; i++)
   {
   ProteinLevelFile<<tissueBacteria.bacteria[i].protein<<"," ;
   }
   ProteinLevelFile<< endl<<endl ;
   FrequencyOfVisit<< "Results of frame "<< tissueBacteria.index1<<" are below : "<<endl ;
   FrequencyOfVisit<< "Surface coverage is equal to:     "<<tissueBacteria.coveragePercentage <<endl ;
   FrequencyOfVisit<< "alfaMin is equal to:    " << alfaMin<<endl ;
   FrequencyOfVisit<<"size of FrequencyOfVisit is equal to:  "<< numOfClasses<<endl ;
   FrequencyOfVisit<<"Frequecny of no visit is equal to:    "<<fNoVisit<<endl ;
   
   for (int i=0 ; i<numOfClasses ; i++)
   {
   FrequencyOfVisit<< frequency[i]<<'\t' ;
   }
   FrequencyOfVisit<<endl<<endl ;
  
}

