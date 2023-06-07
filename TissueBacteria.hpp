
#ifndef TissueBacteria_hpp
#define TissueBacteria_hpp

#include <stdio.h>
#include "Bacteria.hpp"
#include "Diffusion2D.hpp"
#include "ranum2.h"

#endif /* TissueBacteria_hpp */

enum ChemotacticMechanism
{
    classic = 0 ,
    metabolism = 1
    
};
enum InitialCondition
{
    uniform = 0 ,
    circular = 1 ,
    swarm = 2 ,
    center = 3,
    alongNetwork = 4
};




class TissueBacteria
{
public:
    //---------------------------- Parameters and sub-classes ------------------------------
    bacterium bacteria[nbacteria];
    TissueGrid tGrids ;

    int machineID = 3 ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    string animationName = "Bacteria_" ;
    
    bool inLiquid = true ;
    bool PBC = true ;
    ChemotacticMechanism chemotacticMechanism = classic ;
    InitialCondition initialCondition = circular ;
    bool run_calibrated = 1 ;
    double lognormal_run_m = 0.7 ;
    double lognormal_run_s = 0.6 ;
    double lognormal_run_a = 1.0 ;
    double lognormal_wrap_m = -0.4 ;
    double lognormal_wrap_s = 0.9 ;
    double lognormal_wrap_a = 1.0 ;
    std::default_random_engine run_seed ;
    std::default_random_engine wrap_seed ;
    std::lognormal_distribution<double> run_distribution ;
    std::lognormal_distribution<double> wrap_distribution ;

    double length = 2 ;
    double B = 0.5 ;                      // bending constant
    double tetta0=3.1415 ;               // prefered angle, Pi
    double K =  5.0 ;                      // linear spring constant (stiffness)
    double x0 = length/(nnode-1) ;                      // equilibrium distance of spring
    double lj_Energy = 0.0014 ;
    double lj_rMin = 0.6 ;
    bool lj_Interaction = true ;
    double updatingFrequency = 10.0 ;
    
    double fmotor = (nnode) * 0.3/nnode ;                           //   Ft/(N-1)
    double motorEfficiency = 1.0 ;
    double motorEfficiency_Liquid = 1.0 ;
    double motorEfficiency_Agar = 0.4 ;
    
    double Rp = 0.3 ;       // protein exchange rate, used for Myxo study( Abu's paper)
    long  idum = (-799);    // used for random generator in Gaussian distribution( ranum2 file)
    
    int shiftx = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    int shifty = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    
    double domainx = 1000.0 ;
    double domainy = 1000.0 ;
    double dx= 0.5 ;
    double dy= 0.5 ;
    int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
    int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
    int nz = 1 ;
    
    double initialTime = 0.5 ;
    double runTime = 200.0 ;
    double dt=0.000005 ;
    double initialStep = 0.0001 ;
    int index1 = 0 ;    
    
    double sr = 0.25 ;                          // slime rate production
    double kd = 0.05 ;                          // slime decay rate
    double slimeEffectiveness = 0.8 ;           // needed for slime trail following, used to be 5
    double liqBackground = 1.0 ;
    double liqLayer = 10.0  ;                  // liquid layer around fungi network
    double liqHyphae = 0.1 ;                    // liquid value for where hyphae is located ( physical barrier)
    int searchAreaForSlime = static_cast<int> (round(( length )/ min(dx , dy))) ;
    double Slime_CutOff = 1.3 ;
    vector<vector<double> > slime ;
    vector<vector<double> > viscousDamp ;
    vector<vector<double> > gridInMain ;
    vector<vector<double> > sourceChemo ;       //store the source locations in TissueBacteria class
    vector<double> sourceProduction ;           //store the source production rate in TissueBacteria class
    bool sourceAlongHyphae = true ;

    double turnPeriod = 0.1 ;
    double reversalRate = 1.0/ 5.0 ;         //0.025
    double minimumRunTime = 0.3 ;             //this has to be larger than turnPeriod and wrapDuration
    double chemoStrength_run = 5000.0 ;
    double chemoStrength_wrap = 5000.0 ;
    double agarThicknessX = domainx/20 ;
    double agarThicknessY = domainy/20 ;
    
    
    double eta_background = 0.2 , eta_Barrier = 2.0 ;
    
    vector<double> X ;
    vector<double> Y ;
    vector<double> Z ;
    vector<vector<int> > visit ;
    vector<vector<int> > surfaceCoverage ;
    double coveragePercentage = 0 ;
    
    //---------------------------- Functions -----------------------------------------------
    TissueBacteria () ;
    
    void Initialze_AllRandomForce () ;
    void UpdateTissue_FromConfigFile () ;
    void Bacteria_Initialization () ;
    // Bacteria_Initialization Call one of the following initialization functions:
    void Initialization () ;
    void Initialization2 () ;
    void CircularInitialization () ;
    void SwarmingInitialization () ;
    void CenterInitialization () ;
    void AlongNetworkInitialization () ;
    
    
    void Update_LJ_NodePositions () ;
    void Initialize_ReversalTimes () ;
    vector<vector<double> > TB_Cal_ChemoDiffusion2D (double xMin, double xMax, double yMin, double yMax,int nGridX, int nGridY ,vector<vector<double> > sources, vector<double> pSource) ;
    void Cal_AllLinearSpring_Forces() ;
    
    //Does not consider period boundary effect.
    double Cal_NodeToNodeDistance (int b1 ,int n1, int b2, int n2) ; // bacteria, node, bacteria, node
    void Cal_AllBendingSpring_Forces() ;
    double Cal_Cos0ijk_InteriorNodes(int i,int m) ;
    
    //Calculate point to point distance considering period boundary condition
    double Cal_PtoP_Distance_PBC (double x1, double y1, double x2, double y2, int nx , int ny) ;
    //Calculate the minimum distance of 2 points( 2 bacteria or their duplicates)
    double Cal_MinDistance_PBC (double x1 , double y1 , double x2 , double y2) ;
    
    // Update the bacteria-bacteria connection. Used for ProteinExchange() in Myxo study( Abu's paper)
    void Update_BacterialConnection () ;
    double Cal_AllBacteriaLJ_Forces() ;
    //void PiliForce () ;
    void TermalFluctiation_Forces() ;
    double Update_LocalFriction (double, double ) ;
    void Update_ViscousDampingCoeff ();
    
    //Apply equal force to each nodes. Direction is toward the next node.
    void Cal_MotorForce () ;
    //Bacteria follow the highway ( eithier liquid around hyphae of slime produced by other bacteria)
    double Bacterial_HighwayFollowing (int i, double R)  ;      // i th bacteria, radius of the search area
    
    void Reverse_IndividualBacteriaa (int i) ;
    void Check_Perform_AllReversing_andWrapping() ;
    void Update_Bacteria_AllNodes () ;          //Update allNodes with nodes and ljNodes
    void Initialize_BacteriaProteinLevel () ;
    
    void Bacterial_ProteinExchange () ; //Exchange protein between bacteria in Myxo study
    void Update_NodeLevel_ProteinConcentraion () ;
    void BacterialVisualization_ParaView () ;
   // void ParaView2 () ;
    
    // Make a duplicate of the bacteria that are leaving the domain
    void Bacteria_CreateDuplicate () ;
    // Merge the bacteria and its duplicate when it is completely in the domain
    void Bacteria_MergeWithDuplicate () ;
    
    //slime[][] would be calculated based on where bacteria were moving. Used for Myxo study
    //Include "bacterial secretion" and "degradation over time"
    void Update_SlimeConcentration () ;
    
    void Handle_BacteriaTurnOrientation () ;   // Bacteria change orientation after reverse or during the wrap mode
    double Cal_BacteriaOrientation (int i) ;
    
    //Calculate the impact of chemo-gradient on reversal based on observational model
    void UpdateReversalFrequency () ;
    
    //Each bacteria has its own max duration in absence of chemoattractant based on distribution observed in experiments
    void Update_BacteriaMaxDuration () ;
    
    //Calculate the change in chemoattractant for i-th bacteria.
    // It is used for changing the bacteria.reversalPeriod
    double Cal_ChemoGradient (int i) ;
    double Cal_ChemoGradient2 (int i) ;     //Current working version
    
    
    void WriteTrajectoryFile () ;
    void WriteNumberReverse () ;
    
    //The entire hyphae width is filled with liquid. Liquid decreass with distance from hyphae
    void Find_FungalNetworkTrace (Fungi tmpFng) ;
    
    //Update the liquid concentration based on relative location with respect to fungi network
    // This function "determine" the highway for bacterial motion near fungi
    void Find_FungalNetworkTrace2 (Fungi tmpFng) ;         // bacteria follows the surrounding.
    
    //Find the coordinates and secretion rates of chemoattractant sources based on our assumption( tips vs unirform along fungi)
    vector<vector<double> > Find_secretion_Coord_Rate (Fungi tmpFng) ;
    
    void Pass_PointSources_To_Bacteria (vector<vector<double> > sourceP) ;
    
    //write information of the baceria art the current time including its physical and chemical environment
    void WriteBacteria_AllStats () ;
    void LogLinear_RNG () ;     // Used to calibrate durations to experimental data
    
    //---------------------------- Metabolism Functions -----------------------------------------------
    void Update_MotilityMetabolism (double tmpDt) ;
    void Update_MotilityMetabolism_Only (double tmpDt) ;
    void Update_MotilityMetabolism_Only2 (double tmpDt) ;
    void WriteSwitchProbabilities () ;
    void WriteSwitchProbabilitiesByBacteria();
    void Update_MM_Legand () ;

    
};
