
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
    center = 3
};




class TissueBacteria
{
public:
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
    
    double fmotor = (nnode) * 0.3/nnode ;                           //   Ft/(N-1)
    double motorEfficiency = 1.0 ;
    double motorEfficiency_Liquid = 1.0 ;
    double motorEfficiency_Agar = 0.4 ;
    
    double Rp = 0.3 ;                               // protein exchange rate
    long  idum = (-799);                            // used for random generator in Gaussian distribution
    
    int shiftx = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    int shifty = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    
    double domainx = 1000.0 ;
    double domainy = 1000.0 ;
    double dx= 0.5 ;
    double dy= 0.5 ;
    int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
    int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
    int nz = 1 ;
    
    double diffusion ;
    double initialTime = 4.0 ;
    double runTime = 200.0 ;
    double dt=0.000005 ;
    double initialStep = 0.0001 ;
    int index1 = 0 ;    
    
    double sr = 0.25 ;                          // slime rate production
    double kd = 0.05 ;                           // slime decay rate
    double slimeEffectiveness = 0.8 ;           // needed for slime trail following, used to be 5
    double s0 = 1.0 ;
    int searchAreaForSlime = static_cast<int> (round(( length )/ min(dx , dy))) ;
    double Slime_CutOff = 1.3 ;
    //double slime[nx][ny] ;
    vector<vector<double> > slime ;
    vector<vector<double> > gridInMain ;
    vector<vector<double> > sourceChemo ;
    vector<double> sourceProduction ;
    bool sourceAlongHyphae = false ;

    double turnPeriod = 0.1 ;
    double reversalRate = 1.0/ 5.0 ;         //0.025
    double minimumRunTime = 0.3 ;             //this has to be larger than turnPeriod and wrapDuration
    double chemoStrength_run = 5000.0 ;
    double chemoStrength_wrap = 5000.0 ;
    double agarThicknessX = domainx/20 ;
    double agarThicknessY = domainy/20 ;
    
    
    double kblz=1.0 , temp=1.0 , eta1 = 0.02 , eta2 = 4.0 ;
    
    vector<double> X ;
    vector<double> Y ;
    vector<double> Z ;
    vector<vector<int> > visit ;
    vector<vector<int> > surfaceCoverage ;
    double coveragePercentage = 0 ;
    

    TissueBacteria () ;
    
    void Initialze_AllRandomForce () ;
    void UpdateTissue_FromConfigFile () ;
    void Bacteria_Initialization () ;
    void Initialization () ;
    void Initialization2 () ;
    void CircularInitialization () ;
    void SwarmingInitialization () ;
    void CenterInitialization () ;
    void ljNodesPosition () ;
    void InitialReversalTime () ;
    vector<vector<double> > Cal_Diffusion2D (double xMin, double xMax, double yMin, double yMax,int nGridX, int nGridY ,vector<vector<double> > sources, vector<double> pSource) ;
    void Myxo () ;
    void Spring() ;
    double Distance (int ,int, int, int) ;                                     // bacteria, node, bacteria, node
    void Bending() ;
    double Cos0ijk(int i,int m) ;
    double Distance2 (double x1, double y1, double x2, double y2, int nx , int ny) ;        // can be moved somewhere else
    double MinDistance (double x1 , double y1 , double x2 , double y2) ;
    void Connection () ;
    double u_lj() ;
    //void PiliForce () ;
    void RandomForce() ;
    void Diffusion (double, double ) ;
    void Motor () ;
    double SlimeTrailFollowing (int i, double R)  ;      // i th bacteria, radius of the search area
    void Reverse (int i) ;
    void ReversalTime() ;
    void AllNodes () ;
    void InitialProtein () ;
    void ProteinExchange () ;
    void NodeProtein () ;
    void ParaView () ;
   // void ParaView2 () ;
    void Duplicate () ;
    void Merge () ;
    void SlimeTrace () ;
    void TurnOrientation () ;
    double Cal_OrientationBacteria (int i) ;
    void UpdateReversalFrequency () ;
    void Update_BacteriaMaxDuration () ;
    double Cal_ChemoGradient (int i) ;
    double Cal_ChemoGradient2 (int i) ;
    void WriteTrajectoryFile () ;
    void WriteNumberReverse () ;
    void SlimeTraceHyphae (Fungi tmpFng) ;
    void SlimeTraceHyphae2 (Fungi tmpFng) ;         // bacteria follows the surrounding. Incompelete
    vector<vector<double> > GridSources () ;
    void Pass_PointSources_To_Bacteria (vector<vector<double> > sourceP) ;
    
    void Update_MotilityMetabolism (double tmpDt) ;
    void Update_MotilityMetabolism_Only (double tmpDt) ;
    void Update_MotilityMetabolism_Only2 (double tmpDt) ;
    void WriteSwitchProbabilities () ;
    void WriteSwitchProbabilitiesByBacteria();
    void Update_MM_Legand () ;
    
    void WriteBacteria_AllStats () ;
    void LogLinear_RNG () ;

    
};
