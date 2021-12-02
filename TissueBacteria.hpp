
#ifndef TissueBacteria_hpp
#define TissueBacteria_hpp

#include <stdio.h>
#include "Bacteria.hpp"
#include "Diffusion2D.hpp"

#endif /* TissueBacteria_hpp */
class TissueBacteria
{
public:
    bacterium bacteria[nbacteria];
    TissueGrid tGrids ;
    
    double length = 2 ;
    double B = 0.5 ;                      // bending constant
    double tetta0=3.1415 ;               // prefered angle, Pi
    double K =  1.0 ;                      // linear spring constant (stiffness)
    double x0 = length/(nnode-1) ;                      // equilibrium distance of spring
    double fmotor = 0.1 ;                           //   Ft/(N-1)
    double Rp = 0.3 ;                               // protein exchange rate
    double turnPeriod = 50.0 ;
    int shiftx = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    int shifty = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
    const double dx= 0.5 ;
    const double dy= 0.5 ;
    const int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
    const int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
    
    double diffusion ;
    double dt=0.0001 ;
    double initialStep = 0.0001 ;
    int index1 = 0 ;    
    
    double sr = 0.25 ;                          // slime rate production
    double kd = 0.05 ;                           // slime decay rate
    double slimeEffectiveness = 5.0 ;           // needed for slime trail following
    double s0 = 1 ;
    int searchAreaForSlime = static_cast<int> (round((length/2.0)/ min(dx , dy))) ;
    //double slime[nx][ny] ;
    vector<vector<double> > slime ;
    vector<vector<double> > gridInMain ;
    

    double reversalRate = 0.025 ;
    double chemoStrength = 10.0 ;
    
    
    double kblz=1.0 , temp=0.01 , eta1 = 1.0 , eta2 = 5.0 ;
    

    TissueBacteria () ;
    vector<vector<double> > Cal_Diffusion2D (double xMin, double xMax, double yMin, double yMax, vector<vector<double> > sources) ;
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
    double Cal_ChemoGradient (int i) ;
    double Cal_ChemoGradient2 (int i) ;
    void WriteTrajectoryFile () ;
    void WriteNumberReverse () ;
    void SlimeTraceHyphae (Fungi tmpFng) ;

    
};