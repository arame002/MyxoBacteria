
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
    double diffusion ;
    double dt=0.0001 ;
    

    double reversalRate = 0.025 ;
    double chemoStrength = 10.0 ;
    
    
    double kblz=1.0 , temp=0.01 , eta1 = 1.0 , eta2 = 5.0 ;
    

    
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

    
};
