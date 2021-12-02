

#ifndef TissueGrid_hpp
#define TissueGrid_hpp

#include "Grid.hpp"

class TissueGrid
{
public:
    TissueGrid () ;
    vector< vector<Grid> > grids ;
    double xDomainMin = 0 ;
    double xDomainMax = 100 ;
    double yDomainMin = 0 ;
    double yDomainMax = 100 ;
    int numberGridsX = 200 ;
    int numberGridsY = 200 ;
    double grid_dx ;
    double grid_dy ;
    double Diffusion = 200.0 ;
    double deg = 0.0001 ;
    double pro = pow(10, 2) ;
    double grid_dt = 0.05 ;      //timeStep
    vector<double> xSources ;
    vector<double> ySources ;
    vector<int> indexSourceX ;
    vector<int> indexSourceY ;
    
    void DomainBoundaries (double xMin, double xMax, double yMin, double yMax, int nGridX , int nGridY) ;
    void InitializeAllGrids () ;
    void ClearChanges () ;
    void DiffusionChanges () ;
    void FindProductionPoints () ;
    void ProductionChanges () ;
    void DegredationChanges () ;
    void UpdateChanges () ;
    void EulerMethod () ;
    void ParaViewGrids (int) ;
    void RoundToZero () ;
    
    
};
#endif /* TissueGrid_hpp */
