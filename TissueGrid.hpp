

#ifndef TissueGrid_hpp
#define TissueGrid_hpp

#include "Grid.hpp"

enum ChemoBoundaryCondition
{
    NFB = 0 ,
    FluxLeaving = 1
    
};

class TissueGrid
{
public:
    //---------------------------- Parameters and sub-classes ------------------------------
    TissueGrid () ;
    vector< vector<Grid> > grids ;
    ChemoBoundaryCondition chemoBoundary ;
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
    double grad_scale = 1;  //Controls steepnes of Gradient 
    vector<double> xSources ;
    vector<double> ySources ;
    vector<int> indexSourceX ;
    vector<int> indexSourceY ;
    vector<Grid> ghostBottom ;
    vector<Grid> ghostRight ;
    vector<Grid> ghostLeft ;
    vector<Grid> ghostTop ;
    
    int machineID = 1  ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    
    //---------------------------- Functions --------------------------------------------
    
    void DomainBoundaries (double xMin, double xMax, double yMin, double yMax, int nGridX , int nGridY) ;
    void InitializeAllGrids () ;
    void ClearChanges () ;
    void DiffusionChanges () ;
    void FindProductionPoints (vector<double> pSrc) ;
    void ProductionChanges () ;
    void DegredationChanges () ;
    void BoundaryChanges() ;
    void Update_GhostChanges () ;
    void UpdateChanges () ;
    void EulerMethod () ;
    void ParaViewGrids (int) ;
    void RoundToZero () ;
    void UpdateTGridsFolderNames (int ) ;
    void UpdateTGrid_FromConfigFile () ;
    
    
};
#endif /* TissueGrid_hpp */
