

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
    //Size of grif and domain for solving diffusion of the chemical
    double xDomainMin = 0 ;
    double xDomainMax = 100 ;
    double yDomainMin = 0 ;
    double yDomainMax = 100 ;
    int numberGridsX = 200 ;
    int numberGridsY = 200 ;
    double grid_dx ;
    double grid_dy ;
    
    //Parameters for solving diffusion equation using euler method
    double Diffusion = 200.0 ;
    double deg = 0.0001 ;
    double pro = pow(10, 2) ;
    double grid_dt = 0.05 ;      //timeStep
    double grad_scale = 1;  //Controls steepnes of Gradient
    int maxIterator = 100 ;
    
    //List of source's coordinates
    vector<double> xSources ;
    vector<double> ySources ;
    //List of grid index of the sources
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
    
    //Assign the size of the domain
    void DomainBoundaries (double xMin, double xMax, double yMin, double yMax, int nGridX , int nGridY) ;
    
    void InitializeAllGrids () ;
    void ClearChanges () ;
    //Calculate the changes due to diffusion
    void DiffusionChanges () ;
    // Find the grids that are meant to be a source
    void FindProductionPoints (vector<double> pSrc) ;
    //Calculate the changes due to production
    void ProductionChanges () ;
    //Calculate the changes due to degradation over time( propertion to the concentration)
    void DegredationChanges () ;
    //Calculate the changes due to the boundary effect
    void BoundaryChanges() ;
    //Calculate ghost changes, needed for boundary changes
    void Update_GhostChanges () ;
    //Apply all changes all at once
    void UpdateChanges () ;
    //Solve diffusion equation using Euler method.
    void EulerMethod () ;
    //Visualize the chemical concentration using the grids
    void ParaViewGrids (int) ;
    // Paraview can not read small numbers. This function round them to zaro
    void RoundToZero () ;
    //Update parameters based on input file
    void UpdateTGrid_FromConfigFile () ;
    
    
};
#endif /* TissueGrid_hpp */
