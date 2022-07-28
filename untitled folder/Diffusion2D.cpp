#include "Diffusion2D.hpp"
#include "TissueGrid.hpp"


//int main ()
TissueGrid Diffusion2D (double xMin, double xMax, double yMin, double yMax)
{
    
    TissueGrid tissue ;
    tissue.DomainBoundaries(xMin, xMax, yMin, yMax) ;
    tissue.xSources.push_back(tissue.xDomainMax / 2.0) ;
    tissue.ySources.push_back(tissue.yDomainMax / 3.0) ;
    tissue.xSources.push_back(tissue.xDomainMax / 2.0) ;
    tissue.ySources.push_back(tissue.yDomainMax / 3.0 * 2.0) ;
    tissue.InitializeAllGrids() ;
    tissue.FindProductionPoints() ;
    tissue.EulerMethod() ;

    
    return tissue ;
    //return 0 ;
    
    
}
