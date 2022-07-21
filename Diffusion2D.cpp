#include "Diffusion2D.hpp"



TissueGrid Diffusion2D (double xMin, double xMax, double yMin, double yMax,int nGridX , int nGridY ,vector<vector<double> > sources, vector<double> pSrc, TissueGrid tissue)
{
    vector<vector<double> > pointSource = sources ;
    tissue.DomainBoundaries(xMin, xMax, yMin, yMax, nGridX , nGridY) ;
    /*
    tissue.xSources.push_back(tissue.xDomainMax / 2.0) ;
    tissue.ySources.push_back(tissue.yDomainMax / 3.0) ;
    tissue.xSources.push_back(tissue.xDomainMax / 2.0) ;
    tissue.ySources.push_back(tissue.yDomainMax / 3.0 * 2.0) ;
     */
    tissue.xSources = pointSource.at(0) ;
    tissue.ySources = pointSource.at(1) ;
    tissue.InitializeAllGrids() ;
    tissue.FindProductionPoints(pSrc) ;
    tissue.EulerMethod() ;

    
    return tissue ;
    
    
}
