
#ifndef TissueBacteria_hpp
#define TissueBacteria_hpp

#include <stdio.h>
#include "Diffusion2D.hpp"
#include "Bacteria.hpp"

#endif /* TissueBacteria_hpp */
class TissueBacteria
{
public:
    bacterium bacteria[nbacteria];
    TissueGrid tGrids ;
    
    vector<vector<double> > Cal_Diffusion2D (double xMin, double xMax, double yMin, double yMax) ;
    
    
};
