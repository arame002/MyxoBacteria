

#ifndef Diffusion2D_h
#define Diffusion2D_h

#include "TissueGrid.hpp"

TissueGrid Diffusion2D (double xMin, double xMax, double yMin, double yMax ,int nGridX , int nGridY, vector<vector<double> > tips, vector<double> pSrc, TissueGrid tissue) ;

#endif /* Diffusion2D_h */
