//
//  TissueBacteria.cpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 10/8/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#include "TissueBacteria.hpp"

vector<vector<double> > TissueBacteria::Cal_Diffusion2D(double xMin, double xMax, double yMin, double yMax)
{
    tGrids = Diffusion2D(xMin, xMax, yMin, yMax ) ;
    vector<vector<double> > tmpGrid ;
    
    for (unsigned int i=0; i< tGrids.grids.size() ; i++)
    {
        vector<double> tmp ;
        for (unsigned int j=0; j< tGrids.grids.at(i).size(); j++)
        {
            tmp.push_back(tGrids.grids.at(i).at(j).value ) ;
        }
        tmpGrid.push_back(tmp) ;
        tmp.clear() ;
        
    }
    return tmpGrid ;
    
}

