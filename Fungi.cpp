//
//  Fungi.cpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/19/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#include "Fungi.hpp"
void Fungi::Find_Hyphae_Tips()
{
    tips.resize(2) ;
    for (uint i=0 ; i<hyphaeSegments.size() ; i++)
    {
        
        if (hyphaeSegments[i].can_extend == true)
        {
            /*
            for (int j=0; j<4; j++)
            {
            
                tips[0].push_back(hyphaeSegments[i].x2- j*(hyphaeSegments[i].x2-hyphaeSegments[i].x1)/3.0) ;
                tips[1].push_back(hyphaeSegments[i].y2- j*( hyphaeSegments[i].y2-hyphaeSegments[i].y1)/3.0) ;
            }
             */
             // tips only
            tips[0].push_back(hyphaeSegments[i].x2) ;
            tips[1].push_back(hyphaeSegments[i].y2) ;
            tipsID.push_back(i);
             
        }
         
        
    }
    for (uint i=0 ; i<hyphaeSegments.size() ; i++)
    {
        if (hyphaeSegments[i].can_branch == false)
        {
            tips[0].push_back(hyphaeSegments[i].x2) ;
            tips[1].push_back(hyphaeSegments[i].y2) ;
            tipsID.push_back(i);
        }
    }
}

void Fungi:: WriteSourceLoc(vector<vector<double> > pointSource)
{
    ofstream sourceLoc ("sourceLocation.txt") ;
    for (uint i = 0; i< pointSource.at(0).size() ; i++)
    {
        sourceLoc << pointSource.at(0).at(i) <<'\t'<< pointSource.at(1).at(i) <<endl ;
    }
    sourceLoc<< endl ;
}

