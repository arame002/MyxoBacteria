//
//  Fungi.hpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/19/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#ifndef Fungi_hpp
#define Fungi_hpp


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "growthFunctions.h"
#include "hyphaeSegment.h"
class Fungi {
public:
    
    vector<HyphaeSegment> hyphaeSegments ;
    vector<vector<double> > tips ;
    vector<int > tipsID ;
    vector<int> centerConnectedSegments ;
    double hyphaeWidth = 6.0 ;
    double production = pow(10, 2) ;
    double proDecayFactor = 0.5 ;
    
    void Find_Hyphae_Tips() ;
    void Find_Hyphae_Tips2() ;
    void WriteSourceLoc(vector<vector<double> >) ;
    void FindProductionNetwork ();
    void FindProductionNetwork2();
    vector<int> FindNeighborList (int i) ;
    void Rec_CalNeighborProduction (int id, vector<int> ngbrList, int nbrElement, bool connection) ;
    void Rec_CalNeighborProduction2(int id, vector<int> ngbrList, int nbrElement ) ;
    vector<int> FindCenterPointConnection () ;
    
};


#endif /* Fungi_hpp */

