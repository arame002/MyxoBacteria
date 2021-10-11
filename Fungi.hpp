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
    double hyphaeWidth = 6.0 ;
    
    void Find_Hyphae_Tips() ;
    void WriteSourceLoc() ;
    
};


#endif /* Fungi_hpp */

