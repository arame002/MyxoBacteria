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
    
    void Find_Hyphae_Tips() ;
    
};


#endif /* Fungi_hpp */

