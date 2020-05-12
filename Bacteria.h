//
//  Bacteria.hpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 5/12/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//


#ifndef Bacteria_h
#define Bacteria_h

#include "Nodes.h"

class bacterium
{   public:
    node nodes[nnode];
    node ljnodes[nnode-1] ;
    node allnodes[2*nnode-1] ;
    double connection [nbacteria] ;
    double protein ;
    node duplicate[nnode] ;
    bool copy ;
    double reversalTime ;
    double fSTFx ;
    double fSTFy ;
    pilus pili[nPili] ;
    
    bool turnStatus = false ;
    double turnTime = 0.0 ;
    double turnAngle = 0.0 ;
};


#endif /* Bacteria_h */
