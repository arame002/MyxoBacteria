#ifndef HYPHAESEGMENT_H_
#define HYPHAESEGMENT_H_

#include <iostream>
#include <string>
#include <vector>
#include "ConfigParser.h"

using namespace std;

class HyphaeSegment {
    public:
        
    double hyphaeLength = 40 ;      //80
    double x1;
    double y1;
    double angle;
    double x2;
    double y2;
    bool can_branch;
    bool can_extend;
    int from_who;
    int extend_to;
    int branch_to;
    int hyphae_ID;
    
    double p1 = 0.0 ;     //production at the beginning
    double p2 = 0.0 ;     //production at the end
    double tmpP1 = 0.0 ;
    double tmpP2 = 0.0 ;
    HyphaeSegment () ;
    
    
    
};

#endif
