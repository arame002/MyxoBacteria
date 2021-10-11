#ifndef HYPHAESEGMENT_H_
#define HYPHAESEGMENT_H_

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class HyphaeSegment {
    public:
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
    
    
    
};

#endif
