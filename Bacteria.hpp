
#ifndef Bacteria_hpp
#define Bacteria_hpp

#include "Nodes.hpp"

class bacterium
{   public:
    vector<node> nodes ;
    vector<node> ljnodes ;
    vector<node> allnodes ;
    vector<node> duplicate ;
    vector<double> connection ;
    double protein ;
    bool copy ;
    double reversalTime ;
    double fSTFx ;
    double fSTFy ;
    vector<pilus> pili ;
    
    bool turnStatus = false ;
    double turnTime = 0.0 ;
    double turnAngle = 0.0 ;
    
    
    bacterium () ;
};


#endif /* Bacteria_hpp */
