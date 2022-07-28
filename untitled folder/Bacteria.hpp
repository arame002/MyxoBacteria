
#ifndef Bacteria_hpp
#define Bacteria_hpp

#include "Nodes.hpp"

enum MovingDirection
{
    backward = 0 ,
    frward = 1
};

class bacterium
{   public:
    vector<node> nodes ;
    vector<node> ljnodes ;
    vector<node> allnodes ;
    vector<node> duplicate ;
    vector<double> connection ;
    double reversalPeriod = 300.0 ;
    double protein ;
    bool copy ;
    double reversalTime ;
    double fSTFx ;
    double fSTFy ;
    vector<pilus> pili ;
    double orientation ;
    
    bool turnStatus = false ;
    double turnTime = 0.0 ;
    double turnAngle = 0.0 ;
    vector<double> oldLoc ;
    double oldChem ;
    int numberReverse = 0 ;
    bool sourceWithin = false ;
    double timeToSource = 0.0 ;
    bool directionOfMotion = true ;        //true->forward
    
    bool SourceRegion () ;
    
    bacterium () ;
};


#endif /* Bacteria_hpp */
