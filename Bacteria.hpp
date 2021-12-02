
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
    double reversalPeriod = 7.0 ;       //300
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
    double maxTurnAngle = (3.1415)/5.0 ;
    vector<double> oldLoc ;
    double oldChem = 0.0 ;
    int numberReverse = 0 ;
    bool sourceWithin = false ;
    double timeToSource = 0.0 ;
    bool directionOfMotion = true ;        //true->forward
    bool attachedToFungi ;
    bool SourceRegion () ;
    
    bacterium () ;
};


#endif /* Bacteria_hpp */
