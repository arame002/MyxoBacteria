
#ifndef Bacteria_hpp
#define Bacteria_hpp

#include "Nodes.hpp"

enum MovingDirection
{
    backward = 0 ,
    frward = 1
};
class MotilityMetabolism
{
public:
    double a0 = 0.33 ;
    double m0 = 1.0 ;
    double receptorActivity = a0 ;
    double lnA_a = -3.2 ;
    double b = 6.0 ;
    double gamma = 0.5 ;
    double switchProbability ;
    double legand = 20.0 ;              //This is going to be relaced by gridInMain
    double methylation  = 0.93 ;
    double changeRate = 0.0;
    double km = 1.7 ;
    double kI = 18.2 ;       //18.2 or 1.0
    double kA = 3000.0 ;       //3.0 or 1.0
    double N = 6 ;
    double kR = 0.005 ;
    double kB = 0.010 ;
    double MethylEnergy ;
    double LegandEnergy ;
    bool switchMode = false ;
    
    MotilityMetabolism ();
    double Cal_MethylationEnergy ();
    double Cal_LegandEnergy () ;
    double Cal_ReceptorActivity () ;
    double Cal_MethylationLevel (double tmpDt ) ;
    double Cal_SwitchProbability (double tmpDt ) ;
    void UpdateMotility_FromConfigFile() ;
    
    
};

class bacterium
{   public:
    vector<node> nodes ;
    vector<node> ljnodes ;
    vector<node> allnodes ;
    vector<node> duplicate ;
    vector<double> connection ;
    double reversalPeriod = 7.0 ;       //300       //This value updates during the simulations
    double maxRunDuration ;
    double protein ;
    bool copy ;
    double reversalTime ;
    double fSTFx ;
    double fSTFy ;
    vector<pilus> pili ;
    double orientation ;
    vector<vector<double> > chemoPointSources ;
    double locVelocity ;
    double locFriction ;
    
    
    bool turnStatus = false ;
    double turnTime = 0.0 ;
    double turnAngle = 0.0 ;        //the overall turn angle depends on fmotor, turnTime and chosen angle
    double maxTurnAngle = (3.1415)/5.0 ;
    double turnSDV = 3.1415/6.0 * (1.0/3.0) / ( (10.0 * 1.0) * 0.3 ) ;
    
    double chemotaxisPeriod = 0.1 ;             //How often the bacteria responds to the chemo-gradient
    double chemotaxisTime = 0.0 ;
    
    vector<double> oldLoc ;
    double oldChem = 0.0 ;
    int numberReverse = 0 ;
    bool sourceWithin = false ;
    double timeToSource = 0.0 ;
    double timeInSource = 0.0 ;
    bool directionOfMotion = true ;        //true->forward
    bool attachedToFungi ;
    bool SourceRegion () ;
    
    bool wrapMode = false ;
    double wrapDuration = 1.0 ;
    double wrapTime = 0.0 ;
    double wrapAngle = 0.0 ;
    double wrapRate = 1.0 ;
    double wrapAngularVelocity ;
    double wrapProbability = 1.0 ;
    double wrapSlowDown = 15.0/35.0 ;
    double maxWrapAngle = (3.1415)/6.0 ;
    double wrapSDV = (3.1415/2.0) * (1.0/3.0) / ( (10.0 * wrapDuration) * 0.3 ) ; //0.3 is fmotor
    double wrapMeanAngle = (3.1415/2.0)/ ( (10.0 * wrapDuration) * 0.3 ) ;      //0.3 is fmotor
    
    MotilityMetabolism motilityMetabolism ;
    void initialize_randomForce () ;
    void UpdateBacteria_FromConfigFile() ;
    double LogNormalMaxRunDuration (std::lognormal_distribution<double> &dist, std::default_random_engine &generator , double a, bool calib, double runVal) ;
    bacterium () ;
};


#endif /* Bacteria_hpp */
