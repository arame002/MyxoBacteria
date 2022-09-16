

#ifndef Nodes_hpp
#define Nodes_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <random>
#include "SignalCommon.hpp"

using namespace std ;

#define nnode 5                                                // need to be an odd number
#define nbacteria 50                                            // 2* n^2
#define points nbacteria*nnode
//#define domainx  1000.0
//#define domainy  1000.0
#define nPili   1


class node
{   public:
    double x ;
    double y ;
    double fSpringx ;
    double fSpringy ;
    double fBendingx ;
    double fBendingy ;
    double fljx ;
    double fljy ;
    double xdev ;
    double ydev ;
    double fMotorx ;
    double fMotory ;
    double protein ;
    
    void test () ;
    
};

class pilus
{   public:
    double lFree ;                                // pili free length
    double retractionRate ;
    double vRet ;                                // retraction rate for connected pili
    double subAttachmentRate = 1 ;
    double subDetachmentRate = 0.1 ;
    double lContour ;
    double F ;                  //magnitude of the force
    double fx ;
    double fy ;
    double xEnd ;                   //x component of the end of the pili
    double yEnd ;                   //y component of the end of the pili
    bool piliSubs = true ;                 // true for pili-substrate connection and false for pili-pili connection
    // we change this parameter for simulation
    bool attachment ;
    bool retraction ;
    
};

#endif /* Nodes_hpp */

