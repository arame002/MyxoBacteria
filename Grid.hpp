

#ifndef Grid_hpp
#define Grid_hpp

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <sstream>
using namespace std ;
#define pi 3.1415

class Grid
{
public:
    double value = 0.0 ;
    double change = 0.0 ;
    double productionRate = 0.0 ;
    
    
    void UpdateValue () ;
    
};
#endif /* Grid_hpp */
