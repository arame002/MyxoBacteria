
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
#include "ConfigParser.h"
#include "constants.h"
using namespace std ;


double Dist2D (double x1 , double y1 , double x2 , double y2 ) ;
vector<double> Dist_PointToVec2D ( double x1, double y1, vector<double> X2, vector<double> Y2 ) ;
vector<vector <double> > Dist_VecToVec2D ( vector<double> , vector<double> , vector<double>, vector<double>  ) ;
vector<int> Indices_MinMatrix (vector<vector<double > >) ;
vector<int> Indices_MaxMatrix (vector<vector<double> >  ) ;
vector<double> Dist_pointToVec1D ( double x1, vector<double> X2) ;

double DotProduct (double x1 , double y1 , double x2 , double y2) ;
double AngleOfTwoVectors (double x1 , double y1 , double x2 , double y2) ;
double MagnitudeVec (double x , double y) ;
double TriangleArea (double x1, double y1, double x2 , double y2) ;
double sum_over_vec(const vector<vector<double> >& v, int a ) ;
double AngleOfVector (double x1, double y1 , double x2 , double y2 ) ;      // (x1,y1,x2,y2)
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

class plusNum
{
private:
    double num  ;
public:
    plusNum(double n) : num(n) {}
    
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num) const
    {
        return num + arr_num;
    }
} ;
//---------------------------------------------------------------------------------------------
class productNum
{
private:
    double num ;
public:
    productNum(double n) : num(n) {}
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num) const
    {
        return num * arr_num;
    }
} ;

class linearConfig
{
private:
    double a ;
    double b ;
public:
    linearConfig(double n , double m )
    {
        a = n ;
        b = m ;
    }
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num1 , double arr_num2) const
    {
        return a * arr_num1 + b * arr_num2 ;
    }
} ;

class productVec
{
private:
    
public:
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num1 , double arr_num2) const
    {
        return  arr_num1 * arr_num2 ;
    }
} ;

double NormalCDFInverse2(double p);
double RationalApproximation2(double t);


