
#include "SignalCommon.hpp"

using constants::pi;



double Dist2D (double x1 , double y1 , double x2 , double y2 )
{
    return sqrt(pow(x1- x2,2) + pow(y1- y2,2)) ;
}

//---------------------------------------------------------------------------------------------
vector<double> Dist_PointToVec2D (double x1, double y1, vector<double> X2, vector<double> Y2 )
{
    vector<double> point_vecDist ;
    for (unsigned int j =0 ; j < X2.size(); j++)
    {
        point_vecDist.push_back( Dist2D(x1, y1 ,X2.at(j), Y2.at(j)) ) ;
    }
    return point_vecDist ;
}

//---------------------------------------------------------------------------------------------
vector<vector <double> > Dist_VecToVec2D ( vector<double> X1 , vector<double> Y1 , vector<double> X2, vector<double> Y2 )
{
    vector<vector<double> > vec_vecDist ;
    for (unsigned int i=0 ; i < X1.size(); i ++)
    {
        vec_vecDist.push_back(Dist_PointToVec2D( X1.at(i) , Y1.at(i) ,X2 , Y2 ) ) ;
        
    }
    return vec_vecDist ;
}
//---------------------------------------------------------------------------------------------
vector<int> Indices_MinMatrix (vector<vector<double> > A )
{
    vector<int> indexRow ;
    vector<double> minRow ;
    int indexMin_Row ;
    int indexMin_Col ;
    vector<int> result ;
    for (unsigned int l = 0; l < A.size(); l++)
    {
        indexRow.push_back( static_cast<int>(min_element(A.at(l).begin(), A.at(l).end() ) -A.at(l).begin() ) ) ;
        minRow.push_back(*min_element(A.at(l).begin(), A.at(l).end())) ;
    }
    indexMin_Row = static_cast<int>(min_element(minRow.begin(), minRow.end()) - minRow.begin()) ;
    indexMin_Col = indexRow.at( indexMin_Row ) ;
    result.push_back(indexMin_Row) ;
    result.push_back(indexMin_Col) ;
    
    return result ;
}
//---------------------------------------------------------------------------------------------
vector<int> Indices_MaxMatrix (vector<vector<double> > A )
{
    vector<int> indexRow ;
    vector<double> maxRow ;
    int indexMax_Row ;
    int indexMax_Col ;
    vector<int> result ;
    for (unsigned int l = 0; l < A.size(); l++)
    {
        indexRow.push_back( static_cast<int>(max_element(A.at(l).begin(), A.at(l).end() ) -A.at(l).begin() ) ) ;
        maxRow.push_back(*max_element(A.at(l).begin(), A.at(l).end())) ;
    }
    indexMax_Row = static_cast<int>(max_element(maxRow.begin(), maxRow.end()) - maxRow.begin()) ;
    indexMax_Col = indexRow.at( indexMax_Row ) ;
    result.push_back(indexMax_Row) ;
    result.push_back(indexMax_Col) ;
    
    return result ;
}
//---------------------------------------------------------------------------------------------
vector<double> Dist_pointToVec1D (double x1, vector<double> X2 )
{
    vector<double> point_vecDist ;
    for (unsigned int j =0 ; j < X2.size(); j++)
    {
        point_vecDist.push_back( abs(X2.at(j)-x1) ) ;
    }
    return point_vecDist ;
}
//---------------------------------------------------------------------------------------------
double DotProduct (double x1 , double y1 , double x2 , double y2)
{
    return x1 * x2 + y1 * y2 ;
}
//---------------------------------------------------------------------------------------------
double MagnitudeVec (double x , double y)
{
    return sqrt(DotProduct(x, y, x, y)) ;
}
//---------------------------------------------------------------------------------------------
double TriangleArea (double x1, double y1, double x2 , double y2)
{
    double tetta1 = atan2(y1 , x1 ) ;
    double tetta2 = atan2(y2 , x2 ) ;
    double dTetta = fmod( abs( tetta2 -tetta1 ),2 * pi ) ;
    
    double area = 0.5 * MagnitudeVec(x1, y1 ) * MagnitudeVec(x2, y2 ) * abs( sin( dTetta ) ) ;
    return area ;
    
}
//---------------------------------------------------------------------------------------------
double AngleOfTwoVectors (double x1 , double y1 , double x2 , double y2)
{
    double Cos = DotProduct(x1, y1, x2, y2)/(MagnitudeVec(x1, y1) * MagnitudeVec(x2, y2)) ;
    return acos(Cos) ;
}

double sum_over_vec(const vector<vector<double> >& v, int a ) 
{
    return accumulate(v.begin(), v.end(), 0.0 ,
                      [&](double sum, vector<double> curr) { return sum + curr.at(a) ; });
};

double NormalCDFInverse2(double p) {
    
    if (p < 0.5) {
        return -RationalApproximation2( sqrt(-2.0*log(p)));
    }
    else {
        return RationalApproximation2( sqrt(-2.0*log(1-p)));
    }
}

double RationalApproximation2(double t) {
    
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return (t - ((c[2]*t + c[1])*t + c[0]) / (((d[2]*t + d[1])*t + d[0])*t + 1.0));
    
}

double AngleOfVector (double x1,double y1,double x2,double y2)
{
    double ax = x2 - x1 ;
    double ay = y2 - y1 ;
    double Cos = ax/sqrt(ax*ax+ay*ay) ;
    double angle = acos(Cos)* 180.0 / 3.1415 ;
    if(ay<0) angle *= -1 ;
    return angle ;                                   // (-180,180)
    
}



//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

