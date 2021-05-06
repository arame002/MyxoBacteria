//
//  TissueBacteria.cpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 10/8/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#include "TissueBacteria.hpp"

vector<vector<double> > TissueBacteria::Cal_Diffusion2D(double xMin, double xMax, double yMin, double yMax,vector<vector<double> > sources)
{
    tGrids = Diffusion2D(xMin, xMax, yMin, yMax, sources ) ;
    vector<vector<double> > tmpGrid ;
    
    for (unsigned int i=0; i< tGrids.grids.size() ; i++)
    {
        vector<double> tmp ;
        for (unsigned int j=0; j< tGrids.grids.at(i).size(); j++)
        {
            tmp.push_back(tGrids.grids.at(i).at(j).value ) ;
        }
        tmpGrid.push_back(tmp) ;
        tmp.clear() ;
        
    }
    return tmpGrid ;
    
}



void TissueBacteria:: Myxo ()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for ( int j=0 ; j<nnode ; j++)
        {
            //     bacteria[i].nodes[j].x = i* nnode + j ;
            //   bacteria[i].nodes[j].y = 2*i*nnode + 2* j  ;
            
            bacteria[i].nodes[j].x = 20.0 + j ;
            bacteria[i].nodes[j].y = 20.0  ;
            
        }
    }
    
}


void TissueBacteria:: Spring()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for (int j=0; j<nnode; j++)
        {
            if (j==0 )
            {
                bacteria[i].nodes[j].fSpringx= -K * (Distance (i,j,i,j+1 ) - x0) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x) / Distance (i,j,i,j+1) ;
                
                bacteria[i].nodes[j].fSpringy= -K * (Distance (i,j,i,j+1 ) - x0) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y) / Distance (i,j,i,j+1) ;
            }
            
            else if (j== (nnode-1))
            {
                bacteria[i].nodes[j].fSpringx = -K * (Distance (i,j,i,j-1) - x0) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j-1].x) / Distance (i,j,i,j-1) ;
                
                bacteria[i].nodes[j].fSpringy = -K * (Distance (i,j,i,j-1) - x0) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j-1].y) / Distance (i,j,i,j-1) ;
            }
            else
            {
                bacteria[i].nodes[j].fSpringx  = -K * ( Distance(i,j,i,j-1) - x0 ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j-1].x) /Distance(i,j,i,j-1) ;
                bacteria[i].nodes[j].fSpringx += -K * ( Distance(i,j,i,j+1) - x0 ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x) /Distance(i,j,i,j+1) ;
                
                bacteria[i].nodes[j].fSpringy  = -K * ( Distance(i,j,i,j-1) - x0 ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j-1].y) /Distance(i,j,i,j-1) ;
                bacteria[i].nodes[j].fSpringy += -K * ( Distance(i,j,i,j+1) - x0 ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y) /Distance(i,j,i,j+1) ;
            }
        }
    }
    
}

double TissueBacteria:: Distance (int i,int j, int k, int l)
{
    double dis= sqrt ((bacteria[i].nodes[j].x - bacteria[k].nodes[l].x) * (bacteria[i].nodes[j].x-bacteria[k].nodes[l].x) + (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y) * (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y)) ;
    return dis ;
}

void TissueBacteria:: Bending()
{
    double Bend1x[nnode] ;
    double Bend2x[nnode] ;
    double Bend3x[nnode] ;
    
    double Bend1y[nnode] ;
    double Bend2y[nnode] ;
    double Bend3y[nnode] ;
    
    for (int i=0 ; i<nbacteria; i++)
    {   for(int j=0; j<nnode ; j++)
    {
        bacteria[i].nodes[j].fBendingx=0.0 ;
        bacteria[i].nodes[j].fBendingy=0.0 ;
    }
    }
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<nnode; j++)
        {
            Bend1x[j]=0.0 ;
            Bend2x[j]=0.0 ;
            Bend3x[j]=0.0 ;
            
            Bend1y[j]=0.0 ;
            Bend2y[j]=0.0 ;
            Bend3y[j]=0.0 ;
        }
        
        for (int j=0 ;j<nnode-2; j++)
        {
            Bend1x[j] += -B *(bacteria[i].nodes[j+2].x-bacteria[i].nodes[j+1].x)/(Distance(i,j,i,j+1) * Distance(i,j+2,i,j+1)) ;
            Bend1x[j] +=  B *(Cos0ijk(i,j+1)/(Distance(i,j,i,j+1)*Distance(i,j,i,j+1))*(bacteria[i].nodes[j].x-bacteria[i].nodes[j+1].x)) ;
            
            Bend1y[j] += -B *(bacteria[i].nodes[j+2].y-bacteria[i].nodes[j+1].y)/(Distance(i,j,i,j+1) * Distance(i,j+2,i,j+1)) ;
            Bend1y[j] +=  B *(Cos0ijk(i,j+1)/(Distance(i,j,i,j+1)*Distance(i,j,i,j+1))*(bacteria[i].nodes[j].y-bacteria[i].nodes[j+1].y)) ;
            
        }
        for (int j=2; j<nnode ; j++)
        {
            Bend3x[j] += -B *(bacteria[i].nodes[j-2].x-bacteria[i].nodes[j-1].x)/(Distance(i,j-2,i,j-1) * Distance(i,j,i,j-1)) ;
            Bend3x[j] +=  B *(Cos0ijk(i,j-1)/(Distance(i,j,i,j-1)*Distance(i,j,i,j-1))*(bacteria[i].nodes[j].x-bacteria[i].nodes[j-1].x)) ;
            Bend3y[j] += -B *(bacteria[i].nodes[j-2].y-bacteria[i].nodes[j-1].y)/(Distance(i,j-2,i,j-1) * Distance(i,j,i,j-1)) ;
            Bend3y[j] +=  B *(Cos0ijk(i,j-1)/(Distance(i,j,i,j-1)*Distance(i,j,i,j-1))*(bacteria[i].nodes[j].y-bacteria[i].nodes[j-1].y)) ;
        }
        for (int j=1; j<nnode-1; j++)
        {
            Bend2x[j] += -1 * (Bend1x[j-1] + Bend3x[j+1] ) ;
            Bend2y[j] += -1 * (Bend1y[j-1] + Bend3y[j+1] ) ;
        }
        
        for (int j=0; j<nnode; j++)
        {
            bacteria[i].nodes[j].fBendingx=  Bend1x[j] + Bend2x[j] + Bend3x[j] ;
            bacteria[i].nodes[j].fBendingy=  Bend1y[j] + Bend2y[j] + Bend3y[j] ;
        }
    }
    
}

double TissueBacteria:: Cos0ijk(int i,int m)
{
    double dot= (bacteria[i].nodes[m-1].x - bacteria[i].nodes[m].x)*(bacteria[i].nodes[m+1].x - bacteria[i].nodes[m].x);
    dot += (bacteria[i].nodes[m-1].y - bacteria[i].nodes[m].y)*(bacteria[i].nodes[m+1].y -bacteria[i].nodes[m].y) ;
    double cos= dot / ( Distance(i,m-1,i,m) * Distance(i,m+1,i,m) ) ;
    return cos ;
}

double TissueBacteria:: Distance2 (double x1, double y1, double x2, double y2, int nx , int ny)
{
    double dis= sqrt((x1+nx*domainx -x2)*(x1+nx*domainx -x2)+(y1+ny*domainy -y2)*(y1+ny*domainy -y2)) ;
    return dis ;
}


double TissueBacteria::MinDistance (double x1 , double y1 ,double x2 , double y2)
{   shiftx = 0 ;
    shifty = 0 ;
    double min = Distance2(x1, y1, x2, y2 ,0 , 0) ;
    double a ;
    
    for (int nx = -1 ; nx<2; nx++)
    {
        for (int ny= -1 ; ny<2; ny++)
        {   a = Distance2(x1, y1, x2, y2,nx,ny) ;
            if (a< min)
            {
                min= a ;
                shiftx = nx ;
                shifty = ny ;
            }
        }
    }
    return min ;
}


void TissueBacteria::Connection ()
{
    double xmin = sqrt((x0/4)*(x0/4)+ 0.36) ;
    double d  ;
    double min ;
    for (int i=0; i<nbacteria-1; i++)
    {
        for (int j=i+1; j<nbacteria; j++)
        {
            bacteria[i].connection[j] = 0.0 ;
            bacteria[j].connection[i] = 0.0 ;
            
            min = MinDistance(bacteria[i].nodes[(nnode-1)/2].x , bacteria[i].nodes[(nnode-1)/2].y , bacteria[j].nodes[(nnode-1)/2].x , bacteria[j].nodes[(nnode-1)/2].y)  ;
            
            if (min < (1.2 * length) )
            {
                for (int m=0; m<(2*nnode-1); m++)
                {
                    for (int n=0; n<(2*nnode-1); n++)
                    {   d = Distance2 ( bacteria[i].allnodes[m].x , bacteria[i].allnodes[m].y, bacteria[j].allnodes[n].x , bacteria[j].allnodes[n].y ,shiftx , shifty ) ;
                        if ( d < xmin )
                        { //  bacteria[i].connection[j] += 0.5 ;
                            //  bacteria[j].connection[i] += 0.5 ;
                            if (m==0 || n==0 || m== (2*nnode-2)|| n== (2*nnode-2) )
                            {
                                bacteria[i].connection[j] += 0.25 ;
                                bacteria[j].connection[i] += 0.25 ;
                            }
                            else
                            {
                                bacteria[i].connection[j] += 0.5 ;
                                bacteria[j].connection[i] += 0.5 ;
                            }
                        }
                        
                        
                    }
                }
                
            }
        }
    }
}

double TissueBacteria:: u_lj()
{
    double delta_x, delta_y, rcut2;
    for (int i=0; i<nbacteria; i++)
    {   for(int j=0 ; j< nnode ; j++)
    {
        bacteria[i].nodes[j].fljx=0.0;
        bacteria[i].nodes[j].fljy=0.0;
    }
    }
    double rMin = 0.6 ;
    double sigma2= (rMin/1.122)*(rMin/1.222) ;          // Rmin = 1.122 * sigma
    double eps= 0.01 ;
    double ulj=0.0 ;
    rcut2= 6.25 * sigma2  ;                 // rcut=2.5*sigma for Lenard-Jones potential
    double min ;
    for (int i=0; i<(nbacteria -1); i++)
    {
        for (int j=i+1 ; j<nbacteria ; j++)
        {  min = MinDistance(bacteria[i].nodes[(nnode-1)/2].x , bacteria[i].nodes[(nnode-1)/2].y , bacteria[j].nodes[(nnode-1)/2].x , bacteria[j].nodes[(nnode-1)/2].y)  ;
            
            if ( min < (1.2 * length))
            {
                for (int m=0; m<2*nnode-1; m++)
                {
                    for (int n=0 ; n<2*nnode-1 ; n++)
                    {
                        if (m%2==1 && n%2==1) continue ;
                        delta_x=bacteria[i].allnodes[m].x-bacteria[j].allnodes[n].x + ( shiftx * domainx) ;
                        delta_y=bacteria[i].allnodes[m].y-bacteria[j].allnodes[n].y + ( shifty * domainy)  ;
                        double delta2= (delta_x * delta_x) + (delta_y * delta_y) ;
                        if (delta2<rcut2)              //ignore particles having long distance
                        {
                            double r2i=sigma2/delta2 , r6i=r2i*r2i*r2i ;
                            double force=48*eps*r6i*(r6i-0.5)*r2i;
                            ulj=ulj+4*eps*r6i*(r6i-1) ;
                            if (m%2==0)
                            {
                                bacteria[i].nodes[m/2].fljx += force*delta_x ;
                                bacteria[i].nodes[m/2].fljy += force*delta_y;
                            }
                            if (n%2 == 0)
                            {
                                bacteria[j].nodes[n/2].fljx -= force*delta_x ;
                                bacteria[j].nodes[n/2].fljy -= force*delta_y ;
                            }
                        }
                    }
                }
            }
        }
    }
    return ulj ;
}
/*
void TissueBacteria:: PiliForce ()
{
    
    
    double orientationBacteria ;
    double alfa ;                           // dummy name for calculating different angles
    nAttachedPili = 0 ;
    
    for (int i=0 ; i<nbacteria ; i++)
    {
        
        orientationBacteria = AngleOfVector(bacteria[i].nodes[1].x , bacteria[i].nodes[1].y , bacteria[i].nodes[0].x, bacteria[i].nodes[0].y);        //degree
        for ( int j=0 ; j<nPili ; j++)
        {
            if( bacteria[i].pili[j].attachment== false)     //if the pili is detached
            {
                //lambda substrate attachment
                if (rand() / (RAND_MAX + 1.0) < bacteria[i].pili[j].subAttachmentRate* dt)
                {
                    bacteria[i].pili[j].attachment = true ;
                    bacteria[i].pili[j].retraction = true ;
                    nAttachedPili += 1 ;
                    //now it is attached, then we calculate the end point
                    
                    if (bacteria[i].pili[j].piliSubs)       //pili-substrate interaction
                    {
                        alfa = orientationBacteria + 180.0/nPili * (j+1/2.0) -90.0 ;    //angle in degree
                        bacteria[i].pili[j].xEnd = fmod ( bacteria[i].nodes[0].x + bacteria[i].pili[j].lFree * cos(alfa*tetta0/180.0) + domainx ,domainx ) ;
                        bacteria[i].pili[j].yEnd = fmod ( bacteria[i].nodes[0].y + bacteria[i].pili[j].lFree * sin(alfa*tetta0/180.0) + domainy, domainy ) ;
                    }
                    else            //pili-pili interaction
                    {
                        //calculate the intersection interval and intersection point
                    }
                    
                }
            }
            
            else                    //if the pili is attached
            {
                //lambda detachment
                if ( rand() / (RAND_MAX + 1.0) < bacteria[i].pili[j].subDetachmentRate* dt )
                {
                    bacteria[i].pili[j].attachment = false ;
                    bacteria[i].pili[j].fx = 0.0 ;
                    bacteria[i].pili[j].fy = 0.0 ;
                    bacteria[i].pili[j].F = 0.0 ;
                    
                }
                else        //if the pili remains attached
                {
                    nAttachedPili += 1 ;
                    if (bacteria[i].pili[j].piliSubs)
                    {
                        // calculate forces for pili-substrate
                        bacteria[i].pili[j].lContour = MinDistance(bacteria[i].nodes[0].x, bacteria[i].nodes[0].y, bacteria[i].pili[j].xEnd, bacteria[i].pili[j].yEnd ) ;
                        
                        alfa = AngleOfVector(bacteria[i].nodes[0].x + shiftx * domainx , bacteria[i].nodes[0].y + shifty * domainy , bacteria[i].pili[j].xEnd,bacteria[i].pili[j].yEnd) ;
                        
                        bacteria[i].pili[j].F = max(0.0 , kPullPili * (bacteria[i].pili[j].lContour-bacteria[i].pili[j].lFree) ) ;
                        
                        bacteria[i].pili[j].fx = bacteria[i].pili[j].F * cos(alfa*tetta0/180.0)  ;
                        
                        bacteria[i].pili[j].fy = bacteria[i].pili[j].F * sin(alfa*tetta0/180.0)  ;
                        
                    }
                    
                    else
                    {
                        //calculate forces for pili pili
                    }
                }
                
            }
        }
    }
    averageLengthFree = 0 ;
    for (int i=0 ; i<nbacteria ; i++)
    {
        for ( int j=0 ; j<nPili ; j++)
        {
            if(bacteria[i].pili[j].attachment)
            {
                bacteria[i].pili[j].retraction = true ;
            }
            //if the bacteria is not attached, switch state between Protrusion and Retraction
            else if (rand() / (RAND_MAX + 1.0) < bacteria[i].pili[j].retractionRate * dt)
            {
                bacteria[i].pili[j].retraction = ! bacteria[i].pili[j].retraction ;
            }
            
            if (bacteria[i].pili[j].retraction )
            {
                bacteria[i].pili[j].vRet = vRetPili0*(1-bacteria[i].pili[j].F / fStall)  ;
                if (bacteria[i].pili[j].vRet < 0.0)
                {
                    bacteria[i].pili[j].vRet = 0.0 ;
                }
                bacteria[i].pili[j].lFree += -bacteria[i].pili[j].vRet * dt ;
                if (bacteria[i].pili[j].lFree< 0.0)
                {
                    bacteria[i].pili[j].lFree = 0.0 ;
                    bacteria[i].pili[j].attachment = false ;      //the pili is detached when its length is zero
                    bacteria[i].pili[j].retraction = false ;
                    bacteria[i].pili[j].fx = 0.0 ;
                    bacteria[i].pili[j].fy = 0.0 ;
                    bacteria[i].pili[j].F  = 0.0 ;
                }
            }
            else
            {
                bacteria[i].pili[j].lFree += vProPili * dt ;
                if ( bacteria[i].pili[j].lFree > piliMaxLength )     //pili length is limited
                {
                    bacteria[i].pili[j].lFree = piliMaxLength ;
                }
            }
            averageLengthFree += bacteria[i].pili[j].lFree / (nbacteria * nPili) ;
        }
    }
}
*/
/*
void TissueBacteria:: RandomForce()
{
    double sig ;
    for(int i=0; i<nbacteria ; i++)
    {
        for(int j=0 ; j<nnode; j++)
        {
            Diffusion(bacteria[i].nodes[j].x, bacteria[i].nodes[j].y) ;
            sig = sqrt(2.0*diffusion*dt) ;
            bacteria[i].nodes[j].xdev =gasdev()*sig;
            bacteria[i].nodes[j].ydev =gasdev()*sig;
        }
    }
}
*/
void TissueBacteria:: Diffusion (double x , double y)
{
    //int m = (static_cast<int> (round ( fmod (x + domainx , domainx) / dx ) ) ) % nx ;
    //  int n = (static_cast<int> (round ( fmod (y + domainy , domainy) / dy ) ) ) % ny ;
    //  diffusion = kblz * temp/ (eta1/ slime[m][n] ) ;
    diffusion = kblz * temp/ eta1 ;
}

