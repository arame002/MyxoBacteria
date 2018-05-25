
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "ranum2.h"
#include <fstream>
using namespace std ;
/*
 namespace patch{
 template <typename  T> string to_string (const T& n)
 {
 ostringstream stm ;
 stm << n ;
 return stm.str() ;
 }
 }
 */

#define nnode 7                                                // need to be an odd number
#define nbacteria 18                                           // 2* n^2
#define points nbacteria*nnode
#define domainx  40.0
#define domainy  40.0
double initialTime = 4.0 ;
double runTime = 1000.0 ;


long  idum=(-799);

void Spring() ;
void Bending() ;
void Print() ;
double Distance (int ,int, int, int) ;                                     // bacteria, node, bacteria, node
double Distance2 (double, double, double, double, int, int ) ;             //x1,y1,x2,y2, shiftx, shifty
double MinDistance (double , double , double , double) ;                             // x1,y1,x2,y2
double Cos0ijk(int, int) ;
double u_lj () ;
void RandomForce() ;
void Motor () ;
void Reverse () ;                                               // inverse head and tail
void Initialization () ;
void PositionUpdating (double ) ;                               // (time step)
void ljNodesPosition () ;
void InitialProtein () ;
void Connection () ;
void AllNodes () ;
void ProteinExchange () ;
void NodeProtein () ;
void Duplicate () ;
void Merge () ;
void ParaView ( ) ;                                    // (current time, time step)
void ParaView2 () ;
void Track () ;
void Diffusion (double, double ) ;
//-----------------------------------------------------------------------------------------------------
double length = 5.0 ;
double B = 1.0 ;                      // bending constant
double tetta0=3.1415 ;               // prefered angle
double K =  1.0 ;                      // linear spring constant (stiffness)
double x0 = length/(nnode-1) ;                      // equilibrium distance of spring
double kblz=1.0 , temp=0.01 , eta1 = 1.0 , eta2 = 5.0 ;
double diffusion ;
double dt=0.0001 ;
double initialStep = 0.0001 ;
double fmotor = 0.1 ;
double Rp = 0.3 ;
double reversalTime = 1500.0 ;
double sr = 0.25 ;                          // slime rate production
double kd = 0.1 ;
int shiftx = 0 ;
int shifty = 0 ;
int index1 = 0 ;                            // needed for ParaView



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
    
};

class bacterium
{   public:
    node nodes[nnode];
    node ljnodes[nnode-1] ;
    node allnodes[2*nnode-1] ;
    double connection [nbacteria] ;
    double protein ;
    node duplicate[nnode] ;
    bool copy ;
    
};

bacterium bacteria[nbacteria];
const double dx= 0.5 ;
const double dy= 0.5 ;
const int nx = static_cast<int>(domainx/dx) ;
const int ny = static_cast<int> (domainy/dy) ;
const int nz = 1 ;
double X[nx] ;
double Y[ny] ;
double Z[nz]= {0.0} ;
double slime[nx][ny] ;



int main ()
{
    for (int m=0; m<nx; m++)
    {
        for (int n=0; n<ny; n++)
        {
            slime[m][n] = 1.0  ;
        }
    }
    
    for (int i=0; i<nx; i++)
    {
        X[i]= i*dx ;
    }
    for (int j=0; j<ny; j++)
    {
        Y[j]= j * dy ;
    }
    
   
//-----------------------------------------------------------------------------------------------------
    ofstream ProteinLevelFile ;
    ProteinLevelFile.open("ProteinLevelFile.txt") ;
    
 //   cout<<x0<<endl ;
    cout<<"program is running"<<endl  ;
    double nt=runTime/dt + initialTime/initialStep ;                   //total number of steps
    nt =static_cast<int>(nt) ;
    double initialNt =initialTime/initialStep ;                     // run time for initialization
    initialNt =static_cast<int>(initialNt) ;
    int inverseInitialStep = static_cast<int>(initialNt/initialTime) ;  // used for visualization(ParaView)
    int inverseDt =static_cast<int>((nt-initialNt)/(runTime)) ;              // used for visualization(ParaView)
    
    
    int  revnt = static_cast<int>(reversalTime/ dt) ;
    Initialization() ;
    ljNodesPosition() ;
    InitialProtein() ;
    for (int i=0; i<nbacteria; i++)
    {
        bacteria[i].connection[i]=0.0 ;
    }
    


 for (int l=0; l< (nt+1); l++)
    {
        if (l < initialNt)
        {
     
            AllNodes() ;
            Spring () ;
            Bending() ;
            u_lj() ;
            RandomForce() ;
     /*
           if (l%inverseInitialStep==0)
           {
               ParaView() ;
               ParaView2() ;
           }
      */
            PositionUpdating(initialStep) ;
        
        }
        else
        {
            if ( (l % revnt ==0) & (l / revnt != 0))
            {
                Reverse () ;
            }
       
            AllNodes() ;
            Spring () ;
           Bending() ;
            u_lj() ;
            RandomForce() ;
            Motor () ;
            Track() ;
            Connection() ;
            ProteinExchange() ;
            
            
            if (l%inverseDt==0)
            {
       
                for (int i=0; i<nbacteria; i++)
                {
                    ProteinLevelFile<<bacteria[i].protein<<"," ;
                }
                ProteinLevelFile<< endl<<endl ;
         
                ParaView ()  ;
                ParaView2() ;
                cout<<(l-initialNt)/inverseDt<<endl ;
                
            }
            PositionUpdating(dt) ;
        }
        
    }

    cout<<"No Bug"<<endl ;
    return 0 ;
}


/*_____________________________________________________________________________
 ______________________________________________________________________________*/

void PositionUpdating (double t)
{
    for (int i=0; i<nbacteria; i++)
    {
    for(int j=0 ; j<nnode; j++)
    {   Diffusion(bacteria[i].nodes[j].x , bacteria[i].nodes[j].y) ;
        bacteria[i].nodes[j].x += t * (bacteria[i].nodes[j].fSpringx + bacteria[i].nodes[j].fBendingx ) * diffusion / (kblz * temp) ;
        bacteria[i].nodes[j].x += bacteria[i].nodes[j].xdev ;
        bacteria[i].nodes[j].x += t* (bacteria[i].nodes[j].fljx + bacteria[i].nodes[j].fMotorx) * diffusion / (kblz*temp) ;
        
        bacteria[i].nodes[j].y += t * (bacteria[i].nodes[j].fSpringy + bacteria[i].nodes[j].fBendingy ) * diffusion / (kblz * temp) ;
        bacteria[i].nodes[j].y += bacteria[i].nodes[j].ydev ;
        
        bacteria[i].nodes[j].y += t * (bacteria[i].nodes[j].fljy + bacteria[i].nodes[j].fMotory) * diffusion /(kblz*temp);
        
    }
    }
    ljNodesPosition() ;
}

//-----------------------------------------------------------------------------------------------------

void Spring()
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
//-----------------------------------------------------------------------------------------------------
double Distance (int i,int j, int k, int l)
{
    double dis= sqrt ((bacteria[i].nodes[j].x - bacteria[k].nodes[l].x) * (bacteria[i].nodes[j].x-bacteria[k].nodes[l].x) + (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y) * (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y)) ;
    return dis ;
}
//-----------------------------------------------------------------------------------------------------
double Distance2 (double x1, double y1, double x2, double y2, int nx , int ny)
{
    double dis= sqrt((x1+nx*domainx -x2)*(x1+nx*domainx -x2)+(y1+ny*domainy -y2)*(y1+ny*domainy -y2)) ;
    return dis ;
}
//-----------------------------------------------------------------------------------------------------

double Cos0ijk(int i,int m)
{
    double dot= (bacteria[i].nodes[m-1].x - bacteria[i].nodes[m].x)*(bacteria[i].nodes[m+1].x - bacteria[i].nodes[m].x);
    dot += (bacteria[i].nodes[m-1].y - bacteria[i].nodes[m].y)*(bacteria[i].nodes[m+1].y -bacteria[i].nodes[m].y) ;
    double cos= dot / ( Distance(i,m-1,i,m) * Distance(i,m+1,i,m) ) ;
    return cos ;
}
//-----------------------------------------------------------------------------------------------------
void Bending()
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

//-----------------------------------------------------------------------------------------------------

double u_lj()
{
    double delta_x, delta_y, rcut2;
    for (int i=0; i<nbacteria; i++)
    {   for(int j=0 ; j< nnode ; j++)
    {
        bacteria[i].nodes[j].fljx=0.0;
        bacteria[i].nodes[j].fljy=0.0;
    }
    }
    double sigma2= (0.5/1.122)*(0.5/1.222) ;          // Rmin = 1.122 * sigma
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
//-----------------------------------------------------------------------------------------------------

void RandomForce()
{
    double sig ;
    for(int i=0; i<nbacteria ; i++)
    {
        for(int j=0 ; j<nnode; j++)
        {
            Diffusion(bacteria[i].nodes[j].x, bacteria[i].nodes[j].y) ;
            sig = sqrt(2.0*diffusion*dt) ;
            bacteria[i].nodes[j].xdev =gasdev(&idum)*sig;
            bacteria[i].nodes[j].ydev =gasdev(&idum)*sig;
        }
    }
}
//-----------------------------------------------------------------------------------------------------
void Motor()
{
    for (int i=0; i<nbacteria; i++)
    {   for(int j=0 ; j<nnode; j++)
    {
        bacteria[i].nodes[j].fMotorx= 0.0 ;
        bacteria[i].nodes[j].fMotory =0.0 ;
    }
    }
    
    for (int i=0; i<nbacteria; i++)
    {   for(int j=0 ; j<nnode; j++)
    {
        if(j==0)
    {   bacteria[i].nodes[j].fMotorx =  fmotor * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x)/ Distance(i,j,i,j+1) ;
        bacteria[i].nodes[j].fMotory =  fmotor * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y)/ Distance(i,j,i,j+1) ;
        
    }
        else
        {
        bacteria[i].nodes[j].fMotorx = fmotor * (bacteria[i].nodes[j-1].x - bacteria[i].nodes[j].x)/Distance(i,j-1,i,j) ;
        bacteria[i].nodes[j].fMotory = fmotor * (bacteria[i].nodes[j-1].y - bacteria[i].nodes[j].y)/Distance(i,j-1,i,j) ;
        }
    }
    }
}
//-----------------------------------------------------------------------------------------------------
void Reverse ()
{
    for (int i=0; i<nbacteria; i++)
    {
        int start = 0;
        int end = nnode-1 ;
        node temp  ;
        
        while (start < end )
        {
            temp= bacteria[i].nodes[start] ;
            bacteria[i].nodes[start]= bacteria[i].nodes[end] ;
            bacteria[i].nodes[end] = temp ;
            start += 1 ;
            end   -= 1 ;
        }
    }
    ljNodesPosition() ;
}
//-----------------------------------------------------------------------------------------------------
void Initialization ()
{
    int nextLine = static_cast<int>(sqrt(nbacteria/2)) ;
    int row = 0 ;
    int coloum = 0 ;
    double a ;
    double b ;
    double lx = sqrt (2*domainx * domainy / nbacteria ) ;
    double ly = lx ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        if (row%2==0)
        {
            bacteria[i].nodes[j].x = lx * coloum ;
        }
        else   { bacteria[i].nodes[j].x = lx * coloum + lx /2 ; }
        bacteria[i].nodes[j].y=  ly * row /2 ;
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
            bacteria[i].nodes[j+n].x = bacteria[i].nodes[j+n-1].x + x0 * a ; // or nodes[j].x + n * x0 * a
            bacteria[i].nodes[j-n].x = bacteria[i].nodes[j-n+1].x - x0 * a ;
            bacteria[i].nodes[j+n].y = bacteria[i].nodes[j+n-1].y + x0 * b ;
            bacteria[i].nodes[j-n].y = bacteria[i].nodes[j-n+1].y - x0 * b ;
        }
        j=(nnode-1)/2 ;
        if (coloum == nextLine-1)
        {
            row++ ;
            coloum = 0 ;
        }

        else {coloum++ ;}
        
    }
    
}

//-----------------------------------------------------------------------------------------------------
void ljNodesPosition ()
{
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<nnode-1; j++)
        {
            bacteria[i].ljnodes[j].x = (bacteria[i].nodes[j].x + bacteria[i].nodes[j+1].x)/2 ;
            bacteria[i].ljnodes[j].y = (bacteria[i].nodes[j].y + bacteria[i].nodes[j+1].y)/2 ;
        }
    }
}
//-----------------------------------------------------------------------------------------------------
void InitialProtein ()
{
    double allProtein = 0.0 ;
    for (int i=0; i<nbacteria; i++)
    {
       bacteria[i].protein = rand() / (RAND_MAX + 1.0);
        allProtein += bacteria[i].protein ;
        
    }
    cout<<allProtein<<endl ;
    for (int i=0 ; i<nbacteria; i++)
    {
        bacteria[i].protein = bacteria[i].protein/allProtein ;      // normalized protein
    }

}
//-----------------------------------------------------------------------------------------------------
void Connection ()
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

//-----------------------------------------------------------------------------------------------------
void AllNodes ()
{
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j< nnode-1; j++)
        {
            bacteria[i].allnodes[2*j] = bacteria[i].nodes[j] ;
            bacteria[i].allnodes[2*j+1] = bacteria[i].ljnodes[j] ;
            if (j==nnode-2)
            {
                bacteria[i].allnodes[2*(j+1)] = bacteria[i].nodes[j+1] ;
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------------
void ProteinExchange ()
{
    for (int i=0; i<nbacteria-1; i++)
    {
        for (int j=i+1; j<nbacteria; j++)
        {
            if (bacteria[i].connection[j] != 0.0)
            {
                double nbar = (bacteria[i].protein + bacteria[j].protein)/2 ;
                bacteria[i].protein += (nbar-bacteria[i].protein)* Rp * dt * ( bacteria[i].connection[j]/(2*nnode-1) ) ;
                bacteria[j].protein += (nbar-bacteria[j].protein)* Rp * dt * ( bacteria[j].connection[i]/(2*nnode-1) ) ;
            }
        }
    }
    
}


//-----------------------------------------------------------------------------------------------------
void NodeProtein ()
{
    for (int i=0 ; i<nbacteria; i++)
    {
        for (int j=0; j<nnode; j++)
        {
            bacteria[i].nodes[j].protein = bacteria[i].protein ;
        }
    }
}

//-----------------------------------------------------------------------------------------------------
void Duplicate()
{
    for (int i=0; i<nbacteria; i++)
    {   bacteria[i].copy = false ;              // No duplicate is needed
        if (bacteria[i].nodes[(nnode-1)/2].x >5 && bacteria[i].nodes[(nnode-1)/2].x <(domainx-5) && bacteria[i].nodes[(nnode-1)/2].y > 5 && bacteria[i].nodes[(nnode-1)/2].y < (domainy-5) )
        {
            for (int j=0; j<nnode; j++)
            {
                bacteria[i].duplicate[j].x = bacteria[i].nodes[j].x ;
                bacteria[i].duplicate[j].y = bacteria[i].nodes[j].y ;
            }
        }
        else
        {
            for (int j=0; j<nnode; j++)
            {
                if (bacteria[i].nodes[j].x < 0.0)
                {
                    bacteria[i].duplicate[j].x = bacteria[i].nodes[j].x + domainx ;
                    bacteria[i].copy = true ;
                }
                else if (bacteria[i].nodes[j].x > domainx)
                {
                    bacteria[i].duplicate[j].x = bacteria[i].nodes[j].x - domainx ;
                    bacteria[i].copy = true ;
                }
                else
                {
                    bacteria[i].duplicate[j].x = bacteria[i].nodes[j].x ;
                }
                
                if (bacteria[i].nodes[j].y < 0.0)
                {
                    bacteria[i].duplicate[j].y = bacteria[i].nodes[j].y + domainy ;
                    bacteria[i].copy = true ;
                }
                else if (bacteria[i].nodes[j].y > domainy)
                {
                    bacteria[i].duplicate[j].y = bacteria[i].nodes[j].y - domainy ;
                    bacteria[i].copy = true ;
                }
                else
                {
                    bacteria[i].duplicate[j].y = bacteria[i].nodes[j].y ;
                }
               
            }

        }
    }
}
//-----------------------------------------------------------------------------------------------------

void Merge ()
{
    bool ax = true ;
    bool ay = true ;
    for (int i=0 ; i<nbacteria; i++)
    {
        if (bacteria[i].copy== true)
        {
            int j=0 ;
            ax = true ;
            while (j<nnode)
            {
                if ((bacteria[i].nodes[j].x > 0.0 && bacteria[i].nodes[j].x < domainx))
                {   ax = false ;
                    break ;
                }
                j++ ;
            }
            if (ax==true)
            {
                for (j=0 ; j<nnode ; j++)
                {
                    bacteria[i].nodes[j].x = bacteria[i].duplicate[j].x ;
                }
            }
        
            j=0 ;
            ay = true ;
            while (j<nnode)
            {
                if ((bacteria[i].nodes[j].y > 0.0 && bacteria[i].nodes[j].y < domainy))
                {
                    ay = false ;
                    break ;
                }
                j++ ;
            }
            if (ay==true)
            {
                for (j=0 ; j<nnode ; j++)
                {
                    bacteria[i].nodes[j].y = bacteria[i].duplicate[j].y ;
                }
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------
double MinDistance (double x1 , double y1 ,double x2 , double y2)
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

//-----------------------------------------------------------------------------------------------------

void ParaView ()
{
        NodeProtein() ;
        Duplicate() ;
        Merge() ;
        int index = index1 ;
        string vtkFileName = "ECM"+ to_string(index)+ ".vtk" ;
        ofstream ECMOut;
        ECMOut.open(vtkFileName.c_str());
        ECMOut<< "# vtk DataFile Version 3.0" << endl;
        ECMOut<< "Result for paraview 2d code" << endl;
        ECMOut << "ASCII" << endl;
        ECMOut << "DATASET UNSTRUCTURED_GRID" << endl;
        ECMOut << "POINTS " << points << " float" << endl;
        for (uint i = 0; i < nbacteria; i++)
        {
            for (uint j=0; j<nnode ; j++)
            {
                ECMOut << bacteria[i].duplicate[j].x << " " << bacteria[i].duplicate[j].y << " "
                << 0.0 << endl;
                
            }
        }
        ECMOut<< endl;
        ECMOut<< "CELLS " << points-1<< " " << 3 *(points-1)<< endl;
        for (uint i = 0; i < (points-1); i++)
        {
            
            ECMOut << 2 << " " << i << " "
            << i+1 << endl;
            
        }
    
        ECMOut << "CELL_TYPES " << points-1<< endl;
        for (uint i = 0; i < points-1; i++) {
            ECMOut << "3" << endl;
        }
        ECMOut << "POINT_DATA "<<points <<endl ;
        ECMOut << "SCALARS Protein_rate " << "float"<< endl;
        ECMOut << "LOOKUP_TABLE " << "default"<< endl;
        for (uint i = 0; i < nbacteria ; i++)
        {   for ( uint j=0; j<nnode ; j++)
        {
            ECMOut<< bacteria[i].nodes[j].protein <<endl ;
        }
        }
        ECMOut.close();
}

    
    //-----------------------------------------------------------------------------------------------------

   void ParaView2 ()
{
    int index = index1 ;
    string vtkFileName2 = "Grid"+ to_string(index)+ ".vtk" ;
    ofstream SignalOut;
    SignalOut.open(vtkFileName2.c_str());
    SignalOut << "# vtk DataFile Version 2.0" << endl;
    SignalOut << "Result for paraview 2d code" << endl;
    SignalOut << "ASCII" << endl;
    SignalOut << "DATASET RECTILINEAR_GRID" << endl;
    SignalOut << "DIMENSIONS" << " " << nx  << " " << " " << ny << " " << nz  << endl;
    
    SignalOut << "X_COORDINATES " << nx << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int i = 0; i < nx ; i++) {
        SignalOut << X[i] << endl;
    }
    
    SignalOut << "Y_COORDINATES " << ny << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int j = 0; j < ny; j++) {
        SignalOut << Y[j] << endl;
    }
    
    SignalOut << "Z_COORDINATES " << nz << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int k = 0; k < nz ; k++) {
        SignalOut << 0 << endl;
    }
    
    SignalOut << "POINT_DATA " << (nx )*(ny )*(nz ) << endl;
    SignalOut << "SCALARS DPP float 1" << endl;
    SignalOut << "LOOKUP_TABLE default" << endl;
    
    for (int k = 0; k < nz ; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                SignalOut << slime[i][j] << endl;
                
            }
        }
    }
    index1++ ;
}
   
//-----------------------------------------------------------------------------------------------------
void Track ()
{
    int m = 0 ;
    int n = 0 ;

    for (m=0; m<nx; m++)
    {
        for (n=0; n<ny; n++)
        {
            slime[m][n] -= kd * dt * (slime[m][n]-1)  ;
        }
    }
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<2*nnode-1 ; j++)
        {
           // we add domain to x and y in order to make sure m and n would not be negative integers
            m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].x + domainx , domainx) / dx ) ) ) % nx  ;
            n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].y + domainy , domainy) / dy ) ) ) % ny  ;
            slime [m][n] += sr* dt ;
        }
    }
}
//-----------------------------------------------------------------------------------------------------
void Diffusion (double x , double y)
{
    int m = (static_cast<int> (round ( fmod (x + domainx , domainx) / dx ) ) ) % nx ;
    int n = (static_cast<int> (round ( fmod (y + domainy , domainy) / dy ) ) ) % ny ;
  /*
    if (slime[m][n] == 0.0  )
    {
        diffusion = kblz * temp / eta2 ;
    }
    else
    {
        diffusion = kblz * temp / eta1 ;
    }
   */
    diffusion = kblz * temp/ (eta1/ slime[m][n] ) ;
}


