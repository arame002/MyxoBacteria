

#include "TissueBacteria.hpp"
#include "ranum2.h"
#include "SignalCommon.hpp"

long  idum=(-799);



// Inverse function is on inside the InverseTime
// Function: WriteInFile, RunMode ,

//void Myxo () ;
//void Spring() ;
//void Bending() ;
//double Distance (int ,int, int, int) ;                                     // bacteria, node, bacteria, node
//double Distance2 (double x1, double y1 , double x2, double y2, int shiftx, int shifty ) ;
//double MinDistance (double x1 , double y1 , double x2 , double y2) ;
//double Cos0ijk(int, int) ;
//double u_lj () ;
void RandomForce() ;
void Motor () ;
void Reverse (int i) ;                                               // inverse head and tail
void Initialization () ;
void Initialization2 () ;
void CircularInitialization () ;
void PositionUpdating (double  ) ;                               // (time step)
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
void SlimeTrace () ;
void Diffusion (double, double ) ;
void InitialReversalTime () ;
void ReversalTime() ;
double SlimeTrailFollowing (int , double) ;
//void PiliForce () ;
void InitialPili () ;
double AngleOfVector (double x1, double y1 , double x2 , double y2 ) ;      // (x1,y1,x2,y2)
void SurfaceCoverage () ;
void VisitsPerGrid () ;
int PowerLawExponent () ;

void TurnOrientation () ;
double Cal_OrientationBacteria (int i) ;
void UpdateReversalFrequency () ;
double Cal_ChemoGradient (int i ) ;
double Cal_ChemoGradient2 (int i) ;
void WriteTrajectoryFile () ;
void WriteNumberReverse () ;
void SlimeTraceHyphae (Fungi tmpFng) ;


//-----------------------------------------------------------------------------------------------------
//Bacteria properties
double length = 2 ;
double B = 0.5 ;                      // bending constant
double tetta0=3.1415 ;               // prefered angle, Pi
double K =  1.0 ;                      // linear spring constant (stiffness)
double x0 = length/(nnode-1) ;                      // equilibrium distance of spring
double fmotor = 0.1 ;                           //   Ft/(N-1)
double Rp = 0.3 ;                               // protein exchange rate
double turnPeriod = 50.0 ;

double reversalRate = 0.025 ;
double chemoStrength = 10.0 ;

//-----------------------------------------------------------------------------------------------------
//domain and medium properties
double kblz=1.0 , temp=0.01 , eta1 = 1.0 , eta2 = 5.0 ;
double diffusion ;
int shiftx = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
int shifty = 0 ;        //this value changes in the MinDistance function, call MinDistance before using this value
const double dx= 0.5 ;
const double dy= 0.5 ;
const int nx = static_cast<int>(domainx/dx) ;               //number of grids in X-axis
const int ny = static_cast<int> (domainy/dy) ;              //number of grids in Y-axis
const int nz = 1 ;
double X[nx] ;
double Y[ny] ;
double Z[nz]= {0.0} ;
int visit[nx][ny] ;
vector<double> frequency ;
int numOfClasses = 1 ;
int fNoVisit = 0 ;
vector<vector<double> > gridInMain ;

//-----------------------------------------------------------------------------------------------------
//slime properties
double sr = 0.25 ;                          // slime rate production
double kd = 0.05 ;                           // slime decay rate
double slimeEffectiveness = 5.0 ;           // needed for slime trail following
double s0 = 1 ;
double slime[nx][ny] ;
int surfaceCoverage[nx][ny] ;
double coveragePercentage = 0 ;
int searchAreaForSlime = static_cast<int> (round((length/2.0)/ min(dx , dy))) ;


//-----------------------------------------------------------------------------------------------------
//pili properties
double piliMaxLength = 2.5 ;                       //pili maximum length
double vProPili = 0.5 ;                                 // constant protrude velocity
double vRetPili0 = 0.5 ;                                // constant retraction rate when the pili is not connected
double kPullPili = 2000 ;
double fStall = 180.0 ;
double averageLengthFree = 0.0 ;
double nAttachedPili = 0.0 ;

//-----------------------------------------------------------------------------------------------------
//simulation parameters
double dt=0.0001 ;
double initialStep = 0.0001 ;
int index1 = 0 ;                                     // needed for ParaView

//-----------------------------------------------------------------------------------------------------
double initialTime = 4.0 ;
double runTime = 2000.0 ;
    bacterium bacteria[nbacteria];

int main ()
{
    TissueBacteria tissueBacteria ;
    srand(time (0)) ;
    cout<<"Search area for slime is "<<searchAreaForSlime<<endl ;
    for (int m=0; m<nx; m++)
    {
        for (int n=0; n<ny; n++)
        {
            slime[m][n] = s0  ;
            surfaceCoverage[m][n] = 0 ;
            visit[m][n] = 0 ;
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
    ofstream FrequencyOfVisit ;
    FrequencyOfVisit.open("FrequncyOfVisit.txt") ;
    
    cout<<"program is running"<<endl  ;
    double nt=runTime/dt + initialTime/initialStep ;                   //total number of steps
    nt =static_cast<int>(nt) ;
    double initialNt =initialTime/initialStep ;                     // run time for initialization
    initialNt =static_cast<int>(initialNt) ;
    //    int inverseInitialStep = static_cast<int>(initialNt/initialTime) ;  // used for visualization(ParaView)
    int inverseDt =static_cast<int>((nt-initialNt)/(runTime)) ;              // used for visualization(ParaView)
    //Myxo() ;
    //Initialization() ;
    Initialization2() ;
    //CircularInitialization() ;
    ljNodesPosition() ;
    InitialProtein() ;
    InitialReversalTime() ;
    InitialPili () ;
    for (int i=0; i<nbacteria; i++)
    {
        bacteria[i].connection[i]=0.0 ;
        for (int j=0 ; j<nnode ; j++)
        {
            bacteria[i].nodes[j].xdev = 0.0 ;
            bacteria[i].nodes[j].ydev = 0.0 ;
            
        }
    }
    Fungi fungi = driver() ;
    vector<HyphaeSegment> hyphaeSegments_main = fungi.hyphaeSegments ;
    
    vector<vector<double> > pointSources = fungi.tips ;
    gridInMain = tissueBacteria.Cal_Diffusion2D(0, 100, 0, 100,fungi.tips) ;
    SlimeTraceHyphae(fungi) ;
    for (int l=0; l< (nt+1); l++)
    {
        if (l < initialNt)
        {
            
            AllNodes() ;
            tissueBacteria.Spring () ;
            tissueBacteria.Bending() ;
            //      u_lj() ;
            //      RandomForce() ;
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
            ReversalTime() ;
            AllNodes() ;
            tissueBacteria.Spring () ;
            tissueBacteria.Bending() ;
          //  u_lj() ;
            //  RandomForce() ;
            Motor () ;
            TurnOrientation() ;
            //  SlimeTrace() ;
            //  Connection() ;
            //  ProteinExchange() ;
            //  PiliForce() ;
            
            
            if (l%inverseDt==0)
            {
                /*
                 string NumberOfVisits = "NumberOfVisits"+ to_string(index1)+ ".txt" ;
                 ofstream NumberOfVisit;
                 NumberOfVisit.open(NumberOfVisits.c_str());
                 SurfaceCoverage() ;
                 VisitsPerGrid() ;
                 int alfaMin = PowerLawExponent() ;
                 
                 for (int m=0; m< nx; m++)
                 {
                 for (int n=0 ; n < ny ; n++)
                 {
                 NumberOfVisit<< visit[m][n]<<'\t' ;
                 }
                 NumberOfVisit<<endl ;
                 }
                 NumberOfVisit<<endl<<endl ;
                 NumberOfVisit.close() ;
                 
                 for (int i=0; i<nbacteria; i++)
                 {
                 ProteinLevelFile<<bacteria[i].protein<<"," ;
                 }
                 ProteinLevelFile<< endl<<endl ;
                 FrequencyOfVisit<< "Results of frame "<< index1<<" are below : "<<endl ;
                 FrequencyOfVisit<< "Surface coverage is equal to:     "<<coveragePercentage <<endl ;
                 FrequencyOfVisit<< "alfaMin is equal to:    " << alfaMin<<endl ;
                 FrequencyOfVisit<<"size of FrequencyOfVisit is equal to:  "<< numOfClasses<<endl ;
                 FrequencyOfVisit<<"Frequecny of no visit is equal to:    "<<fNoVisit<<endl ;
                 
                 for (int i=0 ; i<numOfClasses ; i++)
                 {
                 FrequencyOfVisit<< frequency[i]<<'\t' ;
                 }
                 FrequencyOfVisit<<endl<<endl ;
                 */
                ParaView ()  ;
                ParaView2() ;
                WriteTrajectoryFile() ;
                cout<<(l-initialNt)/inverseDt<<endl ;
                //   cout << averageLengthFree<<'\t'<<nAttachedPili<<endl ;
                //    cout << coveragePercentage <<endl ;
            }
            PositionUpdating(dt) ;
            //UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
        }
        
    }
    
    WriteNumberReverse() ;
    
    cout<<"No Bug"<<endl ;
    return 0 ;
}


/*_____________________________________________________________________________
 ______________________________________________________________________________*/

void PositionUpdating (double t)
{   double etaInverse = diffusion / (kblz * temp) ;
    double totalForceX = 0 ;
    double totalForceY = 0 ;
    for (int i=0; i<nbacteria; i++)
    {
        bacteria[i].oldLoc.at(0) = bacteria[i].nodes[(nnode-1)/2].x ;
        bacteria[i].oldLoc.at(1) = bacteria[i].nodes[(nnode-1)/2].y ;
        int tmpXIndex = static_cast<int>( round( bacteria[i].nodes[(nnode-1)/2].x / 0.5 ) ) ; //0.5 is grid size, need to change it as a variable
        int tmpYIndex = static_cast<int>( round( bacteria[i].nodes[(nnode-1)/2].y / 0.5 ) ) ;
        //test No chemo
        //bacteria[i].oldChem = gridInMain.at(tmpYIndex).at(tmpXIndex) ;
        for(int j=0 ; j<nnode; j++)
        {   Diffusion(bacteria[i].nodes[j].x , bacteria[i].nodes[j].y) ;
            etaInverse = diffusion / (kblz * temp) ;
            totalForceX = bacteria[i].nodes[j].fSpringx ;       // it is not += because the old value is related to another node
            totalForceX += bacteria[i].nodes[j].fBendingx ;
            //     totalForceX += bacteria[i].nodes[j].fljx ;
            totalForceX += bacteria[i].nodes[j].fMotorx ;
            bacteria[i].nodes[j].x += t * totalForceX * etaInverse ;
            bacteria[i].nodes[j].x += bacteria[i].nodes[j].xdev ;
            
            totalForceY = bacteria[i].nodes[j].fSpringy ;       // it is not += because the old value is related to another node
            totalForceY += bacteria[i].nodes[j].fBendingy ;
            //     totalForceY += bacteria[i].nodes[j].fljy ;
            totalForceY += bacteria[i].nodes[j].fMotory ;
            bacteria[i].nodes[j].y += t * totalForceY * etaInverse ;
            bacteria[i].nodes[j].y += bacteria[i].nodes[j].ydev ;
            
            
            
            // F pili
            
            if (j==0)
            {
                for (int m=0 ; m<nPili ; m++)
                {
                    bacteria[i].nodes[j].x += t * (bacteria[i].pili[m].fx) * etaInverse ;
                    bacteria[i].nodes[j].y += t * (bacteria[i].pili[m].fy) * etaInverse ;
                }
            }
            
        }
    }
    ljNodesPosition() ;
}

//-----------------------------------------------------------------------------------------------------
/*
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
 */
//-----------------------------------------------------------------------------------------------------
double Distance (int i,int j, int k, int l)
{
    double dis= sqrt ((bacteria[i].nodes[j].x - bacteria[k].nodes[l].x) * (bacteria[i].nodes[j].x-bacteria[k].nodes[l].x) + (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y) * (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y)) ;
    return dis ;
}
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double Cos0ijk(int i,int m)
{
    double dot= (bacteria[i].nodes[m-1].x - bacteria[i].nodes[m].x)*(bacteria[i].nodes[m+1].x - bacteria[i].nodes[m].x);
    dot += (bacteria[i].nodes[m-1].y - bacteria[i].nodes[m].y)*(bacteria[i].nodes[m+1].y -bacteria[i].nodes[m].y) ;
    double cos= dot / ( Distance(i,m-1,i,m) * Distance(i,m+1,i,m) ) ;
    return cos ;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
/*
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
 */
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
    {
        double tmpFmotor ;
        if (bacteria[i].directionOfMotion == true)
        {
            tmpFmotor = fmotor ;
        }
        else
        {
            tmpFmotor = 0.5 * fmotor ;
        }
        if (bacteria[i].attachedToFungi == false)
        {
            tmpFmotor *= 0.25 ;
        }
        for(int j=0 ; j<nnode; j++)
        {
            if(j==0)
            {
                double fs = 0.0 ;
                //fungi is modeled by slime
                //Following Slime or hyphae
                 fs = SlimeTrailFollowing(i, length/2) ;
                
                
                bacteria[i].nodes[j].fMotorx =  (tmpFmotor - fs ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x)/ Distance(i,j,i,j+1) ;                          // force to the direction of head
                
                bacteria[i].nodes[j].fMotory =  (tmpFmotor - fs ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y)/ Distance(i,j,i,j+1) ;                          // force to the direction of head
                
                bacteria[i].nodes[j].fMotorx += bacteria[i].fSTFx ;     // Fh + Fs
                bacteria[i].nodes[j].fMotory += bacteria[i].fSTFy ;
            }
            else
            {
                bacteria[i].nodes[j].fMotorx = tmpFmotor * (bacteria[i].nodes[j-1].x - bacteria[i].nodes[j].x)/Distance(i,j-1,i,j) ;
                bacteria[i].nodes[j].fMotory = tmpFmotor * (bacteria[i].nodes[j-1].y - bacteria[i].nodes[j].y)/Distance(i,j-1,i,j) ;
            }
        }
        int m = (static_cast<int> (round ( fmod (bacteria[i].nodes[(nnode-1)/2].x + domainx , domainx) / dx ) ) ) % nx  ;
        int n = (static_cast<int> (round ( fmod (bacteria[i].nodes[(nnode-1)/2].y + domainy , domainy) / dy ) ) ) % ny  ;
        if (slime[m][n] > s0 + 0.9* sr * dt )
        {
            bacteria[i].attachedToFungi = true ;
            bacteria[i].reversalPeriod = 300 ;
            
        }
        else if (bacteria[i].attachedToFungi == true)
        {
            //bacteria[i].attachedToFungi = false ;
            bacteria[i].reversalPeriod = 60 ;
            //Reverse(i) ;
            
        }
        
    }
}
//-----------------------------------------------------------------------------------------------------
void Reverse (int i)
{
   // bacteria[i].directionOfMotion = ! bacteria[i].directionOfMotion ;
    bacteria[i].turnStatus = true ;
    bacteria[i].turnTime = 0.0 ;
  //  bacteria[i].turnAngle =(2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) * (3.1415) ;    //uniform distribution
    if (bacteria[i].attachedToFungi == false)
    {
        bacteria[i].turnAngle = gasdev(&idum) * (3.1415)/15 ;                //guassian distribution
    }
    else
    {
        bacteria[i].turnAngle = 0.0 ;
    }
    //   bacteria[i].turnAngle = 3.1415/3.0 ;
    bacteria[i].SourceRegion() ;
    if (bacteria[i].sourceWithin == false)
    {
        bacteria[i].numberReverse += 1 ;
    }
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
    for (int j=0; j<nnode-1; j++)                                   //updating L-J nodes positions
    {
        bacteria[i].ljnodes[j].x = (bacteria[i].nodes[j].x + bacteria[i].nodes[j+1].x)/2 ;
        bacteria[i].ljnodes[j].y = (bacteria[i].nodes[j].y + bacteria[i].nodes[j+1].y)/2 ;
    }
    
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
void Initialization2 ()
{
    int nextLine = static_cast<int>(ceil( sqrt(nbacteria/2)/2) )  ;
    int row = -nextLine ;
    int coloum = -nextLine ;
    double tmpDomainX = domainx* 0.5 ;
    double tmpDomainY = domainy * 0.5 ;
    double a ;
    double b ;
    double lx = sqrt (2*tmpDomainX * tmpDomainY / nbacteria ) ;
    double ly = lx ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        if (row%2==0)
        {
            bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum ;
        }
        else   { bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum + lx /2.0 ; }
        bacteria[i].nodes[j].y=  domainy/2.0 + ly * row /2.0 ;
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
            coloum = -nextLine ;
        }
        
        else {coloum++ ;}
        
    }
    
}

//-----------------------------------------------------------------------------------------------------
void CircularInitialization ()
{
    double raduis = 0.25 * sqrt(domainx * domainx  )/2.0 ;
    double deltaTetta = 2.0 * 3.1416 / nbacteria ;
    double cntrX = domainx/2.0 ;
    double cntrY = domainy/2.0 ;
    
    double a ;
    double b ;
    double lx = sqrt (2*domainx * domainy / nbacteria ) ;
    double ly = lx ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        bacteria[i].nodes[j].x = cntrX + raduis * cos( static_cast<double>(i* deltaTetta) ) ;
        bacteria[i].nodes[j].y = cntrY + raduis * sin( static_cast<double>(i* deltaTetta) ) ;
        
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
    cout<<"Total amount of Protein is "<< allProtein<<endl ;
    for (int i=0 ; i<nbacteria; i++)
    {
        bacteria[i].protein = bacteria[i].protein/allProtein ;      // normalized protein
    }
    
}
//-----------------------------------------------------------------------------------------------------
/*
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
*/
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
        if (bacteria[i].nodes[(nnode-1)/2].x > length && bacteria[i].nodes[(nnode-1)/2].x <(domainx-length) && bacteria[i].nodes[(nnode-1)/2].y > length && bacteria[i].nodes[(nnode-1)/2].y < (domainy-length) )
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
/*
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
*/
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
    
    for (uint i = 0; i < (points-1); i++)           //number of connections per node
    {
        
        ECMOut << 2 << " " << i << " "
        << i+1 << endl;
        
    }
    
    ECMOut << "CELL_TYPES " << points-1<< endl;             //connection type
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
    //index1++ ;
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
void SlimeTrace ()
{
    int m = 0 ;
    int n = 0 ;
    
    for (m=0; m<nx; m++)
    {
        for (n=0; n<ny; n++)
        {
            slime[m][n] -= kd * dt * (slime[m][n]-s0)  ;
        }
    }
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<2*nnode-1 ; j++)
        {
            // we add domain to x and y in order to make sure m and n would not be negative integers
            m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].x + domainx , domainx) / dx ) ) ) % nx  ;
            n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].y + domainy , domainy) / dy ) ) ) % ny  ;
            slime [m][n] += sr * dt ;
            //   visit[m][n] += 1.0 ;      //undating visits in each time step. if you make it as uncomment, you should make it comment in the VisitsPerGrid function
            
        }
    }
}
//-----------------------------------------------------------------------------------------------------
void Diffusion (double x , double y)
{
    //int m = (static_cast<int> (round ( fmod (x + domainx , domainx) / dx ) ) ) % nx ;
    //  int n = (static_cast<int> (round ( fmod (y + domainy , domainy) / dy ) ) ) % ny ;
    //  diffusion = kblz * temp/ (eta1/ slime[m][n] ) ;
    diffusion = kblz * temp/ eta1 ;
}
//-----------------------------------------------------------------------------------------------------
void InitialReversalTime ()
{
    for( int i=0 ; i<nbacteria ; i++)
    {
        bacteria[i].reversalTime = (rand() / (RAND_MAX + 1.0)) * bacteria[i].reversalPeriod ;
        bacteria[i].reversalTime -= fmod(bacteria[i].reversalTime , dt ) ;
        if (rand() / (RAND_MAX + 1.0) < 0.5)
        {
            //bacteria[i].directionOfMotion = false ;
            bacteria[i].directionOfMotion = true ;
        }
        else
        {
            bacteria[i].directionOfMotion = true ;
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------
void ReversalTime()
{
    for (int i=0; i<nbacteria; i++)
    {
        
        bacteria[i].sourceWithin = bacteria[i].SourceRegion() ;
        if (bacteria[i].sourceWithin == false)
        {
            bacteria[i].timeToSource += dt ;
        }
        
        if (bacteria[i].reversalTime >= bacteria[i].reversalPeriod )
        {
            //bacteria[i].reversalTime -= bacteria[i].reversalPeriod ;
            bacteria[i].reversalTime  = 0.0 ;
            
            Reverse(i) ;
        }
        bacteria[i].reversalTime += dt ;
    }
}
//-----------------------------------------------------------------------------------------------------
double SlimeTrailFollowing (int i, double R)        // i th bacteria, radius of the search area
{
    int m=0   ;
    int n=0   ;
    int gridX ;
    int gridY ;
    double deltax ;
    double deltay ;
    double r ;
    int region = 2 ;
    double totalSlimeInSearchArea = 0.0 ;
    double gridInRegions [5] = { } ;                // number of grids in each region
    double gridWithSlime [5] = { } ;                   // number of grids containing slime
    double slimeRegions[5] = { }   ;                // the amount of slime in each region
    /*
     0:  area between 0 and 36
     1:  area between 36 and 72
     2:  area between 72 and 108
     3:  area between 108 and 144
     4:  area between 144 and 180
     */
    
    double ax = bacteria[i].nodes[0].x - bacteria[i].nodes[1].x ;
    double ay = bacteria[i].nodes[0].y - bacteria[i].nodes[1].y ;
    double Cos = ax/sqrt(ax*ax+ay*ay) ;
    double orientationBacteria = acos(Cos)* 180 / 3.1415 ;        // orientation of the bacteria
    if(ay<0) orientationBacteria *= -1 ;
    double alfa ;
    
    m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[0].x + domainx , domainx) / dx ) ) ) % nx  ;
    n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[0].y + domainy , domainy) / dy ) ) ) % ny  ;
    for (int sx = -1*searchAreaForSlime ; sx <= searchAreaForSlime ; sx++)
    {
        for (int sy = -1*searchAreaForSlime ; sy<= searchAreaForSlime ; sy++ )
        {
            gridX = (m+sx) % nx ;
            gridY = (n+sy) % ny ;
            ax = sx * dx + dx/2 ;
            ay = sy * dy + dy/2 ;
            Cos = ax/ sqrt(ax * ax + ay * ay) ;
            alfa = acos(Cos)*180 / 3.1415 ;
            if (ay<0) alfa *= -1 ;
            double gridAngleInSearchArea =  fmod ( alfa - orientationBacteria, 360 ) ;       // using alfa for searching among grids, Using fmod, because the result should be in the range (-180, 180 ) degree
            
            if (gridAngleInSearchArea >= -90 && gridAngleInSearchArea <= 90 )
            {
                deltax = ((m+sx)*dx + dx/2 ) - bacteria[i].nodes[0].x ;
                deltay = ((n+sy)*dy + dy/2 ) - bacteria[i].nodes[0].y ;
                r = sqrt (deltax*deltax + deltay*deltay) ;
                if ( r < R)
                {
                    if (gridAngleInSearchArea > -90 && gridAngleInSearchArea <= -54)
                    {
                        
                        slimeRegions[0] += slime[gridX][gridY] ;
                        gridInRegions[0] += 1 ;
                        if ( slime[gridX][gridY] - s0 >= 0.0000001)
                        {
                            gridWithSlime[0] += 1 ;
                        }
                    }
                    else if (gridAngleInSearchArea > -54 && gridAngleInSearchArea <= -18)
                    {
                        
                        slimeRegions[1] += slime[gridX][gridY] ;
                        gridInRegions[1] += 1 ;
                        if ( slime[gridX][gridY] - s0  >= 0.0000001)
                        {
                            gridWithSlime[1] += 1 ;
                        }
                    }
                    else if (gridAngleInSearchArea> -18 && gridAngleInSearchArea <= 18)
                    {
                        
                        slimeRegions[2] += slime[gridX][gridY] ;
                        gridInRegions[2] += 1 ;
                        if ( slime[gridX][gridY] - s0  >= 0.0000001)
                        {
                            gridWithSlime[2] += 1 ;
                        }
                    }
                    
                    else if (gridAngleInSearchArea > 18 && gridAngleInSearchArea <= 54)
                    {
                        
                        slimeRegions[3] += slime[gridX][gridY] ;
                        gridInRegions[3] += 1 ;
                        if ( slime[gridX][gridY] - s0  >= 0.0000001)
                        {
                            gridWithSlime[3] += 1 ;
                        }
                    }
                    
                    else if (gridAngleInSearchArea > 54 && gridAngleInSearchArea <= 90)
                    {
                        
                        slimeRegions[4] += slime[gridX][gridY] ;
                        gridInRegions[4] += 1 ;
                        if ( slime[gridX][gridY] - s0 >= 0.0000001)
                        {
                            gridWithSlime[4] += 1 ;
                        }
                    }
                    
                    totalSlimeInSearchArea += slime[gridX][gridY] ;
                    
                }
                
            }
            
        }
    }
    
    if (totalSlimeInSearchArea < 0.000001)
    {
        bacteria[i].fSTFx = 0.0 ;
        bacteria[i].fSTFy = 0.0 ;
        // There is no slime in the search area
        
    }
    /*
     else
     {   double mostFilledRegion= 0.0 ;
     for (int i=0 ; i < 5 ; i++)
     {
     if (( gridWithSlime[i] / gridInRegions[i] ) > 0.8 )
     {
     if ( abs(i-2) < abs(region-2) )             // lower change in angle
     {
     mostFilledRegion = slimeRegions[i] ;
     region = i ;
     }
     else if ( ( abs(i-2) == abs(region-2) ) && rand() % 2 ==0 )     //choosing between two randomly
     {
     mostFilledRegion = slimeRegions[i] ;
     region = i ;
     }
     }
     }
     if ( region ==6) region = 2 ;                       // re-initialization after checking the conditions
     alfa = ( (region -2 ) * 36 + orientationBacteria ) *3.1415/180 ;
     // using alfa(radiant) for direction of the force
     // Fs = EpsilonS *(ft/(N-1))* (S/S0) ( toward slime)
     
     bacteria[i].fSTFx = slimeEffectiveness * fmotor * cos(alfa) ;
     bacteria[i].fSTFy = slimeEffectiveness * fmotor * sin(alfa) ;
     // (S/S0) is not included yet
     }
     */
    else
    {   double maxSlimeRegion= 0.0 ;
        for (int i=0 ; i < 5 ; i++)
        {
            if (slimeRegions[i] > maxSlimeRegion )
            {
                maxSlimeRegion = slimeRegions[i] ;
                region = i ;
            }
        }
        for (int i=0 ; i<5 ; i++)
        {
            if ( (slimeRegions[i] > 0.8 * maxSlimeRegion))
            {
                if ( abs(i-2) < abs(region-2) )
                {
                    region = i ;
                }
                
                //       else if ( abs(i-2) == abs(region-2) && rand() % 2 ==0  )
                else if ( abs(i-2) == abs(region-2) && slimeRegions[i] > slimeRegions[region])
                {
                    region = i ;
                }
                
            }
        }
        alfa = ( (region -2 ) * 36 + orientationBacteria ) *3.1415/180 ;
        // using alfa(radiant) for direction of the force
        // Fs = EpsilonS *(ft/(N-1))* (S/S0) ( toward slime)
        
        bacteria[i].fSTFx = slimeEffectiveness * fmotor * cos(alfa) ;
        bacteria[i].fSTFy = slimeEffectiveness * fmotor * sin(alfa) ;
        // (S/S0) is not included yet
    }
    return  sqrt(bacteria[i].fSTFx * bacteria[i].fSTFx + bacteria[i].fSTFy * bacteria[i].fSTFy) ;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void InitialPili ()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for ( int j=0 ; j<nPili ; j++)
        {
            bacteria[i].pili[j].lFree = rand() / (RAND_MAX + 1.0) * piliMaxLength ;
            bacteria[i].pili[j].attachment = false ;
            bacteria[i].pili[j].retraction = false ;
            bacteria[i].pili[j].fx = 0.0 ;
            bacteria[i].pili[j].fy = 0.0 ;
            bacteria[i].pili[j].F = 0.0 ;
            
            if(rand() / (RAND_MAX + 1.0) < 0.5 )
            {
                bacteria[i].pili[j].retraction = true ;
            }
            else
            {
                bacteria[i].pili[j].retraction = false ;
            }
            if (rand() / (RAND_MAX + 1.0) < 0.5 )
            {
                bacteria[i].pili[j].attachment = true ;
                bacteria[i].pili[j].retraction = true ;
            }
            else
            {
                bacteria[i].pili[j].attachment = false ;
            }
            
        }
    }
}

//-----------------------------------------------------------------------------------------------------
/*
void PiliForce ()
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
//------------------------------------------------------------------------------------------------------------
double AngleOfVector (double x1,double y1,double x2,double y2)
{
    double ax = x2 - x1 ;
    double ay = y2 - y1 ;
    double Cos = ax/sqrt(ax*ax+ay*ay) ;
    double angle = acos(Cos)* 180.0 / 3.1415 ;
    if(ay<0) angle *= -1 ;
    return angle ;                                   // (-180,180)
    
}
//------------------------------------------------------------------------------------------------------------
void SurfaceCoverage ()
{
    double percentage = 0 ;
    
    for (int m=0; m< nx; m++)
    {
        for (int n=0 ; n < ny ; n++)
        {
            if (slime[m][n] > s0 && surfaceCoverage[m][n] == 0 )
            {
                surfaceCoverage[m][n] = 1 ;
                percentage += 1 ;
            }
        }
    }
    percentage = percentage / (nx * ny * 1.0) ;
    coveragePercentage += percentage ;
    
}

//------------------------------------------------------------------------------------------------------------
void VisitsPerGrid ()
{
    //SlimeTrace is needed for this function, it calculate #visit in each grid
    int m = 0 ;
    int n = 0 ;
    //  int newVisit[nx][ny] = {0} ;
    
    
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<2*nnode-1 ; j++)
        {
            // we add domain to x and y in order to make sure m and n would not be negative integers
            m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].x + domainx , domainx) / dx ) ) ) % nx  ;
            n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].y + domainy , domainy) / dy ) ) ) % ny  ;
            visit[m][n] += 1 ;
            //  newVisit[m][n] = 1 ;
            
        }
    }
    /*
     for (int m = 0 ; m< nx ; m++)
     {
     for (int n = 0; n<ny; n++)
     {
     visit[m][n] = visit[m][n] + newVisit[m][n] ;
     }
     
     }
     */
    
}

//------------------------------------------------------------------------------------------------------------
int PowerLawExponent ()
{
    double discrete = 1.5 ;
    double log10discrete = log10(discrete) ;
    fNoVisit = 0 ;
    int minValue = visit[0][0] ;
    int maxValue = visit[0][0] ;
    int alfaMin = 0 ;
    int alfaMax = 0 ;
    
    // alfa and a must be integer, if you want to change them to a double variable you need to make some changes in the code
    int alfa = 0 ;
    for (int n = 0; n<ny; n++)
    {
        for (int m = 0 ; m<nx; m++)
        {
            if(visit[m][n] > maxValue ) maxValue = visit[m][n];
            if(visit[m][n] < minValue ) minValue = visit[m][n];
        }
    }
    alfaMin = max(floor (log10(minValue)/log10discrete) , 0.0 ) ;
    alfaMax = max(floor (log10(maxValue)/log10discrete) , 0.0 );
    numOfClasses = (alfaMax - alfaMin +1 )  ;
    frequency.assign(numOfClasses, 0) ;
    
    for (int n = 0; n<ny; n++)
    {
        for (int m = 0 ; m<nx; m++)
        {   if(visit[m][n] == 0)
        {
            fNoVisit += 1 ;
            continue ;
        }
            alfa = floor(log10(visit[m][n])/log10discrete) ;
            frequency.at(alfa ) += 1 ;
        }
    }
    double normalization = 1 ;
    for (int i = 0 ; i<numOfClasses ; i++)
    {
        normalization = floor(pow(discrete, i+1)) - ceil(pow(discrete, i)) + 1 ;    // it counts integers between this two numbers
        
        if(normalization == 0 )    normalization = 1 ;  // for that i, frequesncy is zero. normalized or unnormalized .
        frequency.at(i) = frequency.at(i) / normalization ;
    }
    
    return alfaMin ;
}

//------------------------------------------------------------------------------------------------------------
void TurnOrientation ()
{
    for (int i=0; i<nbacteria; i++)
    {
        if (bacteria[i].turnStatus)
        {
            // bacteria[i].turnAngle = 0.78 ;
            double fTotal =sqrt(bacteria[i].nodes[0].fMotorx * bacteria[i].nodes[0].fMotorx +
                                bacteria[i].nodes[0].fMotory * bacteria[i].nodes[0].fMotory) ;
            double orientation = Cal_OrientationBacteria(i) ;
            bacteria[i].nodes[0].fMotorx = fTotal * cos(bacteria[i].turnAngle + orientation ) ;
            bacteria[i].nodes[0].fMotory = fTotal * sin(bacteria[i].turnAngle + orientation ) ;
            bacteria[i].turnTime += dt ;
            
            if (bacteria[i].turnTime > turnPeriod)
            {
                bacteria[i].turnStatus = false ;
                bacteria[i].turnTime = 0.0 ;
                bacteria[i].turnAngle = 0.0 ;
            }
            
        }
    }
}

//------------------------------------------------------------------------------------------------------------
double Cal_OrientationBacteria (int i)

{
    
    double ax = 1.0 *(bacteria[i].nodes[1].x - bacteria[i].nodes[0].x) ;
    double ay = 1.0 *(bacteria[i].nodes[1].y - bacteria[i].nodes[0].y) ;
    double Cos = ax/sqrt(ax*ax+ay*ay) ;
    double orientationBacteria = acos(Cos)  ;        // orientation of the bacteria in radian
    if(ay<0.0) orientationBacteria *= -1.0 ;
    bacteria[i].orientation = orientationBacteria ;
    return orientationBacteria ;
}

//------------------------------------------------------------------------------------------------------------
void UpdateReversalFrequency ()
{
    for (int i=0; i< nbacteria ; i++)
    {
        double ds = Cal_ChemoGradient2(i) ;
        Cal_OrientationBacteria(i) ;
        if (ds >= 0.0)
        {
            double preferedAngle = atan2(bacteria[i].nodes[(nnode-1)/2].y - domainy/2.0, bacteria[i].nodes[(nnode-1)/2].x - domainx/2.0) ;
            /*
            if (cos(bacteria[i].orientation - preferedAngle) < 0.0  )
            {
                //cout<< i<< '\t'<< bacteria[i].orientation<< '\t'<< preferedAngle<<'\t'<<cos(bacteria[i].orientation - preferedAngle)<<endl ;
            }
             */
            bacteria[i].reversalPeriod = 1.0 /( reversalRate * exp(chemoStrength * (-1.0 *fmotor  /* *cos(bacteria[i].orientation - preferedAngle )  */ )) ) ;
            
            
        }
        else
        {
            bacteria[i].reversalPeriod = 1.0/ reversalRate ;
        }
    }
}
//------------------------------------------------------------------------------------------------------------
double Cal_ChemoGradient (int i)
{
    double c0 = 1.0 ;   //amplitude of chemical
   // double ds = bacteria[i].nodes[(nnode-1)/2].x - bacteria[i].oldLoc.at(0) ;
    double r2 = sqrt( pow( (bacteria[i].nodes[(nnode-1)/2].x - domainx/2.0 ) ,2 ) + pow( (bacteria[i].nodes[(nnode-1)/2].y - domainy/2.0 ) ,2 ) );
    double r1 = sqrt( pow( (bacteria[i].oldLoc.at(0) - domainx/2.0 ) ,2 ) + pow( (bacteria[i].oldLoc.at(1) - domainy/2.0 ) ,2 ) ) ;
    double ds = -1.0* c0 * (r2 - r1) ;
    //cout<< ds<<endl ;
    return ds;
}
//------------------------------------------------------------------------------------------------------------
void WriteTrajectoryFile ()
{
    ofstream trajectories ("trajectories.txt", ofstream::app) ;
    for (uint i = 0; i< nbacteria; i++)
    {
        trajectories << bacteria[i].nodes.at((nnode-1)/2).x <<'\t'<<bacteria[i].nodes.at((nnode-1)/2).y <<'\t' ;
    }
    trajectories<< endl ;
}
//------------------------------------------------------------------------------------------------------------

void WriteNumberReverse ()
{
    ofstream histogramReversal ("HistogramReversal.txt") ;
    for (int i = 0; i < nbacteria ; i++)
    {
        histogramReversal << bacteria[i].numberReverse << '\t' << bacteria[i].timeToSource << endl ;
    }
}
//------------------------------------------------------------------------------------------------------------
double Cal_ChemoGradient2 (int i)
{
    double c0 = 1.0 ;   //amplitude of chemical
    // double ds = bacteria[i].nodes[(nnode-1)/2].x - bacteria[i].oldLoc.at(0) ;
    int tmpXIndex = static_cast<int>( round( bacteria[i].nodes[(nnode-1)/2].x / 0.5 ) ) ; //0.5 is grid size, need to change it as a variable
    int tmpYIndex = static_cast<int>( round( bacteria[i].nodes[(nnode-1)/2].y / 0.5 ) ) ;
    gridInMain.at(tmpYIndex).at(tmpXIndex) ;
    double ds = gridInMain.at(tmpYIndex).at(tmpXIndex) - bacteria[i].oldChem ;
   // cout<<ds<<endl ;
    return ds;
}
//------------------------------------------------------------------------------------------------------------
void SlimeTraceHyphae (Fungi tmpFng)
{
    int m = 0 ;
    int n = 0 ;
    double xMin;
    double xMax ;
    double yMin ;
    double yMax ;
    int mMin = 0 ;
    int mMax= 200 ;
    int nMin = 0 ;
    int nMax = 200 ;
    double tmpH ;
    double tmpS ;
    double tmpL ;
    double vec1x , vec1y, vec2x , vec2y ;
    
    
    for (uint i=0; i< tmpFng.hyphaeSegments.size(); i++)
    {
        xMin = min(tmpFng.hyphaeSegments.at(i).x1,tmpFng.hyphaeSegments.at(i).x2) ;
        xMax = max(tmpFng.hyphaeSegments.at(i).x1,tmpFng.hyphaeSegments.at(i).x2) ;
        yMin = min(tmpFng.hyphaeSegments.at(i).y1,tmpFng.hyphaeSegments.at(i).y2) ;
        yMax = max(tmpFng.hyphaeSegments.at(i).y1,tmpFng.hyphaeSegments.at(i).y2) ;
        mMin = (static_cast<int> (floor ( fmod (xMin + domainx , domainx) / dx ) ) ) % nx  ;
        mMax = (static_cast<int> (ceil ( fmod (xMax + domainx , domainx) / dx ) ) ) % nx  ;
        nMin = (static_cast<int> (floor ( fmod (yMin + domainy , domainy) / dy ) ) ) % ny  ;
        nMax = (static_cast<int> (ceil ( fmod (yMax + domainy , domainy) / dy ) ) ) % ny  ;
        
        tmpL = Dist2D(tmpFng.hyphaeSegments.at(i).x1, tmpFng.hyphaeSegments.at(i).y1, tmpFng.hyphaeSegments.at(i).x2, tmpFng.hyphaeSegments.at(i).y2) ;
        
        for (int m = mMin ; m <= mMax ; m++)
            
        {
            for (int n = nMin; n <= nMax; n++)
            {
                vec1x = tmpFng.hyphaeSegments.at(i).x1 - m * dx ;
                vec1y = tmpFng.hyphaeSegments.at(i).y1 - n * dy ;
                vec2x = tmpFng.hyphaeSegments.at(i).x2 - m * dx ;
                vec2y = tmpFng.hyphaeSegments.at(i).y2 - n * dy ;
                tmpS = TriangleArea(vec1x, vec1y, vec2x, vec2y) ;
                // 0.5 * h * l  = S
                tmpH = 2.0 * tmpS / tmpL ;
                if (tmpH < dx )
                {
                    //simple hyphae
                    slime [m][n] = 2*s0 ;
                    double decay = 10.0;
                    double tmpDis = 0.0 ;
                    
                    for (int j=0; j<tmpFng.tips.at(0).size(); j++)
                    {
                        int id = tmpFng.tipsID.at(j) ;
                        vec1x= tmpFng.hyphaeSegments.at(id).x2 - m*dx ;
                        vec1y= tmpFng.hyphaeSegments.at(id).y2 - n*dy ;
                        vec2x= tmpFng.hyphaeSegments.at(id).x2 - tmpFng.hyphaeSegments.at(id).x1 ;
                        vec2y= tmpFng.hyphaeSegments.at(id).y2 - tmpFng.hyphaeSegments.at(id).y1 ;
                        double tmplength = MagnitudeVec(vec2x, vec2y) ;
                        double tmpAngle = AngleOfTwoVectors(vec1x, vec1y, vec2x, vec2y) ;
                        tmpDis = Dist2D(tmpFng.tips.at(0).at(j), tmpFng.tips.at(1).at(j), m * dx, n * dy) ;
                        if (tmpAngle < pi/ 20.0 && tmpDis < tmplength   )
                        {
                            
                        slime[m][n] += exp(-tmpDis/ decay) ;
                            
                        }
                    }
                }
    
                
            }
            
        }
    }
}
