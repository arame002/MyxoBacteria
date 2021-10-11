

#include "TissueBacteria.hpp"
#include "ranum2.h"
//#include "SignalCommon.hpp"
#include <chrono>

long  idum=(-799);



// Inverse function is on inside the InverseTime
// Function: WriteInFile, RunMode ,


void RandomForce() ;
                                            // inverse head and tail
void Initialization () ;
void Initialization2 () ;
void CircularInitialization () ;
void SwarmingInitialization () ;
void PositionUpdating (double  ) ;                               // (time step)
void ljNodesPosition () ;
 void ParaView2 () ;
void InitialReversalTime () ;
void InitialPili () ;
void SurfaceCoverage () ;
void VisitsPerGrid () ;
int PowerLawExponent () ;


//-----------------------------------------------------------------------------------------------------
//Bacteria properties
double length = 2 ;
double B = 0.5 ;                      // bending constant
double tetta0=3.1415 ;               // prefered angle, Pi
double K =  1.0 ;                      // linear spring constant (stiffness)
double x0 = length/(nnode-1) ;                      // equilibrium distance of spring
double fmotor = 0.3 ;                           //   Ft/(N-1)
double Rp = 0.3 ;                               // protein exchange rate
double turnPeriod = 5.0 ;

double reversalRate = .25 ;
double chemoStrength = 10.0 ;

//-----------------------------------------------------------------------------------------------------
//domain and medium properties
//double kblz=1.0 , temp=1.0 , eta1 = 0.02 , eta2 = 5.0 ;
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
//vector<vector<double> > gridInMain ;

//-----------------------------------------------------------------------------------------------------
//slime properties
double sr = 0.25 ;                          // slime rate production
double kd = 0.05 ;                           // slime decay rate
//double slimeEffectiveness = 0.6 ;           // needed for slime trail following
double s0 = 1 ;
//double slime[nx][ny] ;
int surfaceCoverage[nx][ny] ;
double coveragePercentage = 0 ;
int searchAreaForSlime = static_cast<int> (round((4.0 * length)/ min(dx , dy))) ;


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
double runTime = 200.0 ;

    TissueBacteria tissueBacteria ;
    //bacterium bacteria[nbacteria];

int main ()
{

   auto start = std::chrono::high_resolution_clock::now() ;
    
    srand(time (0)) ;
    cout<<"Search area for slime is "<<searchAreaForSlime<<endl ;
    for (int m=0; m<nx; m++)
    {
        for (int n=0; n<ny; n++)
        {
           tissueBacteria.slime[m][n] = s0  ;
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
    // each frame is 0.1 second
    inverseDt = inverseDt/10 ;
   
    //Myxo() ;
    //Initialization() ;
    //Initialization2() ;
    //CircularInitialization() ;
   SwarmingInitialization () ;
    ljNodesPosition() ;
    tissueBacteria.InitialProtein() ;
    InitialReversalTime() ;
    InitialPili () ;
    for (int i=0; i<nbacteria; i++)
    {
       tissueBacteria.bacteria[i].connection[i]=0.0 ;
        for (int j=0 ; j<nnode ; j++)
        {
           tissueBacteria.bacteria[i].nodes[j].xdev = 0.0 ;
           tissueBacteria.bacteria[i].nodes[j].ydev = 0.0 ;
            
        }
    }
    Fungi fungi = driver() ;
    vector<HyphaeSegment> hyphaeSegments_main = fungi.hyphaeSegments ;
    
    vector<vector<double> > pointSources ;
    pointSources = fungi.tips ;
    fungi.WriteSourceLoc() ;
   //source for simulation 1
    //pointSources = tissueBacteria.GridSources() ;
    tissueBacteria.gridInMain = tissueBacteria.Cal_Diffusion2D(0, domainx, 0, domainy ,nx , ny ,pointSources) ;
   
    //Bacteria would try to follow hyphae as a highway
     tissueBacteria.SlimeTraceHyphae(fungi) ;
    for (int l=0; l< (nt+1); l++)
    {
        if (l < initialNt)
        {
            
            tissueBacteria.AllNodes() ;
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
            tissueBacteria.ReversalTime() ;
            tissueBacteria.AllNodes() ;
            tissueBacteria.Spring () ;
            tissueBacteria.Bending() ;
            //tissueBacteria.u_lj() ;
            //  RandomForce() ;
            tissueBacteria.Motor () ;
            tissueBacteria.TurnOrientation() ;
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
                tissueBacteria.ParaView ()  ;
                ParaView2() ;
               tissueBacteria.WriteTrajectoryFile() ;
               cout<<(l-initialNt)/inverseDt<<endl ;
                //   cout << averageLengthFree<<'\t'<<nAttachedPili<<endl ;
                //    cout << coveragePercentage <<endl ;
            }
            PositionUpdating(dt) ;
           tissueBacteria.UpdateReversalFrequency() ;       // test Effect of chemoattacrant on reversal motion
        }
        
    }
    
   tissueBacteria.WriteNumberReverse() ;
   std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
   std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
   cout << "Time to simulate "<<nbacteria<<" bacteria is : " << duration.count() << " seconds" << endl;
    cout<<"No Bug"<<endl ;
    return 0 ;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------


void PositionUpdating (double t)
{   double etaInverse = tissueBacteria.diffusion / (tissueBacteria.kblz * tissueBacteria.temp) ;
    double totalForceX = 0 ;
    double totalForceY = 0 ;
    for (int i=0; i<nbacteria; i++)
    {
       tissueBacteria.bacteria[i].oldLoc.at(0) = tissueBacteria.bacteria[i].nodes[(nnode-1)/2].x ;
       tissueBacteria.bacteria[i].oldLoc.at(1) = tissueBacteria.bacteria[i].nodes[(nnode-1)/2].y ;
       int tmpXIndex = static_cast<int>( fmod( round(fmod( tissueBacteria.bacteria[i].nodes[(nnode-1)/2].x + domainx,domainx ) / dx ), nx) ) ;
       int tmpYIndex = static_cast<int>( fmod( round(fmod( tissueBacteria.bacteria[i].nodes[(nnode-1)/2].y + domainy,domainy ) / dy ), ny )) ;
       if (tmpXIndex< 0 || tmpXIndex > nx-1 || tmpYIndex< 0 || tmpYIndex > ny -1 )
       {
          cout<<"positionUpdating, index out of range "<<tmpXIndex<<'\t'<<tmpYIndex<<endl ;
       }
       tissueBacteria.bacteria[i].oldChem = tissueBacteria.gridInMain.at(tmpYIndex).at(tmpXIndex) ;
        for(int j=0 ; j<nnode; j++)
        {
            tissueBacteria.Diffusion(tissueBacteria.bacteria[i].nodes[j].x , tissueBacteria.bacteria[i].nodes[j].y) ;
            etaInverse = tissueBacteria.diffusion / (tissueBacteria.kblz * tissueBacteria.temp) ;
            totalForceX = tissueBacteria.bacteria[i].nodes[j].fSpringx ;       // it is not += because the old value is related to another node
            totalForceX += tissueBacteria.bacteria[i].nodes[j].fBendingx ;
               //  totalForceX += tissueBacteria.bacteria[i].nodes[j].fljx ;
            totalForceX += tissueBacteria.bacteria[i].nodes[j].fMotorx ;
           tissueBacteria.bacteria[i].nodes[j].x += t * totalForceX * etaInverse ;
           tissueBacteria.bacteria[i].nodes[j].x += tissueBacteria.bacteria[i].nodes[j].xdev ;
            
            totalForceY = tissueBacteria.bacteria[i].nodes[j].fSpringy ;       // it is not += because the old value is related to another node
            totalForceY += tissueBacteria.bacteria[i].nodes[j].fBendingy ;
              //   totalForceY += tissueBacteria.bacteria[i].nodes[j].fljy ;
            totalForceY += tissueBacteria.bacteria[i].nodes[j].fMotory ;
           tissueBacteria.bacteria[i].nodes[j].y += t * totalForceY * etaInverse ;
           tissueBacteria.bacteria[i].nodes[j].y += tissueBacteria.bacteria[i].nodes[j].ydev ;
            
            
            
            // F pili
            
            if (j==0)
            {
                for (int m=0 ; m<nPili ; m++)
                {
                   tissueBacteria.bacteria[i].nodes[j].x += t * (tissueBacteria.bacteria[i].pili[m].fx) * etaInverse ;
                   tissueBacteria.bacteria[i].nodes[j].y += t * (tissueBacteria.bacteria[i].pili[m].fy) * etaInverse ;
                }
            }
            
        }
    }
    ljNodesPosition() ;
}

//-----------------------------------------------------------------------------------------------------
double Distance (int i,int j, int k, int l)
{
    double dis= sqrt ((tissueBacteria.bacteria[i].nodes[j].x - tissueBacteria.bacteria[k].nodes[l].x) * (tissueBacteria.bacteria[i].nodes[j].x- tissueBacteria.bacteria[k].nodes[l].x) + (tissueBacteria.bacteria[i].nodes[j].y - tissueBacteria.bacteria[k].nodes[l].y) * (tissueBacteria.bacteria[i].nodes[j].y - tissueBacteria.bacteria[k].nodes[l].y)) ;
    return dis ;
}
//-----------------------------------------------------------------------------------------------------
double Cos0ijk(int i,int m)
{
    double dot= (tissueBacteria.bacteria[i].nodes[m-1].x - tissueBacteria.bacteria[i].nodes[m].x)*(tissueBacteria.bacteria[i].nodes[m+1].x - tissueBacteria.bacteria[i].nodes[m].x);
    dot += (tissueBacteria.bacteria[i].nodes[m-1].y - tissueBacteria.bacteria[i].nodes[m].y)*(tissueBacteria.bacteria[i].nodes[m+1].y - tissueBacteria.bacteria[i].nodes[m].y) ;
    double cos= dot / ( Distance(i,m-1,i,m) * Distance(i,m+1,i,m) ) ;
    return cos ;
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
           tissueBacteria.bacteria[i].nodes[j].x = lx * coloum ;
        }
        else   { tissueBacteria.bacteria[i].nodes[j].x = lx * coloum + lx /2 ; }
       tissueBacteria.bacteria[i].nodes[j].y=  ly * row /2 ;
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
           tissueBacteria.bacteria[i].nodes[j+n].x = tissueBacteria.bacteria[i].nodes[j+n-1].x + x0 * a ; // or nodes[j].x + n * x0 * a
           tissueBacteria.bacteria[i].nodes[j-n].x = tissueBacteria.bacteria[i].nodes[j-n+1].x - x0 * a ;
           tissueBacteria.bacteria[i].nodes[j+n].y = tissueBacteria.bacteria[i].nodes[j+n-1].y + x0 * b ;
           tissueBacteria.bacteria[i].nodes[j-n].y = tissueBacteria.bacteria[i].nodes[j-n+1].y - x0 * b ;
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
    int nextLine = static_cast<int>(round( sqrt(nbacteria/2.0)/2.0) )  ;
    cout<< "next line is: "<<nextLine<< endl ;
    int row = -2.0 * nextLine ;
    int coloum = -nextLine  ;
    double tmpDomainX = domainx * 0.8 ;
    double tmpDomainY = domainy * 0.8 ;
    double a ;
    double b ;
    double lx = sqrt (2.0 * tmpDomainX * tmpDomainY / nbacteria ) ;
    double ly = lx ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        if (row%2==0)
        {
           tissueBacteria.bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum ;
        }
        else   { tissueBacteria.bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum + lx /2.0 ; }
       tissueBacteria.bacteria[i].nodes[j].y=  domainy/2.0 + ly * row /2.0 ;
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
           tissueBacteria.bacteria[i].nodes[j+n].x = tissueBacteria.bacteria[i].nodes[j+n-1].x + x0 * a ; // or nodes[j].x + n * x0 * a
           tissueBacteria.bacteria[i].nodes[j-n].x = tissueBacteria.bacteria[i].nodes[j-n+1].x - x0 * a ;
           tissueBacteria.bacteria[i].nodes[j+n].y = tissueBacteria.bacteria[i].nodes[j+n-1].y + x0 * b ;
           tissueBacteria.bacteria[i].nodes[j-n].y = tissueBacteria.bacteria[i].nodes[j-n+1].y - x0 * b ;
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
       tissueBacteria.bacteria[i].nodes[j].x = cntrX + raduis * cos( static_cast<double>(i* deltaTetta) ) ;
       tissueBacteria.bacteria[i].nodes[j].y = cntrY + raduis * sin( static_cast<double>(i* deltaTetta) ) ;
        
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
           tissueBacteria.bacteria[i].nodes[j+n].x = tissueBacteria.bacteria[i].nodes[j+n-1].x + x0 * a ; // or nodes[j].x + n * x0 * a
           tissueBacteria.bacteria[i].nodes[j-n].x = tissueBacteria.bacteria[i].nodes[j-n+1].x - x0 * a ;
           tissueBacteria.bacteria[i].nodes[j+n].y = tissueBacteria.bacteria[i].nodes[j+n-1].y + x0 * b ;
           tissueBacteria.bacteria[i].nodes[j-n].y = tissueBacteria.bacteria[i].nodes[j-n+1].y - x0 * b ;
        }
    }
    
}

void SwarmingInitialization ()
{
   
   double a ;
   double b ;
   double minX =  domainx - 11.0 * tissueBacteria.agarThicknessX ;    //3
   double maxX = domainx - 10.0 * tissueBacteria.agarThicknessX ;     //2
   double tmpDomainSizeX = maxX - minX ;
   double minY = 2.0 * tissueBacteria.agarThicknessY ;
   double maxY = domainy - 2.0 * tissueBacteria.agarThicknessY ;
   double tmpDomainSizeY = maxY - minY ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
       a = (rand() / (RAND_MAX + 1.0)) * tmpDomainSizeX ;
       b = (rand() / (RAND_MAX + 1.0)) * tmpDomainSizeY ;
       tissueBacteria.bacteria[i].nodes[j].x = minX + a ;
       tissueBacteria.bacteria[i].nodes[j].y = minY + b ;
        
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
           // or nodes[j].x + n * x0 * a
           tissueBacteria.bacteria[i].nodes[j+n].x = tissueBacteria.bacteria[i].nodes[j+n-1].x + x0 * a ;
           tissueBacteria.bacteria[i].nodes[j-n].x = tissueBacteria.bacteria[i].nodes[j-n+1].x - x0 * a ;
           tissueBacteria.bacteria[i].nodes[j+n].y = tissueBacteria.bacteria[i].nodes[j+n-1].y + x0 * b ;
           tissueBacteria.bacteria[i].nodes[j-n].y = tissueBacteria.bacteria[i].nodes[j-n+1].y - x0 * b ;
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
           tissueBacteria.bacteria[i].ljnodes[j].x = (tissueBacteria.bacteria[i].nodes[j].x + tissueBacteria.bacteria[i].nodes[j+1].x)/2 ;
           tissueBacteria.bacteria[i].ljnodes[j].y = (tissueBacteria.bacteria[i].nodes[j].y + tissueBacteria.bacteria[i].nodes[j+1].y)/2 ;
        }
    }
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
    SignalOut << "SCALARS liquid float 1" << endl;
    SignalOut << "LOOKUP_TABLE default" << endl;
    
    for (int k = 0; k < nz ; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                SignalOut << tissueBacteria.slime[i][j] << endl;
            }
        }
    }
    
    index1++ ;
}

//-----------------------------------------------------------------------------------------------------
void InitialReversalTime ()
{
    for( int i=0 ; i<nbacteria ; i++)
    {
       tissueBacteria.bacteria[i].reversalTime = (rand() / (RAND_MAX + 1.0)) * tissueBacteria.bacteria[i].reversalPeriod ;
       tissueBacteria.bacteria[i].reversalTime -= fmod(tissueBacteria.bacteria[i].reversalTime , dt ) ;
        if (rand() / (RAND_MAX + 1.0) < 0.5)
        {
            //bacteria[i].directionOfMotion = false ;
           tissueBacteria.bacteria[i].directionOfMotion = true ;
        }
        else
        {
           tissueBacteria.bacteria[i].directionOfMotion = true ;
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------
void InitialPili ()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for ( int j=0 ; j<nPili ; j++)
        {
           tissueBacteria.bacteria[i].pili[j].lFree = rand() / (RAND_MAX + 1.0) * piliMaxLength ;
           tissueBacteria.bacteria[i].pili[j].attachment = false ;
           tissueBacteria.bacteria[i].pili[j].retraction = false ;
           tissueBacteria.bacteria[i].pili[j].fx = 0.0 ;
           tissueBacteria.bacteria[i].pili[j].fy = 0.0 ;
           tissueBacteria.bacteria[i].pili[j].F = 0.0 ;
            
            if(rand() / (RAND_MAX + 1.0) < 0.5 )
            {
               tissueBacteria.bacteria[i].pili[j].retraction = true ;
            }
            else
            {
               tissueBacteria.bacteria[i].pili[j].retraction = false ;
            }
            if (rand() / (RAND_MAX + 1.0) < 0.5 )
            {
               tissueBacteria.bacteria[i].pili[j].attachment = true ;
               tissueBacteria.bacteria[i].pili[j].retraction = true ;
            }
            else
            {
               tissueBacteria.bacteria[i].pili[j].attachment = false ;
            }
            
        }
    }
}

//-----------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
void SurfaceCoverage ()
{
    double percentage = 0 ;
    
    for (int m=0; m< nx; m++)
    {
        for (int n=0 ; n < ny ; n++)
        {
            if (tissueBacteria.slime[m][n] > s0 && surfaceCoverage[m][n] == 0 )
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
            m = (static_cast<int> (round ( fmod (tissueBacteria.bacteria[i].allnodes[j].x + domainx , domainx) / dx ) ) ) % nx  ;
            n = (static_cast<int> (round ( fmod (tissueBacteria.bacteria[i].allnodes[j].y + domainy , domainy) / dy ) ) ) % ny  ;
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

