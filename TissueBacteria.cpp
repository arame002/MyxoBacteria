//
//  TissueBacteria.cpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 10/8/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#include "TissueBacteria.hpp"

using constants::pi;

TissueBacteria:: TissueBacteria()
{
    slime.resize(nx,vector<double> (ny) ) ;
    viscousDamp.resize(nx, vector<double> (ny) ) ;
    sourceChemo.resize(2) ;
    X.resize(nx) ;
    Y.resize(ny) ;
    Z.push_back(0.0) ;
    visit.resize(nx, vector<int> (ny, 0.0) ) ;
    surfaceCoverage.resize(nx, vector<int> (ny, 0.0) ) ;
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria::Initialze_AllRandomForce()
{
    for (int i=0; i<nbacteria; i++)
    {
        bacteria[i].connectionToOtherB[i]=0.0 ;
        bacteria[i].initialize_randomForce() ;
    }
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria::UpdateTissue_FromConfigFile()
{
    //Input and output parameters
    folderName = globalConfigVars.getConfigValue("AnimationFolder").toString() ;
    statsFolder = globalConfigVars.getConfigValue("StatFolderName").toString() ;
    animationName = globalConfigVars.getConfigValue("AnimationName").toString() ;
    
    //Timimg control parametes
    initialTime = globalConfigVars.getConfigValue("InitTimeStage").toDouble() ;
    runTime = globalConfigVars.getConfigValue("SimulationTotalTime").toDouble() ;
    dt = globalConfigVars.getConfigValue("SimulationTimeStep").toDouble() ;
    
    //Bacteria mechanical parameters
    length = globalConfigVars.getConfigValue("BacteiaLength").toDouble() ;
    B = globalConfigVars.getConfigValue("BacteriaBenCoeff").toDouble() ;
    K = globalConfigVars.getConfigValue("BacteriaLinearStiff").toDouble() ;
    x0 = length/(nnode - 1) ;
    
    //Lennard-Jones interaction parameters
    lj_Energy = globalConfigVars.getConfigValue("ELJ").toDouble() ;
    lj_rMin = globalConfigVars.getConfigValue("rMinLJ").toDouble() ;
    lj_Interaction = static_cast<bool>(globalConfigVars.getConfigValue("interactingLJ").toInt() ) ;
    
    //Bacteria motor parametes
    fmotor = globalConfigVars.getConfigValue("Bacteira_motorForce").toDouble() ;
    motorEfficiency_Liquid = globalConfigVars.getConfigValue("Bacteria_motorEfficiency_Liquid").toDouble() ;
    motorEfficiency_Agar = globalConfigVars.getConfigValue("Bacteria_motorEfficiency_Agar").toDouble() ;
    motorEfficiency = motorEfficiency_Liquid ;
    
    //Bacteria reversing parametes
    turnPeriod = globalConfigVars.getConfigValue("Bacteria_turnPeriod").toDouble() ;
    minimumRunTime = globalConfigVars.getConfigValue("Bacteria_minRunTime").toDouble() ;
    reversalRate = globalConfigVars.getConfigValue("Bacteria_reversalRate").toDouble() ;
    chemoStrength_run = globalConfigVars.getConfigValue("Bacteria_chemotaxisSensitivity_run").toDouble() ;
    chemoStrength_wrap = globalConfigVars.getConfigValue("Bacteria_chemotaxisSensitivity_wrap").toDouble() ;
    
    
    //Medium related parameters
    domainx =   globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble() -
                globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble() ;
    domainy =   globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble() -
                globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble() ;
    eta_background = globalConfigVars.getConfigValue("damping1").toDouble() ;
    eta_Barrier = globalConfigVars.getConfigValue("damping2").toDouble() ;
    agarThicknessX = domainx * globalConfigVars.getConfigValue("boundary_agarThicknessX_Coeff").toDouble() ;
    agarThicknessY = domainy * globalConfigVars.getConfigValue("boundary_agarThicknessY_Coeff").toDouble() ;
    
    //Slime( liquid) related parameters
    sr = globalConfigVars.getConfigValue("slime_SecretionRate").toDouble() ;
    kd = globalConfigVars.getConfigValue("slime_DecayRate").toDouble() ;
    slimeEffectiveness = globalConfigVars.getConfigValue("slime_Effectiveness").toDouble() ;
    liqBackground = globalConfigVars.getConfigValue("slime_background").toDouble() ;
    liqLayer = globalConfigVars.getConfigValue("slime_Layer").toDouble() ;
    liqHyphae = globalConfigVars.getConfigValue("slime_Hyphae").toDouble() ;
    Slime_CutOff = globalConfigVars.getConfigValue("slime_Attachment_CutOff").toDouble() ;
    searchAreaForSlime = static_cast<int> (round(( length )/ min(dx , dy))) ;
    
    inLiquid = globalConfigVars.getConfigValue("Bacteria_inLiquid").toInt() ;
    PBC = globalConfigVars.getConfigValue("Bacteria_PBC").toInt() ;
    chemotacticMechanism = static_cast<ChemotacticMechanism>( globalConfigVars.getConfigValue("Bacteria_chemotacticModel").toInt() ) ;
    initialCondition =static_cast<InitialCondition>( globalConfigVars.getConfigValue("Bacteria_InitialCondition").toInt() ) ;
    run_calibrated = static_cast<bool>( globalConfigVars.getConfigValue("Run_Calibrated").toInt() ) ;
    
    for (int i = 0; i< nbacteria; i++)
    {
        bacteria[i].UpdateBacteria_FromConfigFile() ;
    }
    tGrids.UpdateTGrid_FromConfigFile() ;
    
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria::Bacteria_Initialization()
{
    if (initialCondition == uniform)
    {
        Initialization2() ;
    }
    else if (initialCondition == circular)
    {
        CircularInitialization() ;
    }
    else if (initialCondition == swarm)
    {
        SwarmingInitialization () ;
    }
    else if (initialCondition == center)
    {
        CenterInitialization() ;
    }
    else if (initialCondition == alongNetwork)
    {
        AlongNetworkInitialization () ;
    }
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria::Initialization ()
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
void TissueBacteria::Initialization2 ()
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
           bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum ;
        }
        else   { bacteria[i].nodes[j].x = domainx/2.0 + lx * coloum + lx /2.0 ; }
       bacteria[i].nodes[j].y =  domainy/2.0 + ly * row /2.0 ;
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
void TissueBacteria::CircularInitialization ()
{
    double raduis = 0.75 * sqrt(domainx * domainx  )/2.0 ;
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
void TissueBacteria::SwarmingInitialization ()
{
   
   double a ;
   double b ;
   double minX =  domainx - 4.0 * agarThicknessX ;    //3
   double maxX = domainx - 3.0 * agarThicknessX ;     //2
   double tmpDomainSizeX = maxX - minX ;
   double minY = 2.0 * agarThicknessY ;               //2
   double maxY = domainy - 2.0 * agarThicknessY ;     //2
   double tmpDomainSizeY = maxY - minY ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
       //a = (rand() / (RAND_MAX + 1.0)) * tmpDomainSizeX ;
       a = fmod(i, 2) * tmpDomainSizeX ;
       //b = (rand() / (RAND_MAX + 1.0)) * tmpDomainSizeY ;
       b = static_cast<double>(i)/ static_cast<double>(nbacteria) * tmpDomainSizeY  ;
       bacteria[i].nodes[j].x = minX + a ;
       bacteria[i].nodes[j].y = minY + b ;
        
        a= gasdev(&idum) ;
        b=gasdev(&idum) ;
        a= a/ sqrt(a*a+b*b) ;
        b = b/ sqrt(a*a+b*b) ;
        for (int n=1; n<=(nnode-1)/2; n++)
        {
           // or nodes[j].x + n * x0 * a
           bacteria[i].nodes[j+n].x = bacteria[i].nodes[j+n-1].x + x0 * a ;
           bacteria[i].nodes[j-n].x = bacteria[i].nodes[j-n+1].x - x0 * a ;
           bacteria[i].nodes[j+n].y = bacteria[i].nodes[j+n-1].y + x0 * b ;
           bacteria[i].nodes[j-n].y = bacteria[i].nodes[j-n+1].y - x0 * b ;
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria::CenterInitialization()
{
    double raduis = 0.2 * sqrt(domainx * domainx  )/2.0 ;
    double cntrX = domainx/2.0 ;
    double cntrY = domainy/2.0 ;
    
    double a ;
    double b ;
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        bacteria[i].nodes[j].x = cntrX + raduis * (2.0 * rand() / (RAND_MAX + 1.0) - 1.0 )  ;
        bacteria[i].nodes[j].y = cntrY + raduis * (2.0 * rand() / (RAND_MAX + 1.0) - 1.0 ) ;
        
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

void TissueBacteria::AlongNetworkInitialization()
{

    double tmpX ;
    double tmpY ;
    int m_x ;
    int n_y ;
    double a ;
    double b ;
    
    for (int i=0 , j=(nnode-1)/2 ; i<nbacteria; i++)
    {
        do {
            tmpX = (rand() / (RAND_MAX + 1.0) ) * domainx ;
            tmpY = (rand() / (RAND_MAX + 1.0) ) * domainy ;
            m_x = (static_cast<int> (round ( fmod (tmpX + domainx , domainx) / dx ) ) ) % nx ;
            n_y = (static_cast<int> (round ( fmod (tmpY + domainy , domainy) / dy ) ) ) % ny ;
        } while (slime[m_x][n_y] != liqLayer  );
        
        bacteria[i].nodes[j].x = tmpX ;
        bacteria[i].nodes[j].y = tmpY ;
        
        //cout<< tmpX<<'\t'<<tmpY<< '\t'<< m_x<< '\t'<< n_y<<endl ;
        
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

void TissueBacteria::Update_LJ_NodePositions ()
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
void TissueBacteria::Initialize_ReversalTimes ()
{
    for( int i=0 ; i<nbacteria ; i++)
    {
       bacteria[i].reversalPeriod =  bacteria[i].maxRunDuration ;
       bacteria[i].reversalTime = (rand() / (RAND_MAX + 1.0)) * bacteria[i].maxRunDuration ;
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
    /*
   for( int i=0 ; i<nbacteria ; i++)
   {
      cout<< "initial reversal time is "<<bacteria[i].reversalTime << endl ;
   }
     */
}

//-----------------------------------------------------------------------------------------------------
vector<vector<double> > TissueBacteria::TB_Cal_ChemoDiffusion2D(double xMin, double xMax, double yMin, double yMax,int nGridX , int nGridY,vector<vector<double> > sources, vector<double> pSource)
{
    tGrids = Diffusion2D(xMin, xMax, yMin, yMax,nGridX , nGridY ,sources, pSource, tGrids ) ;
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
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: Cal_AllLinearSpring_Forces()
{
    for (int i=0 ; i<nbacteria ; i++)
    {
        for (int j=0; j<nnode; j++)
        {
            if (j==0 )
            {
                bacteria[i].nodes[j].fSpringx= -K * (Cal_NodeToNodeDistance (i,j,i,j+1 ) - x0) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x) / Cal_NodeToNodeDistance (i,j,i,j+1) ;
                
                bacteria[i].nodes[j].fSpringy= -K * (Cal_NodeToNodeDistance (i,j,i,j+1 ) - x0) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y) / Cal_NodeToNodeDistance (i,j,i,j+1) ;
            }
            
            else if (j== (nnode-1))
            {
                bacteria[i].nodes[j].fSpringx = -K * (Cal_NodeToNodeDistance (i,j,i,j-1) - x0) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j-1].x) / Cal_NodeToNodeDistance (i,j,i,j-1) ;
                
                bacteria[i].nodes[j].fSpringy = -K * (Cal_NodeToNodeDistance (i,j,i,j-1) - x0) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j-1].y) / Cal_NodeToNodeDistance (i,j,i,j-1) ;
            }
            else
            {
                bacteria[i].nodes[j].fSpringx  = -K * ( Cal_NodeToNodeDistance(i,j,i,j-1) - x0 ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j-1].x) /Cal_NodeToNodeDistance(i,j,i,j-1) ;
                bacteria[i].nodes[j].fSpringx += -K * ( Cal_NodeToNodeDistance(i,j,i,j+1) - x0 ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x) /Cal_NodeToNodeDistance(i,j,i,j+1) ;
                
                bacteria[i].nodes[j].fSpringy  = -K * ( Cal_NodeToNodeDistance(i,j,i,j-1) - x0 ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j-1].y) /Cal_NodeToNodeDistance(i,j,i,j-1) ;
                bacteria[i].nodes[j].fSpringy += -K * ( Cal_NodeToNodeDistance(i,j,i,j+1) - x0 ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y) /Cal_NodeToNodeDistance(i,j,i,j+1) ;
            }
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------
double TissueBacteria:: Cal_NodeToNodeDistance (int i,int j, int k, int l)
{
    double dis= sqrt ((bacteria[i].nodes[j].x - bacteria[k].nodes[l].x) * (bacteria[i].nodes[j].x-bacteria[k].nodes[l].x) + (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y) * (bacteria[i].nodes[j].y - bacteria[k].nodes[l].y)) ;
    return dis ;
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria:: Cal_AllBendingSpring_Forces()
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
            Bend1x[j] += -B *(bacteria[i].nodes[j+2].x-bacteria[i].nodes[j+1].x)/(Cal_NodeToNodeDistance(i,j,i,j+1) * Cal_NodeToNodeDistance(i,j+2,i,j+1)) ;
            Bend1x[j] +=  B *(Cal_Cos0ijk_InteriorNodes(i,j+1)/(Cal_NodeToNodeDistance(i,j,i,j+1)*Cal_NodeToNodeDistance(i,j,i,j+1))*(bacteria[i].nodes[j].x-bacteria[i].nodes[j+1].x)) ;
            
            Bend1y[j] += -B *(bacteria[i].nodes[j+2].y-bacteria[i].nodes[j+1].y)/(Cal_NodeToNodeDistance(i,j,i,j+1) * Cal_NodeToNodeDistance(i,j+2,i,j+1)) ;
            Bend1y[j] +=  B *(Cal_Cos0ijk_InteriorNodes(i,j+1)/(Cal_NodeToNodeDistance(i,j,i,j+1)*Cal_NodeToNodeDistance(i,j,i,j+1))*(bacteria[i].nodes[j].y-bacteria[i].nodes[j+1].y)) ;
            
        }
        for (int j=2; j<nnode ; j++)
        {
            Bend3x[j] += -B *(bacteria[i].nodes[j-2].x-bacteria[i].nodes[j-1].x)/(Cal_NodeToNodeDistance(i,j-2,i,j-1) * Cal_NodeToNodeDistance(i,j,i,j-1)) ;
            Bend3x[j] +=  B *(Cal_Cos0ijk_InteriorNodes(i,j-1)/(Cal_NodeToNodeDistance(i,j,i,j-1)*Cal_NodeToNodeDistance(i,j,i,j-1))*(bacteria[i].nodes[j].x-bacteria[i].nodes[j-1].x)) ;
            Bend3y[j] += -B *(bacteria[i].nodes[j-2].y-bacteria[i].nodes[j-1].y)/(Cal_NodeToNodeDistance(i,j-2,i,j-1) * Cal_NodeToNodeDistance(i,j,i,j-1)) ;
            Bend3y[j] +=  B *(Cal_Cos0ijk_InteriorNodes(i,j-1)/(Cal_NodeToNodeDistance(i,j,i,j-1)*Cal_NodeToNodeDistance(i,j,i,j-1))*(bacteria[i].nodes[j].y-bacteria[i].nodes[j-1].y)) ;
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
double TissueBacteria:: Cal_Cos0ijk_InteriorNodes(int i,int m)
{
    double dot= (bacteria[i].nodes[m-1].x - bacteria[i].nodes[m].x)*(bacteria[i].nodes[m+1].x - bacteria[i].nodes[m].x);
    dot += (bacteria[i].nodes[m-1].y - bacteria[i].nodes[m].y)*(bacteria[i].nodes[m+1].y -bacteria[i].nodes[m].y) ;
    double cos= dot / ( Cal_NodeToNodeDistance(i,m-1,i,m) * Cal_NodeToNodeDistance(i,m+1,i,m) ) ;
    return cos ;
}
//-----------------------------------------------------------------------------------------------------
double TissueBacteria:: Cal_PtoP_Distance_PBC (double x1, double y1, double x2, double y2, int nx , int ny)
{
    double dis= sqrt((x1+nx*domainx -x2)*(x1+nx*domainx -x2)+(y1+ny*domainy -y2)*(y1+ny*domainy -y2)) ;
    return dis ;
}

//-----------------------------------------------------------------------------------------------------
double TissueBacteria::Cal_MinDistance_PBC (double x1 , double y1 ,double x2 , double y2)
{   shiftx = 0 ;
    shifty = 0 ;
    double min = Cal_PtoP_Distance_PBC(x1, y1, x2, y2 ,0 , 0) ;
    double a ;
    
    for (int nx = -1 ; nx<2; nx++)
    {
        for (int ny= -1 ; ny<2; ny++)
        {   a = Cal_PtoP_Distance_PBC(x1, y1, x2, y2,nx,ny) ;
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

void TissueBacteria::Update_BacterialConnection ()
{
    double xmin = sqrt((x0/4)*(x0/4)+ 0.36) ;
    double d  ;
    double min ;
    for (int i=0; i<nbacteria-1; i++)
    {
        for (int j=i+1; j<nbacteria; j++)
        {
            bacteria[i].connectionToOtherB[j] = 0.0 ;
            bacteria[j].connectionToOtherB[i] = 0.0 ;
            
            min = Cal_MinDistance_PBC(bacteria[i].nodes[(nnode-1)/2].x , bacteria[i].nodes[(nnode-1)/2].y , bacteria[j].nodes[(nnode-1)/2].x , bacteria[j].nodes[(nnode-1)/2].y)  ;
            
            if (min < (1.2 * length) )
            {
                for (int m=0; m<(2*nnode-1); m++)
                {
                    for (int n=0; n<(2*nnode-1); n++)
                    {   d = Cal_PtoP_Distance_PBC ( bacteria[i].allnodes[m].x , bacteria[i].allnodes[m].y, bacteria[j].allnodes[n].x , bacteria[j].allnodes[n].y ,shiftx , shifty ) ;
                        if ( d < xmin )
                        { //  bacteria[i].connectionToOtherB[j] += 0.5 ;
                            //  bacteria[j].connectionToOtherB[i] += 0.5 ;
                            if (m==0 || n==0 || m== (2*nnode-2)|| n== (2*nnode-2) )
                            {
                                bacteria[i].connectionToOtherB[j] += 0.25 ;
                                bacteria[j].connectionToOtherB[i] += 0.25 ;
                            }
                            else
                            {
                                bacteria[i].connectionToOtherB[j] += 0.5 ;
                                bacteria[j].connectionToOtherB[i] += 0.5 ;
                            }
                        }
                        
                        
                    }
                }
                
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------
double TissueBacteria:: Cal_AllBacteriaLJ_Forces()
{
    double delta_x, delta_y, rcut2;
    for (int i=0; i<nbacteria; i++)
    {   for(int j=0 ; j< nnode ; j++)
    {
        bacteria[i].nodes[j].fljx=0.0;
        bacteria[i].nodes[j].fljy=0.0;
    }
    }
    double sigma2= (lj_rMin/1.122)*(lj_rMin/1.222) ;          // Rmin = 1.122 * sigma
    double ulj=0.0 ;
    rcut2= 6.25 * sigma2  ;                 // rcut=2.5*sigma for Lenard-Jones potential
    double min ;
    if (lj_Interaction == false)
    {
        return ulj ;
    }
    
    for (int i=0; i<(nbacteria -1); i++)
    {
        for (int j=i+1 ; j<nbacteria ; j++)
        {  min = Cal_MinDistance_PBC(bacteria[i].nodes[(nnode-1)/2].x , bacteria[i].nodes[(nnode-1)/2].y , bacteria[j].nodes[(nnode-1)/2].x , bacteria[j].nodes[(nnode-1)/2].y)  ;
            
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
                            //double force = 48.0 * eps * r6i * (r6i-0.5) * r2i ;
                            double force = 48.0 * lj_Energy * r6i * (r6i /*- 0.5*/ ) / sqrt(delta2) ;
                            ulj = ulj + 4.0 * lj_Energy * r6i * (r6i - 1) ;
                            if (m%2==0)
                            {
                               // bacteria[i].nodes[m/2].fljx += force * delta_x ;
                               // bacteria[i].nodes[m/2].fljy += force * delta_y;
                                
                                bacteria[i].nodes[m/2].fljx += force * delta_x / sqrt(delta2) ;
                                bacteria[i].nodes[m/2].fljy += force * delta_y / sqrt(delta2) ;
                            }
                            if (n%2 == 0)
                            {
                              //  bacteria[j].nodes[n/2].fljx -= force * delta_x  ;
                              //  bacteria[j].nodes[n/2].fljy -= force * delta_y  ;
                                
                                bacteria[j].nodes[n/2].fljx -= force * delta_x / sqrt(delta2) ;
                                bacteria[j].nodes[n/2].fljy -= force * delta_y / sqrt(delta2) ;
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
void TissueBacteria:: TermalFluctiation_Forces()
{
    double sig ;
    double vis ;
    for(int i=0; i<nbacteria ; i++)
    {
        for(int j=0 ; j<nnode; j++)
        {
            vis = Update_LocalFriction(bacteria[i].nodes[j].x, bacteria[i].nodes[j].y) ;
            sig = sqrt(2.0*dt/vis) ;
            bacteria[i].nodes[j].xdev =gasdev()*sig;
            bacteria[i].nodes[j].ydev =gasdev()*sig;
        }
    }
}
*/
//-----------------------------------------------------------------------------------------------------
double TissueBacteria:: Update_LocalFriction (double x , double y)
{
    int m = (static_cast<int> (round ( fmod (x + domainx , domainx) / dx ) ) ) % nx ;
    int n = (static_cast<int> (round ( fmod (y + domainy , domainy) / dy ) ) ) % ny ;
    double vis ;
    return viscousDamp[m][n] ;
    //motorEfficiency no longer in use.
    
    /*
    if (PBC == false)
    {
        if (x< agarThicknessX || x> domainx - agarThicknessX
            // || y< agarThicknessY || y> domainy - agarThicknessY
            )
        {
            vis = eta_Barrier ;
            motorEfficiency = motorEfficiency_Agar ;     //0.4
        }
        else
        {
            vis = eta_background ;
            motorEfficiency = motorEfficiency_Liquid ;
        }
    }
    else
    {
        if (inLiquid == true)
        {
            vis = eta_background ;
            motorEfficiency = motorEfficiency_Liquid ;
        }
        else
        {
            vis = (eta_background/ slime[m][n] ) ;
            //vis = eta1 ;
            motorEfficiency = motorEfficiency_Agar ;
        }
    }
     */
}
//-----------------------------------------------------------------------------------------------------
//Calculate the damping coefficient based on liquid value and boundary conditions
void TissueBacteria::Update_ViscousDampingCoeff()
{
    for (int m=0; m< nx; m++)
    {
        for (int n=0 ; n < ny ; n++)
        {
            if (inLiquid == true)
            {
                viscousDamp[m][n] = eta_background * (liqBackground / liqLayer ) ;
                if ( (m*dx < agarThicknessX || m*dx > domainx - agarThicknessX /* || n*dy< agarThicknessY || n*dy > domainy - agarThicknessY */ ) && PBC == false )
                {
                    viscousDamp[m][n] = eta_Barrier ;
                }
            }
            else
            {
                viscousDamp[m][n] = eta_background* (liqBackground / slime[m][n])  ;
                if ( (m*dx < agarThicknessX || m*dx > domainx - agarThicknessX /* || n*dy< agarThicknessY || n*dy > domainy - agarThicknessY */ ) && PBC == false )
                {
                    viscousDamp[m][n] = eta_Barrier ;
                }
            }
        }
    }
}



//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: Cal_MotorForce()
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
        double tmpFmotor = fmotor * motorEfficiency ;
        if (bacteria[i].wrapMode)
        {
            tmpFmotor *= bacteria[i].wrapSlowDown ;
        }
        /*
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
         */
        for(int j=0 ; j<nnode; j++)
        {
            if(j==0)
            {
                double fs = 0.0 ;
                
                if (inLiquid == false)
                {
                    //fungi is modeled by slime
                    //Following Slime or hyphae
                    fs = Bacterial_HighwayFollowing(i, 2.0 * length) ;        // length/2
                }
                
                bacteria[i].nodes[j].fMotorx =  (tmpFmotor - fs ) * (bacteria[i].nodes[j].x - bacteria[i].nodes[j+1].x)/ Cal_NodeToNodeDistance(i,j,i,j+1) ;                          // force to the direction of head
                
                bacteria[i].nodes[j].fMotory =  (tmpFmotor - fs ) * (bacteria[i].nodes[j].y - bacteria[i].nodes[j+1].y)/ Cal_NodeToNodeDistance(i,j,i,j+1) ;                          // force to the direction of head
                
                bacteria[i].nodes[j].fMotorx += bacteria[i].fSTFx ;     // Fh + Fs
                bacteria[i].nodes[j].fMotory += bacteria[i].fSTFy ;
            }
            else
            {
                bacteria[i].nodes[j].fMotorx = tmpFmotor * (bacteria[i].nodes[j-1].x - bacteria[i].nodes[j].x)/Cal_NodeToNodeDistance(i,j-1,i,j) ;
                bacteria[i].nodes[j].fMotory = tmpFmotor * (bacteria[i].nodes[j-1].y - bacteria[i].nodes[j].y)/Cal_NodeToNodeDistance(i,j-1,i,j) ;
            }
        }
        int m = (static_cast<int> (round ( fmod (bacteria[i].nodes[(nnode-1)/2].x + domainx , domainx) / dx ) ) ) % nx  ;
        int n = (static_cast<int> (round ( fmod (bacteria[i].nodes[(nnode-1)/2].y + domainy , domainy) / dy ) ) ) % ny  ;
        
        // Thre is more liquid near the fungi network. Slime_CutOff is not calibrated
        if (slime[m][n] > Slime_CutOff * liqBackground )
        {
            bacteria[i].attachedToFungi = true ;
            
        }
        else if (bacteria[i].attachedToFungi == true)
        {
            bacteria[i].attachedToFungi = false ;
            
        }
         
        
    }
}

//-----------------------------------------------------------------------------------------------------
// Break down the environment in front of bacteria into 5 sectors.
// Rotate the motor force of the frontier node to the region with  alteast 80% of max liquid with minumum angle change
double TissueBacteria:: Bacterial_HighwayFollowing (int i, double R)        // i th bacteria, radius of the search area
{
    
    double tmpFmotor = fmotor * motorEfficiency ;
    if (bacteria[i].wrapMode)
    {
        tmpFmotor *= bacteria[i].wrapSlowDown ;
    }
    int m=0   ;
    int n=0   ;
    int gridX ;
    int gridY ;
    double deltax ;
    double deltay ;
    double r ;
    int region = 2 ;
    double lowThreshold = 0.8 ;                     // 0.8
    double smallValue = 0.0000001 ;
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
    orientationBacteria = atan2(ay, ax) * 180 / 3.1415 ;
    
    double alfa ;
    
    m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[0].x + domainx , domainx) / dx ) ) ) % nx  ;
    n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[0].y + domainy , domainy) / dy ) ) ) % ny  ;
    for (int sx = -1*searchAreaForSlime ; sx <= searchAreaForSlime ; sx++)
    {
        for (int sy = -1*searchAreaForSlime ; sy<= searchAreaForSlime ; sy++ )
        {
            //nx added to gridX to avoid termitaing
            gridX = (m+sx + nx) % nx ;
            gridY = (n+sy + ny) % ny ;
            if (gridX<0 || gridX> nx-1 || gridY<0 || gridY > ny-1 )
            {
                cout<< m<< '\t'<<sx<<'\t'<<nx<<endl ;
               // gridX = (m+sx+nx) % nx ;
               // gridY = (n+sy+ny) % ny ;
            }
            ax = sx * dx + dx/2.0 ;
            ay = sy * dy + dy/2.0 ;
            Cos = ax/ sqrt(ax * ax + ay * ay) ;
            alfa = acos(Cos)* 180.0 / 3.1415 ;
            if (ay<0) alfa *= -1.0 ;
            alfa = atan2(ay, ax) * 180 / 3.1415 ;
            
            double gridAngleInSearchArea =  fmod (alfa - orientationBacteria , 360 ) ;       // using alfa for searching among grids, Using fmod, because the result should be in the range (-180, 180 ) degree
            if (gridAngleInSearchArea > 180)
            {
                gridAngleInSearchArea -= 360.0 ;
            }
            else if (gridAngleInSearchArea < -180)
            {
                gridAngleInSearchArea += 360.0 ;
            }
            
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
                        if ( slime[gridX][gridY] - liqBackground >= smallValue)
                        {
                            gridWithSlime[0] += 1 ;
                        }
                    }
                    else if (gridAngleInSearchArea > -54 && gridAngleInSearchArea <= -18)
                    {
                        
                        slimeRegions[1] += slime[gridX][gridY] ;
                        gridInRegions[1] += 1 ;
                        if ( slime[gridX][gridY] - liqBackground  >= smallValue)
                        {
                            gridWithSlime[1] += 1 ;
                        }
                    }
                    else if (gridAngleInSearchArea> -18 && gridAngleInSearchArea <= 18)
                    {
                        
                        slimeRegions[2] += slime[gridX][gridY] ;
                        gridInRegions[2] += 1 ;
                        if ( slime[gridX][gridY] - liqBackground  >= smallValue)
                        {
                            gridWithSlime[2] += 1 ;
                        }
                    }
                    
                    else if (gridAngleInSearchArea > 18 && gridAngleInSearchArea <= 54)
                    {
                        
                        slimeRegions[3] += slime[gridX][gridY] ;
                        gridInRegions[3] += 1 ;
                        if ( slime[gridX][gridY] - liqBackground  >= smallValue)
                        {
                            gridWithSlime[3] += 1 ;
                        }
                    }
                    
                    else if (gridAngleInSearchArea > 54 && gridAngleInSearchArea <= 90)
                    {
                        
                        slimeRegions[4] += slime[gridX][gridY] ;
                        gridInRegions[4] += 1 ;
                        if ( slime[gridX][gridY] - liqBackground >= smallValue)
                        {
                            gridWithSlime[4] += 1 ;
                        }
                    }
                    
                    totalSlimeInSearchArea += slime[gridX][gridY] ;
                    
                }
                
            }
            
        }
    }
    
    if (totalSlimeInSearchArea < smallValue)
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
            if ( (slimeRegions[i] > lowThreshold * maxSlimeRegion))
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
        alfa = ( (region -2 ) * 36.0 + orientationBacteria ) * 3.1415/180.0 ;
        // using alfa(radiant) for direction of the force
        // Fs = EpsilonS *(ft/(N-1))* (S/S0) ( toward slime)
        
        //cout<< orientationBacteria<<'\t'<< alfa<< '\t'<<(region -2 ) * 36.0<<endl;
        
        bacteria[i].fSTFx = slimeEffectiveness * tmpFmotor * cos(alfa) ;
        bacteria[i].fSTFy = slimeEffectiveness * tmpFmotor * sin(alfa) ;
        // (S/S0) is not included yet
    }
    return  sqrt(bacteria[i].fSTFx * bacteria[i].fSTFx + bacteria[i].fSTFy * bacteria[i].fSTFy) ;
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria:: Reverse_IndividualBacteriaa (int i)
{
   // bacteria[i].directionOfMotion = ! bacteria[i].directionOfMotion ;
    bacteria[i].turnStatus = true ;
    bacteria[i].turnTime = 0.0 ;
   
  //  bacteria[i].turnAngle =(2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) *  bacteria[i].maxTurnAngle ;    //uniform distribution
    if (bacteria[i].attachedToFungi == false || inLiquid == true)
    {
        //bacteria[i].turnAngle = gasdev(&idum) *  bacteria[i].maxTurnAngle ;                //guassian distribution
        bacteria[i].turnAngle = (2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) * bacteria[i].maxTurnAngle ;
        //test turnAngle
        //bacteria[i].turnAngle = 0.0 ;
        /*
        std::default_random_engine generator ;
        std::normal_distribution<double> distribution( 0.0, bacteria[i].turnSDV  ) ;
        bacteria[i].turnAngle = distribution(generator);
         */
        //cout<< "turning angle for "<<i<<"th bacteria is: "<<bacteria[i].turnAngle << endl ;
        
    }
    else
    {
        bacteria[i].turnAngle = 0.0 ;
        //cout<< i<<" is attached to fungi "<<endl;
        //test attachtoFungi=false
        //bacteria[i].turnAngle = (2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) *  bacteria[i].maxTurnAngle ;
        /*
        std::default_random_engine generator ;
        std::normal_distribution<double> distribution( 0.0, bacteria[i].turnSDV ) ;
        bacteria[i].turnAngle = distribution(generator);
         */
    }
    
    //bacteria[i].turnAngle =  bacteria[i].maxTurnAngle ;
    if (bacteria[i].sourceWithin == false)
    {
        //bacteria[i].sourceWithin = bacteria[i].SourceRegion() ;
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
    Cal_BacteriaOrientation (i) ;
    
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria:: Check_Perform_AllReversing_andWrapping()
{
    for (int i=0; i<nbacteria; i++)
    {
        if (chemotacticMechanism == classic)
        {
            if (bacteria[i].sourceWithin == false)
            {
                bacteria[i].sourceWithin = bacteria[i].SourceRegion() ;
                bacteria[i].timeToSource += dt ;
            }
            if (bacteria[i].SourceRegion() )
            {
                bacteria[i].timeInSource += dt ;
            }
            
            if (bacteria[i].reversalTime >= bacteria[i].reversalPeriod && bacteria[i].wrapMode == false )
            {
                bacteria[i].reversalTime  = 0.0 ;
                if ( rand() / (RAND_MAX + 1.0) < bacteria[i].wrapProbability )
                {
                    bacteria[i].turnStatus = false ;
                    bacteria[i].turnTime = 0.0 ;
                    
                    bacteria[i].wrapMode = true ;
                    bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(wrap_distribution,wrap_seed, lognormal_wrap_a, run_calibrated, 1.0/bacteria[i].wrapRate ) ;
                  //  bacteria[i].wrapAngle = (2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) * bacteria[i].maxWrapAngle ;
                    
                    std::default_random_engine generator ;
                    std::normal_distribution<double> distribution(bacteria[i].wrapMeanAngle , bacteria[i].wrapSDV ) ;
                    bacteria[i].wrapAngle = distribution(generator);
                    if ( rand() / (RAND_MAX + 1.0) < 0.5)
                    {
                        bacteria[i].wrapAngle *= -1.0 ;
                    }
                }
                else
                {
                    /*  // check to see whether these are neccessary or not
                    bacteria[i].wrapMode = false ;
                    bacteria[i].wrapTime = 0.0 ;
                    bacteria[i].wrapAngle = 0.0 ;
                     */
                    bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(run_distribution, run_seed, lognormal_run_a, run_calibrated, 1.0/reversalRate) ;
                    Reverse_IndividualBacteriaa(i) ;
                }
            }
            if (bacteria[i].wrapTime > bacteria[i].wrapDuration)
            {
                bacteria[i].wrapMode = false ;
                bacteria[i].wrapTime = 0.0 ;
                bacteria[i].wrapAngle = 0.0 ;
                
                bacteria[i].reversalTime  = 0.0 ;
                bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(run_distribution, run_seed, lognormal_run_a, run_calibrated, 1.0/reversalRate) ;
                Reverse_IndividualBacteriaa(i) ;
                
            }
            bacteria[i].reversalTime += dt ;
        }
        else        // chemotacticMechanism==true
        {
            if (bacteria[i].sourceWithin == false)
            {
                bacteria[i].sourceWithin = bacteria[i].SourceRegion() ;
                bacteria[i].timeToSource += dt ;
            }
            if (bacteria[i].SourceRegion() )
            {
                bacteria[i].timeInSource += dt ;
            }
            
            if (bacteria[i].wrapMode == false && bacteria[i].motilityMetabolism.switchMode == true )
            {
                bacteria[i].reversalTime  = 0.0 ;
                if ( rand() / (RAND_MAX + 1.0) < bacteria[i].wrapProbability )
                {
                    bacteria[i].turnStatus = false ;
                    bacteria[i].turnTime = 0.0 ;
                    
                    bacteria[i].wrapMode = true ;
                    bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(wrap_distribution,wrap_seed, lognormal_wrap_a, run_calibrated, 1.0/bacteria[i].wrapRate) ;
                  //  bacteria[i].wrapAngle = (2.0*(rand() / (RAND_MAX + 1.0))-1.0 ) * bacteria[i].maxWrapAngle ;
                    
                    std::default_random_engine generator ;
                    std::normal_distribution<double> distribution(bacteria[i].wrapMeanAngle , bacteria[i].wrapSDV ) ;
                    bacteria[i].wrapAngle = distribution(generator);
                    if ( rand() / (RAND_MAX + 1.0) < 0.5)
                    {
                        bacteria[i].wrapAngle *= -1.0 ;
                    }
                }
                else
                {
                    /*  // check to see whether these are neccessary or not
                    bacteria[i].wrapMode = false ;
                    bacteria[i].wrapTime = 0.0 ;
                    bacteria[i].wrapAngle = 0.0 ;
                     */
                    bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(run_distribution, run_seed, lognormal_run_a, run_calibrated, 1.0/reversalRate) ;
                    Reverse_IndividualBacteriaa(i) ;
                    bacteria[i].motilityMetabolism.switchMode = false ;
                }
            }
            else if ( bacteria[i].wrapMode == true && bacteria[i].motilityMetabolism.switchMode == true)
            {
                bacteria[i].wrapMode = false ;
                bacteria[i].wrapTime = 0.0 ;
                bacteria[i].wrapAngle = 0.0 ;
                
                bacteria[i].reversalTime  = 0.0 ;
                bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(run_distribution, run_seed, lognormal_run_a, run_calibrated, 1.0/reversalRate) ;
                Reverse_IndividualBacteriaa(i) ;
                bacteria[i].motilityMetabolism.switchMode = false ;
                
            }
            bacteria[i].reversalTime += dt ;
            
        }
    }
}


//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: Update_Bacteria_AllNodes ()
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


void TissueBacteria:: Initialize_BacteriaProteinLevel ()
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
void TissueBacteria:: Bacterial_ProteinExchange ()
{
    for (int i=0; i<nbacteria-1; i++)
    {
        for (int j=i+1; j<nbacteria; j++)
        {
            if (bacteria[i].connectionToOtherB[j] != 0.0)
            {
                double nbar = (bacteria[i].protein + bacteria[j].protein)/2 ;
                bacteria[i].protein += (nbar-bacteria[i].protein)* Rp * dt * ( bacteria[i].connectionToOtherB[j]/(2*nnode-1) ) ;
                bacteria[j].protein += (nbar-bacteria[j].protein)* Rp * dt * ( bacteria[j].connectionToOtherB[i]/(2*nnode-1) ) ;
            }
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------
void TissueBacteria:: Update_NodeLevel_ProteinConcentraion ()
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
void TissueBacteria:: BacterialVisualization_ParaView ()
{
    Update_NodeLevel_ProteinConcentraion() ;
    Bacteria_CreateDuplicate() ;
    Bacteria_MergeWithDuplicate() ;
    int index = index1 ;
    string vtkFileName = folderName + animationName + to_string(index)+ ".vtk" ;
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
    index1++ ;
}

//-----------------------------------------------------------------------------------------------------
/*
void TissueBacteria:: ParaView2 ()
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

*/

//-----------------------------------------------------------------------------------------------------
void TissueBacteria:: Bacteria_CreateDuplicate()
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

void TissueBacteria:: Bacteria_MergeWithDuplicate ()
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

void TissueBacteria:: Update_SlimeConcentration ()
{
    int m = 0 ;
    int n = 0 ;
    //Calculate the decay over time
    for (m=0; m<nx; m++)
    {
        for (n=0; n<ny; n++)
        {
            cout<< " WARNING: Decay is not possible when slime < liqBackground " ;
            slime[m][n] -= kd * dt * (slime[m][n]-liqBackground)  ;
        }
    }
    //calculate the secretion of slime based on where bacteria are located
    for (int i=0; i<nbacteria; i++)
    {
        for (int j=0; j<2*nnode-1 ; j++)
        {
            // we add domain to x and y in order to make sure m and n would not be negative integers
            m = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].x + domainx , domainx) / dx ) ) ) % nx  ;
            n = (static_cast<int> (round ( fmod (bacteria[i].allnodes[j].y + domainy , domainy) / dy ) ) ) % ny  ;
            slime[m][n] += sr * dt ;
            //   visit[m][n] += 1.0 ;      //undating visits in each time step. if you make it as uncomment, you should make it comment in the VisitsPerGrid function
            
        }
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: Handle_BacteriaTurnOrientation ()
{
    for (int i=0; i<nbacteria; i++)
    {
        if (bacteria[i].turnStatus)
        {
            double fTotal =sqrt(bacteria[i].nodes[0].fMotorx * bacteria[i].nodes[0].fMotorx +
                                bacteria[i].nodes[0].fMotory * bacteria[i].nodes[0].fMotory) ;
            double orientation = Cal_BacteriaOrientation(i) ;
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
        
        if (bacteria[i].wrapMode)
        {
            double fTotal =sqrt(bacteria[i].nodes[0].fMotorx * bacteria[i].nodes[0].fMotorx +
                                bacteria[i].nodes[0].fMotory * bacteria[i].nodes[0].fMotory) ;
            double orientation = Cal_BacteriaOrientation(i) ;
            bacteria[i].nodes[0].fMotorx = fTotal * cos(bacteria[i].wrapAngle + orientation ) ;
            bacteria[i].nodes[0].fMotory = fTotal * sin(bacteria[i].wrapAngle + orientation ) ;
            bacteria[i].wrapTime += dt ;
            
        }
        if (bacteria[i].turnStatus && bacteria[i].wrapMode )
        {
            cout<< i <<'\t'<<bacteria[i].turnTime<<'\t'<< bacteria[i].wrapTime<<endl ;
        }
    }
}
//-----------------------------------------------------------------------------------------------------

double TissueBacteria:: Cal_BacteriaOrientation (int i)

{
    
    double ax = 1.0 *(bacteria[i].nodes[1].x - bacteria[i].nodes[0].x) ;
    double ay = 1.0 *(bacteria[i].nodes[1].y - bacteria[i].nodes[0].y) ;
    double Cos = ax/sqrt(ax*ax+ay*ay) ;
    double orientationBacteria = acos(Cos)  ;        // orientation of the bacteria in radian
    if(ay<0.0) orientationBacteria *= -1.0 ;
    bacteria[i].orientation = orientationBacteria ;
    return orientationBacteria ;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: UpdateReversalFrequency ()
{
   
    for (int i=0; i< nbacteria ; i++)
    {
        bacteria[i].chemotaxisTime += 1.0/ updatingFrequency ;
        //Controlls how often the bacteria do chemotaxis( lazy bacteria)
        if (bacteria[i].chemotaxisTime < bacteria[i].chemotaxisPeriod )
        {
            continue ;
        }
        else
        {
            bacteria[i].chemotaxisTime = 0.0 ;
            double ds = Cal_ChemoGradient2(i) ;
            double tmpChemoStrength = chemoStrength_run ;
            if (bacteria[i].wrapMode)
            {
                tmpChemoStrength = chemoStrength_wrap ;
            }
            /*
             if (i=0)
             {
             cout<< ds<<endl ;
             }
             */
            Cal_BacteriaOrientation(i) ;
            if (ds >= 0.0)
            {
                double preferedAngle = atan2(bacteria[i].nodes[(nnode-1)/2].y - domainy/2.0, bacteria[i].nodes[(nnode-1)/2].x - domainx/2.0) ;
                /*
                 if (cos(bacteria[i].orientation - preferedAngle) < 0.0  )
                 {
                 //cout<< i<< '\t'<< bacteria[i].orientation<< '\t'<< preferedAngle<<'\t'<<cos(bacteria[i].orientation - preferedAngle)<<endl ;
                 }
                 */
                //optimistic
                //double tmpPeriod = 1.0 /( reversalRate * exp(tmpChemoStrength * (-1.0 *fmotor  /* *cos(bacteria[i].orientation - preferedAngle )  */ )) ) ;
                //bacteria[i].reversalPeriod = max(tmpPeriod, minimumRunTime ) ;
                //pessimistic
                //bacteria[i].reversalPeriod = 1.0/ reversalRate ;      //constant run durations
                if (bacteria[i].wrapMode == false)
                {
                    bacteria[i].reversalPeriod = max(bacteria[i].maxRunDuration, minimumRunTime) ;
                }
                else
                {
                    bacteria[i].wrapDuration = max(bacteria[i].maxRunDuration, minimumRunTime) ;
                }
                
            }
            else
            {
                //optimistic
                //bacteria[i].reversalPeriod = 1.0/ reversalRate ;
                
                //pessimistic
                //double tmpPeriod =1.0 /( reversalRate * exp(tmpChemoStrength * (-1.0 *fmotor * ds   /* *cos(bacteria[i].orientation - preferedAngle )  */ )) ) ;         // constant run duration
                double tmpPeriod =  bacteria[i].maxRunDuration /( exp(tmpChemoStrength * (-1.0 *fmotor * ds   /* *cos(bacteria[i].orientation - preferedAngle )  */ )) ) ;
                //test
                if (bacteria[i].wrapMode == false)
                {
                    bacteria[i].reversalPeriod = max(tmpPeriod, minimumRunTime ) ;
                }
                else
                {
                    bacteria[i].wrapDuration = max(tmpPeriod, minimumRunTime ) ;
                }
                
            }
        }
    }
    
}
//-----------------------------------------------------------------------------------------------------

double TissueBacteria:: Cal_ChemoGradient (int i)
{
    double c0 = 1.0 ;   //amplitude of chemical
   // double ds = bacteria[i].nodes[(nnode-1)/2].x - bacteria[i].oldLoc.at(0) ;
    double r2 = sqrt( pow( (bacteria[i].nodes[(nnode-1)/2].x - domainx/2.0 ) ,2 ) + pow( (bacteria[i].nodes[(nnode-1)/2].y - domainy/2.0 ) ,2 ) );
    double r1 = sqrt( pow( (bacteria[i].oldLoc.at(0) - domainx/2.0 ) ,2 ) + pow( (bacteria[i].oldLoc.at(1) - domainy/2.0 ) ,2 ) ) ;
    double ds = -1.0* c0 * (r2 - r1) ;
    //cout<< ds<<endl ;
    return ds;
}
//-----------------------------------------------------------------------------------------------------

double TissueBacteria:: Cal_ChemoGradient2 (int i)
{
    double c0 = 1.0 ;   //amplitude of chemical
    // double ds = bacteria[i].nodes[(nnode-1)/2].x - bacteria[i].oldLoc.at(0) ;
    int tmpXIndex = static_cast<int>( fmod( round(fmod( bacteria[i].nodes[(nnode-1)/2].x + 10.0 * domainx, domainx ) / tGrids.grid_dx ), tGrids.numberGridsX) ) ;
    int tmpYIndex = static_cast<int>( fmod( round(fmod( bacteria[i].nodes[(nnode-1)/2].y + 10.0 * domainy, domainy ) / tGrids.grid_dx ), tGrids.numberGridsX )) ;
    if (tmpXIndex < 0 || tmpXIndex > tGrids.numberGridsX - 1 || tmpYIndex < 0 || tmpYIndex > tGrids.numberGridsY - 1)
    {
        cout <<"Cal_ChemoGradient, index out of range "<<tmpXIndex<<'\t'<<tmpYIndex<<endl ;
    }
    double ds = gridInMain.at(tmpYIndex).at(tmpXIndex) - bacteria[i].oldChem ;
   // cout<<ds<<endl ;
    bacteria[i].oldChem = gridInMain.at(tmpYIndex).at(tmpXIndex) ;
    return ds;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: WriteTrajectoryFile ()
{
    ofstream trajectories (statsFolder +"trajectories.txt", ofstream::app) ;
    for (uint i = 0; i< nbacteria; i++)
    {
        trajectories << bacteria[i].nodes.at((nnode-1)/2).x <<'\t'<<bacteria[i].nodes.at((nnode-1)/2).y <<endl ;
    }
    trajectories<< endl ;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: WriteNumberReverse ()
{
    ofstream histogramReversal (statsFolder + "HistogramReversal.txt") ;
    for (int i = 0; i < nbacteria ; i++)
    {
        histogramReversal << i << '\t' << bacteria[i].numberReverse << '\t' << bacteria[i].timeToSource << '\t'<<bacteria[i].timeInSource<<'\t'<<bacteria[i].sourceWithin<<'\t'<< bacteria[i].maxRunDuration << endl ;
    }
}
//-----------------------------------------------------------------------------------------------------

//The entire hyphae width is filled with liquid. Liquid decreass with distance from hyphae
void TissueBacteria:: Find_FungalNetworkTrace (Fungi tmpFng)
{
    int m = 0 ;
    int n = 0 ;
    double xMin;
    double xMax ;
    double yMin ;
    double yMax ;
    int mMin = 0 ;
    int mMax= nx ;
    int nMin = 0 ;
    int nMax = ny ;
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
        mMin -= 4; mMax += 4; nMin -= 4; nMax += 4;
        
        tmpL = Dist2D(tmpFng.hyphaeSegments.at(i).x1, tmpFng.hyphaeSegments.at(i).y1, tmpFng.hyphaeSegments.at(i).x2, tmpFng.hyphaeSegments.at(i).y2) ;
        //going through nearby grids
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
                
                vec2x= tmpFng.hyphaeSegments.at(i).x1 - tmpFng.hyphaeSegments.at(i).x2 ;
                vec2y= tmpFng.hyphaeSegments.at(i).y1 - tmpFng.hyphaeSegments.at(i).y2 ;
                double tmpAngle1 = AngleOfTwoVectors(vec1x, vec1y, vec2x, vec2y) ;
                
                vec1x = tmpFng.hyphaeSegments.at(i).x2 - m * dx ;
                vec1y = tmpFng.hyphaeSegments.at(i).y2 - n * dy;
                vec2x *= -1.0 ;
                vec2y *= -1.0 ;
                
                double tmpAngle2 = AngleOfTwoVectors(vec1x, vec1y, vec2x, vec2y) ;
                //angles are used to define the rectangle as a hyphae segment
                if (tmpH < dx * tmpFng.hyphaeWidth && tmpAngle1 < pi/2.0 && tmpAngle2 < pi/2.0 )
                {
                    slime[m][n] = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH),slime[m][n]) ;
                    
                    /*
                    double decay = 2.0;
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
                        if (tmpAngle < pi/ 20.0 &&  tmpDis < tmplength   )
                        {
                            
                        //slime[m][n] += exp(-tmpDis/ decay) ;
                            slime [m][n] = max(s0 + s0* exp(-tmpH),slime[m][n]) ;
                            
                        }
                     }
                     */
                }
                //detect the grids close to the hyphae to make them a source of chemoattractant if we want secretion along the whole network
                if (tmpH < dx * tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq && tmpAngle1 < pi/2.0 && tmpAngle2 < pi/2.0 && sourceAlongHyphae== true )
                {
                    sourceChemo.at(0).push_back(m*dx) ;
                    sourceChemo.at(1).push_back(n*dy) ;
                    double tmpXtoX1 = Dist2D(tmpFng.hyphaeSegments.at(i).x1, tmpFng.hyphaeSegments.at(i).y1, m*dx, n*dy ) ;
                    
                    //linear secretion change along the hyphae segment. ( constant if p2==p1)
                    double pSource = tmpFng.hyphaeSegments.at(i).p1 + (tmpFng.hyphaeSegments.at(i).p2 - tmpFng.hyphaeSegments.at(i).p1)*(tmpXtoX1/tmpL) ;
                    sourceProduction.push_back(pSource) ;
                    
                }
    
                
            }
            
        }
    }
}

//-----------------------------------------------------------------------------------------------------

//Update the liquid concentration based on relative location with respect to fungi network
void TissueBacteria:: Find_FungalNetworkTrace2 (Fungi tmpFng)
{
    int m = 0 ;
    int n = 0 ;
    double xMin;
    double xMax ;
    double yMin ;
    double yMax ;
    int mMin = 0 ;
    int mMax= nx ;
    int nMin = 0 ;
    int nMax = ny ;
    double tmpH ;
    double tmpS ;
    double tmpL ;
    double vec1x , vec1y, vec2x , vec2y ;
    vector<vector<pair<bool, double>> > tmpHyphaeEdge ;
    //false when the grid is on top of at least one hyphae segment
    tmpHyphaeEdge.resize(nx, vector<pair<bool, double> > (ny, make_pair(true, liqBackground) ) ) ;
    
    
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
        int tmpNghbrhood = static_cast<int>( 1.5 * tmpFng.hyphaeWidth / dx) ;
        mMin -= tmpNghbrhood ; mMax += tmpNghbrhood ; nMin -= tmpNghbrhood ; nMax += tmpNghbrhood;
        
        tmpL = Dist2D(tmpFng.hyphaeSegments.at(i).x1, tmpFng.hyphaeSegments.at(i).y1, tmpFng.hyphaeSegments.at(i).x2, tmpFng.hyphaeSegments.at(i).y2) ;
        //going through nearby grids
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
                
                vec2x= tmpFng.hyphaeSegments.at(i).x1 - tmpFng.hyphaeSegments.at(i).x2 ;
                vec2y= tmpFng.hyphaeSegments.at(i).y1 - tmpFng.hyphaeSegments.at(i).y2 ;
                double tmpAngle1 = AngleOfTwoVectors(vec1x, vec1y, vec2x, vec2y) ;
                
                vec1x = tmpFng.hyphaeSegments.at(i).x2 - m * dx ;
                vec1y = tmpFng.hyphaeSegments.at(i).y2 - n * dy;
                vec2x *= -1.0 ;
                vec2y *= -1.0 ;
                
                double tmpAngle2 = AngleOfTwoVectors(vec1x, vec1y, vec2x, vec2y) ;
                //angles are used to define the rectangle as a hyphae segment
                if (tmpH < dx * tmpFng.hyphaeWidth && tmpAngle1 < pi/2.0 && tmpAngle2 < pi/2.0 )
                {
                    //the grid is inside of the hyphae segment
                    if (tmpH < dx * tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq)
                    {
                        tmpHyphaeEdge[m][n].first = false ;
                    }
                    else
                    {
                        /* concentration exponetially decays with distance from hyphae ( but there is a cut off
                        tmpHyphaeEdge[m][n].second = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                        slime[m][n] = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                         */
                        tmpHyphaeEdge[m][n].second = liqLayer ;
                        slime[m][n] = liqLayer ;
                    }
                    
                }
                //detect the grids close to the hyphae to make them a source of chemoattractant if we want secretion along the whole network
                if (tmpH < dx * tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq && tmpAngle1 < pi/2.0 && tmpAngle2 < pi/2.0 && sourceAlongHyphae== true )
                {
                    sourceChemo.at(0).push_back(m*dx) ;
                    sourceChemo.at(1).push_back(n*dy) ;
                    double tmpXtoX1 = Dist2D(tmpFng.hyphaeSegments.at(i).x1, tmpFng.hyphaeSegments.at(i).y1, m*dx, n*dy ) ;
                    //linear secretion change along the hyphae segment. ( constant if p2==p1)
                    double pSource = tmpFng.hyphaeSegments.at(i).p1 + (tmpFng.hyphaeSegments.at(i).p2 - tmpFng.hyphaeSegments.at(i).p1)*(tmpXtoX1/tmpL) ;
                    sourceProduction.push_back(pSource) ;
                    
                }
    
                
            }
            
        }
    }
    //Add semi-circle at the endpoint of each hyphae segment
    for (uint i=0; i<tmpFng.hyphaeSegments.size(); i++ )
    {
        xMin = tmpFng.hyphaeSegments.at(i).x2 ;
        yMin = tmpFng.hyphaeSegments.at(i).y2 ;
        
        mMin = (static_cast<int> (floor ( fmod (xMin - tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq  + domainx , domainx) / dx ) ) ) % nx  ;
        mMax = (static_cast<int> (ceil ( fmod (xMin + tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainx , domainx) / dx ) ) ) % nx  ;
        nMin = (static_cast<int> (floor ( fmod (yMin - tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainy , domainy) / dy ) ) ) % ny  ;
        nMax = (static_cast<int> (ceil ( fmod (yMin + tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainy , domainy) / dy ) ) ) % ny  ;
        for (int m = mMin ; m <= mMax ; m++)
        {
            for (int n = nMin; n <= nMax; n++)
            {
                tmpH = Dist2D(xMin, yMin, m * dx, n * dy) ;
                if ( tmpH < 0.5 * tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq  )
                {
                    tmpHyphaeEdge[m][n].first = false ;
                }
                else if ( tmpH < tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq )
                {
                    /* concentration exponetially decays with distance from hyphae ( but there is a cut off
                    tmpHyphaeEdge[m][n].second = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                    slime [m][n] = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                     */
                    tmpHyphaeEdge[m][n].second = liqLayer ;
                    slime[m][n] = liqLayer ;
                }
            }
        }
        //Add the semi-circle for the very first segments
        if (tmpFng.hyphaeSegments.at(i).from_who == -1)
            {
                xMin = tmpFng.hyphaeSegments.at(i).x1 ;
                yMin = tmpFng.hyphaeSegments.at(i).y1 ;
        
                //xMin = tmpFng.hyphaeSegments.at(i).x2 ;
                //yMin = tmpFng.hyphaeSegments.at(i).y2 ;
        
                mMin = (static_cast<int> (floor ( fmod (xMin - tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq  + domainx , domainx) / dx ) ) ) % nx  ;
                mMax = (static_cast<int> (ceil ( fmod (xMin + tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainx , domainx) / dx ) ) ) % nx  ;
                nMin = (static_cast<int> (floor ( fmod (yMin - tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainy , domainy) / dy ) ) ) % ny  ;
                nMax = (static_cast<int> (ceil ( fmod (yMin + tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq + domainy , domainy) / dy ) ) ) % ny  ;
                for (int m = mMin ; m <= mMax ; m++)
                {
                for (int n = nMin; n <= nMax; n++)
                    {
                        tmpH = Dist2D(xMin, yMin, m * dx, n * dy) ;
                        if ( tmpH < 0.5 * tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq  )
                        {
                            tmpHyphaeEdge[m][n].first = false ;
                        }
                        else if ( tmpH < tmpFng.hyphaeWidth * tmpFng.hyphaeOverLiq )
                        {
                            /* concentration exponetially decays with distance from hyphae ( but there is a cut off
                            tmpHyphaeEdge[m][n].second = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                            slime [m][n] = max(5.0*liqBackground + 5.0 *liqBackground * exp(-tmpH/tmpFng.hyphaeWidth ),slime[m][n]) ;
                             */
                            tmpHyphaeEdge[m][n].second = liqLayer ;
                            slime[m][n] = liqLayer ;
                        }
                    }
                 }
            }
    }
    //Update values for the grids based on where with respect to the fungi network
    for (int i=0 ; i< tmpHyphaeEdge.size(); i++)
    {
        for (int j=0 ; j< tmpHyphaeEdge.at(i).size() ; j++)
        {
            if (tmpHyphaeEdge.at(i).at(j).first == true)
            {
                slime[i][j] = tmpHyphaeEdge.at(i).at(j).second ;
            }
            else
            {
                slime[i][j] = liqHyphae ;
            }
        }
    }
            
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Update_MotilityMetabolism_Only(double tmpDt)
{
    Update_MM_Legand() ;
    for (int i=0 ; i< nbacteria ; i++)
    {
        bacteria[i].motilityMetabolism.Cal_SwitchProbability(tmpDt) ;
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Update_MotilityMetabolism_Only2(double tmpDt)
{
    Update_MM_Legand() ;
    for (int j=0; j < 10; j++)
    {
        for (int i=0 ; i< nbacteria ; i++)
        {
            bacteria[i].motilityMetabolism.Cal_SwitchProbability(tmpDt) ;
        }
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Update_MotilityMetabolism(double tmpDt)
{
    Update_MM_Legand() ;
    
    for (int i=0 ; i< nbacteria ; i++)
    {
        bacteria[i].motilityMetabolism.Cal_SwitchProbability(tmpDt) ;
        if ( bacteria[i].motilityMetabolism.switchMode == false &&
            rand() / (RAND_MAX + 1.0) < bacteria[i].motilityMetabolism.switchProbability)
        {
            bacteria[i].motilityMetabolism.switchMode = true ;
        }
    }
    WriteSwitchProbabilities() ;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: WriteSwitchProbabilities()
{
    ofstream strSwitchP (statsFolder + "SwitchProbablities.txt", ofstream::app ) ;
    for (int i = 0; i < 1 ; i++)
    {
       // strSwitchP << bacteria[i].motilityMetabolism.switchProbability << '\t'  ;
        strSwitchP  << bacteria[i].motilityMetabolism.receptorActivity << '\t'
                    << bacteria[i].motilityMetabolism.legand << '\t'
                    << bacteria[i].motilityMetabolism.methylation << '\t'
                    << bacteria[i].motilityMetabolism.LegandEnergy << '\t'
                    << bacteria[i].motilityMetabolism.MethylEnergy << '\t'
                    << bacteria[i].motilityMetabolism.switchProbability  ;
                    
    }
    strSwitchP<< endl ;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria:: WriteSwitchProbabilitiesByBacteria()
{   for (int i=0; i<nbacteria; i++)
    {
    ofstream strSwitchP2;
    strSwitchP2.open(statsFolder + "WriteSwitchProbabilities_Bacteria" + to_string(i) + ".txt", ios::app);
    {
       // strSwitchP << bacteria[i].motilityMetabolism.switchProbability << '\t'  ;
        strSwitchP2 << setw(10) << bacteria[i].maxRunDuration << '\t' 
                    << setw(10) << bacteria[i].motilityMetabolism.receptorActivity << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.legand << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.methylation << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.changeRate << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.LegandEnergy << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.MethylEnergy << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.switchProbability  << '\t'
                    << setw(10) << bacteria[i].motilityMetabolism.switchMode  << '\t';
                    
    }
    strSwitchP2<< endl ;
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Update_MM_Legand()
{
    for (int i = 0; i< nbacteria; i++)
    {
        int tmpXIndex = static_cast<int>( fmod( round(fmod( bacteria[i].nodes[(nnode-1)/2].x + 10.0 * domainx, domainx ) / tGrids.grid_dx ), tGrids.numberGridsX) ) ;
        int tmpYIndex = static_cast<int>( fmod( round(fmod( bacteria[i].nodes[(nnode-1)/2].y + 10.0 * domainy, domainy ) / tGrids.grid_dy ), tGrids.numberGridsY )) ;
        if (tmpXIndex < 0 || tmpXIndex > tGrids.numberGridsX - 1 || tmpYIndex < 0 || tmpYIndex > tGrids.numberGridsY - 1)
        {
            cout <<"Cal_ChemoGradient, index out of range "<<tmpXIndex<<'\t'<<tmpYIndex<<endl<<flush ;
        }
        double s = gridInMain.at(tmpYIndex).at(tmpXIndex) ;
        bacteria[i].motilityMetabolism.legand = 1.0 * s ;
    }
    
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Pass_PointSources_To_Bacteria(vector<vector<double> > sourceP)
{
    for (int i=0; i< nbacteria; i++)
    {
        bacteria[i].chemoPointSources = sourceP ;
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::Update_BacteriaMaxDuration()
{
    if (run_calibrated)
    {
        for (int i=0; i<nbacteria; i++)
        {
            bacteria[i].maxRunDuration = bacteria[i].LogNormalMaxRunDuration(run_distribution, run_seed, lognormal_run_a, run_calibrated, 1.0/reversalRate ) ;
            
        }
    }
    else
    {
        for (int i=0; i<nbacteria; i++)
        {
            bacteria[i].maxRunDuration = 1.0 / reversalRate ;
        }
    }
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::WriteBacteria_AllStats()
{
    int index = index1 ;
    string statFileName = statsFolder + animationName + to_string(index)+ ".txt" ;
    ofstream stat;
    stat.open(statFileName.c_str());
    for (uint i = 0 ; i< nbacteria; i++)
    {
        stat<< bacteria[i].nodes[(nnode-1)/2].x <<'\t'<< bacteria[i].nodes[(nnode-1)/2].y<<'\t'
        <<bacteria[i].locVelocity << '\t' << bacteria[i].locFriction << '\t' << Cal_BacteriaOrientation(i) << '\t' << bacteria[i].oldChem  ;
        
        stat<<endl ;
    }
    stat.close() ;
}
//-----------------------------------------------------------------------------------------------------

void TissueBacteria::LogLinear_RNG()
{
    std::lognormal_distribution<double> distribution(lognormal_run_m , lognormal_run_s ) ;
    std::lognormal_distribution<double> distribution2(lognormal_wrap_m , lognormal_wrap_s ) ;
    
    run_distribution = distribution ;
    wrap_distribution = distribution2 ;
}
//-----------------------------------------------------------------------------------------------------

vector<vector<double> > TissueBacteria::Find_secretion_Coord_Rate(Fungi tmpFng)
{
    vector<vector<double> > pointSources ;
    double equiv_prodRate = 0.0 ;
    if (sourceAlongHyphae)
    {
        pointSources = sourceChemo ;
        equiv_prodRate = tmpFng.production * tmpFng.tips.at(0).size() / sourceChemo.at(0).size() ;
        transform(sourceProduction.begin(), sourceProduction.end(), sourceProduction.begin(), productNum(equiv_prodRate) ) ;
        
    }
    else
    {
        // production at the tips
        pointSources = tmpFng.tips ;
        sourceChemo = tmpFng.tips ;
        equiv_prodRate = tmpFng.production ;
        sourceProduction.resize(pointSources.at(0).size()) ;
        fill(sourceProduction.begin(), sourceProduction.end(), equiv_prodRate) ;
    }
    return pointSources ;
}
//-----------------------------------------------------------------------------------------------------

 
