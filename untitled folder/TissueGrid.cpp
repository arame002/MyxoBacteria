
#include "TissueGrid.hpp"

TissueGrid::TissueGrid ()
{
    dx = static_cast<double> ( ( xDomainMax - xDomainMin ) / numberGridsX ) ;
    dy = static_cast<double> ( ( yDomainMax - yDomainMin ) / numberGridsY ) ;
    dt *= 1.0 /Diffusion ;
}

void TissueGrid::DomainBoundaries(double xMin, double xMax, double yMin, double yMax)
{
    xDomainMin = xMin ;
    xDomainMax = xMax ;
    yDomainMin = yMin ;
    yDomainMax = yMax ;
}



void TissueGrid::InitializeAllGrids()
{
    vector<vector <Grid> > tmp (numberGridsY,vector<Grid> ( numberGridsX) ) ;
    grids = tmp ;
    
    return ;
}

void TissueGrid::ClearChanges()
{
    for (int i =0; i < numberGridsY; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            
            grids.at(i).at(j).change = 0.0 ;
            
        }
    }
    return ;
}

void TissueGrid::DiffusionChanges()
{
    double tmpChange = 0.0 ;
    
    for (int i =0; i < numberGridsY; i++)
    {
        for (int j = 0; j< numberGridsX -1 ; j++)
        {
            tmpChange = Diffusion * ( grids.at(i).at(j+1).value - grids.at(i).at(j).value) * dt / (dx * dx) ;
            grids.at(i).at(j).change += tmpChange ;
            grids.at(i).at(j+1).change += -tmpChange ;

        }
    }
    for (int i =0; i < numberGridsY -1; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            tmpChange = Diffusion * (grids.at(i+1).at(j).value - grids.at(i).at(j).value) * dt / (dy * dy) ;
            grids.at(i).at(j).change += tmpChange ;
            grids.at(i+1).at(j).change += -tmpChange ;
            
        }
    }
    return ;
}

void TissueGrid::FindProductionPoints()
{
    int tmpIndexX ;
    int tmpIndexY ;
    for (unsigned int i = 0; i < xSources.size(); i++)
    {
        tmpIndexX = static_cast<int>(round ( ( xSources.at(i) - xDomainMin ) / dx ) ) ;
        tmpIndexY = static_cast<int>(round ( ( ySources.at(i) - yDomainMin ) / dy ) ) ;
        grids.at(tmpIndexY).at(tmpIndexX).productionRate = pro ;
        indexSourceX.push_back(tmpIndexX) ;
        indexSourceY.push_back(tmpIndexY) ;
        
    }
    return ;
}

void TissueGrid::ProductionChanges()
{
    double tmpChange ;
    int tmpIdX ;
    int tmpIdY ;
    for (unsigned int i =0; i < indexSourceX.size() ; i++)
    {
        tmpIdX = indexSourceX.at(i) ;
        tmpIdY = indexSourceY.at(i) ;
        
        tmpChange = ( grids.at(tmpIdY).at(tmpIdX).productionRate) * dt ;
        grids.at(tmpIdY).at(tmpIdX).change += tmpChange ;
        
    }
    return ;
}

void TissueGrid::DegredationChanges()
{
    double tmpChange ;
    for (int i =0; i < numberGridsY ; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            tmpChange = - deg * (grids.at(i).at(j).value ) * dt ;
            grids.at(i).at(j).change += tmpChange ;
            
        }
    }
}

void TissueGrid::UpdateChanges()
{
    for (int i =0; i < numberGridsY ; i++)
    {
        for (int j = 0; j< numberGridsX ; j++)
        {
            grids.at(i).at(j).UpdateValue() ;
        }
    }
    return ;
}

void TissueGrid::EulerMethod()
{
    int l = 0 ;
    bool status = false ;
    while (status == false)
    {
        status = true ;
        double smallValue = 0.0001 ;
        ClearChanges() ;
        DiffusionChanges() ;
        DegredationChanges() ;
        ProductionChanges() ;
        for (int i=0; i< numberGridsY ; i++)
        {
            for (int j = 0; j< numberGridsX ; j++)
            {
                if (grids.at(i).at(j).change/( (grids.at(i).at(j).value + smallValue) ) > smallValue * dt  )
                {
                    status = false ;
                    break ;
                    
                }
            }
        }
        UpdateChanges() ;
        if (l%100 == 0)
        {
            ParaViewGrids(l/100) ;
        }
        l++ ;
        
    }
    cout<<l<<endl ;
    return ;
}

void TissueGrid::ParaViewGrids(int index)
{
    string vtkFileName2 = "GridChem"+ to_string(index)+ ".vtk" ;
    ofstream SignalOut;
    SignalOut.open(vtkFileName2.c_str());
    SignalOut << "# vtk DataFile Version 2.0" << endl;
    SignalOut << "Result for paraview 2d code" << endl;
    SignalOut << "ASCII" << endl;
    SignalOut << "DATASET RECTILINEAR_GRID" << endl;
    SignalOut << "DIMENSIONS" << " " << numberGridsX  << " " << " " << numberGridsY << " " << 1  << endl;
    
    SignalOut << "X_COORDINATES " << numberGridsX << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int i = 0; i < numberGridsX ; i++) {
        SignalOut << i * dx << endl;
    }
    
    SignalOut << "Y_COORDINATES " << numberGridsY << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int j = 0; j < numberGridsY; j++) {
        SignalOut << j * dy << endl;
    }
    
    SignalOut << "Z_COORDINATES " << 1 << " float" << endl;
    //write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
    for (int k = 0; k < 1 ; k++) {
        SignalOut << 0 << endl;
    }
    
    SignalOut << "POINT_DATA " << (numberGridsX )*( numberGridsY ) << endl;
    SignalOut << "SCALARS DPP float 1" << endl;
    SignalOut << "LOOKUP_TABLE default" << endl;
    
    for (int k = 0; k < 1 ; k++) {
        for (int j = 0; j < numberGridsX; j++) {
            for (int i = 0; i < numberGridsY ; i++) {
                SignalOut << grids.at(i).at(j).value << endl;
            }
        }
    }
    

}
