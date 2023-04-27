//
//  Fungi.cpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/19/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#include "Fungi.hpp"
//end points and branches are tips
Fungi::Fungi ()
{
    //UpdateFungi_FromConfigFile() ;
}


void Fungi::Find_Hyphae_Tips()
{
    tips.resize(2) ;
    for (uint i=0 ; i<hyphaeSegments.size() ; i++)
    {
        
        if (hyphaeSegments[i].can_extend == true)
        {
            /*
            for (int j=0; j<4; j++)
            {
            
                tips[0].push_back(hyphaeSegments[i].x2- j*(hyphaeSegments[i].x2-hyphaeSegments[i].x1)/3.0) ;
                tips[1].push_back(hyphaeSegments[i].y2- j*( hyphaeSegments[i].y2-hyphaeSegments[i].y1)/3.0) ;
            }
             */
             // tips only
            tips[0].push_back(hyphaeSegments[i].x2) ;
            tips[1].push_back(hyphaeSegments[i].y2) ;
            tipsID.push_back(i);
             
        }
         
        
    }
    for (uint i=0 ; i<hyphaeSegments.size() ; i++)
    {
        if (hyphaeSegments[i].can_branch == false)
        {
            tips[0].push_back(hyphaeSegments[i].x2) ;
            tips[1].push_back(hyphaeSegments[i].y2) ;
            tipsID.push_back(i);
        }
    }
}
//---------------------------------------------------------------------------------------------
//Only end points are tips
void Fungi::Find_Hyphae_Tips2()
{
    tips.resize(2) ;
    for (uint i=0 ; i<hyphaeSegments.size() ; i++)
    {
        
        //if (hyphaeSegments[i].can_extend == true && hyphaeSegments[i].can_branch == true )
        if (hyphaeSegments[i].can_extend == true )
        {
            //for motion inLiquid tip is beginning of element, so we can have the source in the center of domain
            tips[0].push_back(hyphaeSegments[i].x2) ;
            tips[1].push_back(hyphaeSegments[i].y2) ;
            //tips[0].push_back(hyphaeSegments[i].x1) ;
            //tips[1].push_back(hyphaeSegments[i].y1) ;
            tipsID.push_back(i);
        }
    }
}
//---------------------------------------------------------------------------------------------

void Fungi:: WriteSourceLoc(vector<vector<double> > pointSource)
{
    ofstream sourceLoc (statsFolder + "sourceLocation.txt") ;
    for (uint i = 0; i< pointSource.at(0).size() ; i++)
    {
        sourceLoc << pointSource.at(0).at(i) <<'\t'<< pointSource.at(1).at(i) <<endl ;
    }
    sourceLoc<< endl ;
}
//---------------------------------------------------------------------------------------------

void Fungi::FindProductionNetwork()
{
    vector<int> ngbrList ;
    for (int i = 0; i< tipsID.size() ; i++)
    {
        bool tmpConnection = false ;
        int id = tipsID.at(i) ;
        hyphaeSegments.at(id).p2 += production ;
        hyphaeSegments.at(id).p1 += production * proDecayFactor ;
        hyphaeSegments.at(id).tmpP1 = production * proDecayFactor ;
        hyphaeSegments.at(id).tmpP2 = production ;
        
        ngbrList = FindNeighborList(id) ;
        for (int nbrElement = 0; nbrElement< ngbrList.size(); nbrElement++)
        {
            Rec_CalNeighborProduction(id, ngbrList, nbrElement, tmpConnection) ;
        }
        
    }
    
}
//---------------------------------------------------------------------------------------------

void Fungi::FindProductionNetwork2()
{
    vector<int> ngbrList ;
    for (int i = 0; i< hyphaeSegments.size() ; i++)
    {
        if (hyphaeSegments.at(i).from_who == -1 )
        {
            hyphaeSegments.at(i).p1 = production ;
            hyphaeSegments.at(i).p2 = production * proDecayFactor ;
            hyphaeSegments.at(i).tmpP2 = production * proDecayFactor ;
            hyphaeSegments.at(i).tmpP1 = production ;
        
            ngbrList = FindNeighborList(i) ;
            for (int nbrElement = 0; nbrElement< ngbrList.size(); nbrElement++)
            {
                Rec_CalNeighborProduction2(i, ngbrList, nbrElement) ;
            }
        }
        
    }
    
}
//---------------------------------------------------------------------------------------------

vector<int> Fungi::FindNeighborList(int i)
{
    vector<int> ngbrList ;
    
    if( hyphaeSegments.at(i).from_who !=  -1 )
        ngbrList.push_back(hyphaeSegments.at(i).from_who) ;
    if( hyphaeSegments.at(i).can_branch != true)
        ngbrList.push_back(hyphaeSegments.at(i).branch_to) ;
    if( hyphaeSegments.at(i).can_extend != true)
        ngbrList.push_back(hyphaeSegments.at(i).extend_to) ;
    
    return ngbrList ;
}
//---------------------------------------------------------------------------------------------

void Fungi::Rec_CalNeighborProduction(int id, vector<int> ngbrList, int nbrElement, bool connection )
{
        int neighbor = ngbrList.at(nbrElement) ;
        if (neighbor < id)
        {
        hyphaeSegments.at(neighbor).p2 += hyphaeSegments.at(id).tmpP1 ;
        hyphaeSegments.at(neighbor).p1 += hyphaeSegments.at(id).tmpP1 * proDecayFactor ;
        hyphaeSegments.at(neighbor).tmpP1 = hyphaeSegments.at(id).tmpP1 * proDecayFactor ;
        hyphaeSegments.at(neighbor).tmpP2 = hyphaeSegments.at(id).tmpP1 ;
            
        }
        else
        {
            hyphaeSegments.at(neighbor).p1 += hyphaeSegments.at(id).tmpP2 ;
            hyphaeSegments.at(neighbor).p2 += hyphaeSegments.at(id).tmpP2 * proDecayFactor ;
            hyphaeSegments.at(neighbor).tmpP2 = hyphaeSegments.at(id).tmpP2 * proDecayFactor ;
            hyphaeSegments.at(neighbor).tmpP1 = hyphaeSegments.at(id).tmpP2 ;
        }
        vector<int> tmpList = FindNeighborList(neighbor) ;
        
        //remove id from tmpList
        for (int k =0 ; k<tmpList.size(); k++)
        {
            if (tmpList.at(k) == id )
            {
                tmpList.erase(tmpList.begin()+k ) ;
            }
        }
        for (int k =0 ; k<tmpList.size(); k++)
        {
            Rec_CalNeighborProduction(neighbor, tmpList, k, connection) ;
        }
    
    if (hyphaeSegments.at(neighbor).from_who == -1 &&  connection == false)
    {
        connection = true ;
        for (int i=0; i< centerConnectedSegments.size() ; i++)
        {
            int centerNeighbor = centerConnectedSegments.at(i) ;
            
            if (centerConnectedSegments.at(i) != neighbor )
            {
                hyphaeSegments.at(centerNeighbor ).p1 += hyphaeSegments.at(neighbor).tmpP1 ;
                hyphaeSegments.at(centerNeighbor ).p2 += hyphaeSegments.at(neighbor).tmpP1 * proDecayFactor ;
                hyphaeSegments.at(centerNeighbor ).tmpP2 = hyphaeSegments.at(neighbor).tmpP1 * proDecayFactor ;
                hyphaeSegments.at(centerNeighbor ).tmpP1 = hyphaeSegments.at(neighbor).tmpP1 ;
                
                
                vector<int> tmpListCenter = FindNeighborList(centerNeighbor) ;
                
                //remove id from tmpList
                for (int k =0 ; k<tmpListCenter.size(); k++)
                {
                    if (tmpListCenter.at(k) == neighbor )
                    {
                        tmpListCenter.erase(tmpListCenter.begin()+k ) ;
                    }
                }
                for (int k =0 ; k<tmpListCenter.size(); k++)
                {
                    Rec_CalNeighborProduction(centerNeighbor, tmpListCenter, k, connection) ;
                }
                
                
            }
        }
    }
    
    
    
   
}
//---------------------------------------------------------------------------------------------

void Fungi::Rec_CalNeighborProduction2(int id, vector<int> ngbrList, int nbrElement )
{
    
        int neighbor = ngbrList.at(nbrElement) ;
        
        hyphaeSegments.at(neighbor).p1 = hyphaeSegments.at(id).tmpP2 ;
        hyphaeSegments.at(neighbor).p2 = hyphaeSegments.at(id).tmpP2 * proDecayFactor ;
        hyphaeSegments.at(neighbor).tmpP2 = hyphaeSegments.at(id).tmpP2 * proDecayFactor ;
        hyphaeSegments.at(neighbor).tmpP1 = hyphaeSegments.at(id).tmpP2 ;
    
        vector<int> tmpList = FindNeighborList(neighbor) ;
        
        //remove id from tmpList
        for (int k =0 ; k<tmpList.size(); k++)
        {
            if (tmpList.at(k) == id )
            {
                tmpList.erase(tmpList.begin()+k ) ;
            }
        }
        for (int k =0 ; k<tmpList.size(); k++)
        {
            Rec_CalNeighborProduction2(neighbor, tmpList, k) ;
        }
}

//---------------------------------------------------------------------------------------------
vector<int> Fungi::FindCenterPointConnection()
{
    vector<int> cntNeighbor ;
    for (int i=0; i< hyphaeSegments.size(); i++)
    {
        if( hyphaeSegments.at(i).from_who ==  -1 )
        {
            cntNeighbor.push_back(i) ;
        }
    }
    centerConnectedSegments = cntNeighbor ;
    return cntNeighbor ;
}

//---------------------------------------------------------------------------------------------
void Fungi::UpdateFungiFolderNames(int id)
{
    machineID = id ;
    folderName = "./animation/machine" + to_string(machineID) + "/" ;
    statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    UpdateFungi_FromConfigFile() ;
    
}

void Fungi::UpdateFungi_FromConfigFile()
{
    folderName = globalConfigVars.getConfigValue("AnimationFolder").toString() ;
    statsFolder = globalConfigVars.getConfigValue("StatFolderName").toString() ;
    
    init_Count = globalConfigVars.getConfigValue("hyphae_initCount").toDouble() ;
    initX = globalConfigVars.getConfigValue("hyphae_InitPosX").toDouble() ;
    initY = globalConfigVars.getConfigValue("hyphae_InitPosY").toDouble() ;
    production = globalConfigVars.getConfigValue("hyphae_production").toDouble() ;
    proDecayFactor = globalConfigVars.getConfigValue("hyphae_proDecayFactor").toDouble() ;
    hyphaeWidth = globalConfigVars.getConfigValue("hyphae_Width").toDouble() ;
    
}
