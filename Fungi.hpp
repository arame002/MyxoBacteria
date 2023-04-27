//
//  Fungi.hpp
//  Myxobacteria
//
//  Created by Alireza Ramezani on 11/19/20.
//  Copyright © 2020 Alireza Ramezani. All rights reserved.
//

#ifndef Fungi_hpp
#define Fungi_hpp


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#include "growthFunctions.h"

class Fungi {
public:
    
    vector<HyphaeSegment> hyphaeSegments ;
    vector<vector<double> > tips ;
    vector<int > tipsID ;
    vector<int> centerConnectedSegments ;
    vector<double> Hypahe_X1Values;
    vector<double> Hypahe_X2Values;
    vector<double> Hypahe_Y1Values;
    vector<double> Hypahe_Y2Values;
    double extendingAngle = 13.0 * (constants::pi / 180.0) ;
    double hyphaeWidth = 6.0 ;
    double production = pow(10, 2) ;
    double proDecayFactor = 0.5 ;
    double initX = 540.0 ;
    double initY = 500 ;
    int init_Count = 4 ;
    
    int machineID = 1 ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    
    Fungi () ;
    void Find_Hyphae_Tips() ;
    void Find_Hyphae_Tips2() ;
    void WriteSourceLoc(vector<vector<double> >) ;
    void FindProductionNetwork ();
    void FindProductionNetwork2();
    vector<int> FindNeighborList (int i) ;
    void Rec_CalNeighborProduction (int id, vector<int> ngbrList, int nbrElement, bool connection) ;
    void Rec_CalNeighborProduction2(int id, vector<int> ngbrList, int nbrElement ) ;
    vector<int> FindCenterPointConnection () ;
    void UpdateFungiFolderNames (int ) ;
    void UpdateFungi_FromConfigFile () ;
    
};


#endif /* Fungi_hpp */

