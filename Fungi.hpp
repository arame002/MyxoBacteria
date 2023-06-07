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
    //---------------------------- Parameters and sub-classes ------------------------------
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
    double hyphaeOverLiq = 0.5 ;
    double production = pow(10, 2) ;
    double proDecayFactor = 0.5 ;
    double initX = 500.0 ;
    double initY = 500 ;
    int init_Count = 4 ;
    bool loading_Network = true ;
    bool branchIsTip = false ;
    
    int machineID = 1 ;
    string folderName = "./animation/machine" + to_string(machineID) + "/" ;
    string statsFolder = "./dataStats/machine" + to_string(machineID) + "/" ;
    
    //---------------------------- Functions ---------------------------------------
    Fungi () ;
    void Find_Hyphae_Tips() ;               // Store the branch points and the tips in the "tips vector"
    void Find_Hyphae_Tips2() ;              //Only the tips of fungi network are stored in the "tips vector"
    void WriteSourceLoc(vector<vector<double> >) ;
    //Need calibration if you want to use it. (proDecayFactor, HyphaeSegment.p1 & p2, and etc. )
    //Scale production rates based on the network-distance from:
    void FindProductionNetwork ();          // the tips
    void FindProductionNetwork2();          //the Center of the network
    vector<int> FindNeighborList (int i) ;
    void Rec_CalNeighborProduction (int id, vector<int> ngbrList, int nbrElement, bool connection) ;
    void Rec_CalNeighborProduction2(int id, vector<int> ngbrList, int nbrElement ) ;
    vector<int> FindCenterPointConnection () ;
    void UpdateFungiFolderNames (int ) ;
    void UpdateFungi_FromConfigFile () ;
    
};


#endif /* Fungi_hpp */

