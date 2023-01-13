
#include "driver.h"

using namespace std;
using constants::pi;


// ---------------------------------------------------
//   Main Driver Function
// ---------------------------------------------------
Fungi driver (Fungi fungi){

    // ---------------------------------------------------
    //   Set Up Static Mycelia
    // ---------------------------------------------------
    // Initialize starting parameters
    fungi.UpdateFungi_FromConfigFile() ;
    double init_hyphae_count = fungi.init_Count;   // Number of initial hyphae segments
    
    //vector<HyphaeSegment> fungi.hyphaeSegments;       // Vector of classes where hyphae info is stored
    HyphaeSegment new_hy;

    // NOTE:
    // The i^th hyphae segment is considered a tip if 
    //      hy[i].can_extend == true
    // The location of the tip is given by 
    //      hy[i].x2 & hy[i].y2 

    // Set up initial mycelia
    cout << " " << endl;
    cout << "INITIAL HYPHAE:" << endl;
    cout << "---------------" << endl;
    double initAngle =  pi ; // pi / 4.0 ;
    for ( int i = 0; i < init_hyphae_count; i++ ) {
        
        // Starting points and angles
        
        new_hy.x1 = fungi.initX ;            //60                                          // One endpoint of i^th segment: x-coordinate
        new_hy.y1 = fungi.initY;             //60                                         // One endpoint of i^th segment: y-coordinate
        new_hy.angle = i * ( ( 2 * pi ) / init_hyphae_count ) + initAngle ;  // Angle of i^th segment
        new_hy.x2 = new_hy.x1 + new_hy.hyphaeLength * cos( new_hy.angle );                      // Other endpoint of i^th segment: x-coordinate
        new_hy.y2 = new_hy.y1 + new_hy.hyphaeLength * sin( new_hy.angle );                      // Other endpoint of i^th segment: y-coordinate
        new_hy.can_branch = true;                                            // The i^th segment can extend
        new_hy.can_extend = true;                                            // The i^th segment can branch
        new_hy.from_who = -1 ;
        new_hy.extend_to = -1 ;
        new_hy.branch_to = -1 ;
        
        fungi.hyphaeSegments.push_back( new_hy );                             // Append hy vector for new segment

        // Print out info:
        cout << "Segment " << i << endl;
        cout << "   Endpoints of new segment:       (" << fungi.hyphaeSegments[i].x1 << ", " << fungi.hyphaeSegments[i].y1 << ") to (" << fungi.hyphaeSegments[i].x2 << ", " << fungi.hyphaeSegments[i].y2 << ")" << endl;
        cout << "   Angle of new segment (radians): " << fungi.hyphaeSegments[i].angle << endl;
    }
      
    // Define extra segments (make the shape non-uniform/non-circular)
    cout << " " << endl;
    cout << "EXTENDED/BRANCHED HYPHAE:" << endl;
    cout << "-------------------------" << endl;
    //add_new_hyphae( fungi.hyphaeSegments, 0, fungi.extendingAngle, "extend");      //4th
    //add_new_hyphae( fungi.hyphaeSegments, 1, fungi.extendingAngle, "extend");      //4th
    //add_new_hyphae( fungi.hyphaeSegments, 1, ( pi / 3.0 ), "branch");             //7th
    /*
    add_new_hyphae( fungi.hyphaeSegments, 0, fungi.extendingAngle, "extend");      //4th
    add_new_hyphae( fungi.hyphaeSegments, 2, 0, "extend");                      //5th
    add_new_hyphae( fungi.hyphaeSegments, 5, 0, "extend");                      //6th
    add_new_hyphae( fungi.hyphaeSegments, 5, ( pi / 3.0 ), "branch");             //7th
    
    add_new_hyphae(fungi.hyphaeSegments, 1, 0, "extend") ;                   //8th
    add_new_hyphae(fungi.hyphaeSegments, 1, pi/3.0, "branch") ;                  //9th
    add_new_hyphae(fungi.hyphaeSegments, 9, 0, "extend") ;                      //10th
   // add_new_hyphae(fungi.hyphaeSegments, 10, 0, "extend") ;                      //11th
    
    add_new_hyphae(fungi.hyphaeSegments, 3, -pi/3.0, "branch") ;                  //12th
    
   // add_new_hyphae( fungi.hyphaeSegments, 7, 0, "extend");                      //13th
    
    add_new_hyphae( fungi.hyphaeSegments, 2, -pi/3.0 , "branch");                      //14th
    add_new_hyphae( fungi.hyphaeSegments, 12, 0, "extend");                      //15th
   // add_new_hyphae( fungi.hyphaeSegments, 15, pi/3.0 , "branch");                      //16th
   // add_new_hyphae( fungi.hyphaeSegments, 6, 0, "extend");                      //17th
    */
    cout << " " << endl;
    /*
    add_new_hyphae( fungi.hyphaeSegments, 0, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 1, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 2, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 3, 0, "extend");
    */
    

    //fungi.Find_Hyphae_Tips() ;
    fungi.Find_Hyphae_Tips2() ;         //branch points are not tip
    fungi.FindCenterPointConnection() ;
    
    fungi.FindProductionNetwork() ;

    return fungi ;
    
}
