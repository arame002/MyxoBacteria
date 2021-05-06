
#include "driver.h"
#include "constants.h"

using namespace std;
using constants::dx;
using constants::pi;


// ---------------------------------------------------
//   Main Driver Function
// ---------------------------------------------------
Fungi driver (){

    // ---------------------------------------------------
    //   Set Up Static Mycelia
    // ---------------------------------------------------
    // Initialize starting parameters
    double init_hyphae_count = 4;   // Number of initial hyphae segments
    
    Fungi fungi ;                   //Alireza
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
    for ( int i = 0; i < init_hyphae_count; i++ ) {
        
        // Starting points and angles
        fungi.hyphaeSegments.push_back( new_hy );                                             // Append hy vector for new segment
        fungi.hyphaeSegments[i].x1 = 50;                                                      // One endpoint of i^th segment: x-coordinate
        fungi.hyphaeSegments[i].y1 = 50;                                                      // One endpoint of i^th segment: y-coordinate
        fungi.hyphaeSegments[i].angle = i * ( ( 2 * pi ) / init_hyphae_count ) + ( pi / 4 );  // Angle of i^th segment
        fungi.hyphaeSegments[i].x2 = fungi.hyphaeSegments[i].x1 + dx * cos( fungi.hyphaeSegments[i].angle );                      // Other endpoint of i^th segment: x-coordinate
        fungi.hyphaeSegments[i].y2 = fungi.hyphaeSegments[i].y1 + dx * sin( fungi.hyphaeSegments[i].angle );                      // Other endpoint of i^th segment: y-coordinate
        fungi.hyphaeSegments[i].can_branch = true;                                            // The i^th segment can extend
        fungi.hyphaeSegments[i].can_extend = true;                                            // The i^th segment can branch

        // Print out info:
        cout << "Segment " << i << endl;
        cout << "   Endpoints of new segment:       (" << fungi.hyphaeSegments[i].x1 << ", " << fungi.hyphaeSegments[i].y1 << ") to (" << fungi.hyphaeSegments[i].x2 << ", " << fungi.hyphaeSegments[i].y2 << ")" << endl;
        cout << "   Angle of new segment (radians): " << fungi.hyphaeSegments[i].angle << endl;
    }
        
    // Define extra segments (make the shape non-uniform/non-circular)
    cout << " " << endl;
    cout << "EXTENDED/BRANCHED HYPHAE:" << endl;
    cout << "-------------------------" << endl;
    add_new_hyphae( fungi.hyphaeSegments, 0, constants::dTheta, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 2, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 5, 0, "extend");
    add_new_hyphae( fungi.hyphaeSegments, 5, ( pi / 3 ), "branch");
    cout << " " << endl;
    
    fungi.Find_Hyphae_Tips() ;

    return fungi ;
    
}
