#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <tuple>
#include <vector>

#include "growthFunctions.h"
using namespace std;
using constants::pi;

void add_new_hyphae( vector<HyphaeSegment> &hy, int from_idx, double angle_new, string growth_type ) {

    /* 
     * Purpose of ADD_NEW_HYPHAE:
     *    Append the vector hy with information regarding a newly extended or branched hyphae segment
     * 
     * @param[in] hy : vector of class intances containing info about each hyphae segment
     * @param[in] from_idx : index of hyphae segment from which new segment extended or branched
     * @param[in] angle_new : how the angle of the new segment is modified
     * @param[in] growth_type : "extend" if the growth is due to extension, "branch" if growth is due to branching
    */

    // Add new segment
    HyphaeSegment new_hy;
    hy.push_back( new_hy );                     // Append vector of hyphae instances
    int new_idx = hy.size() - 1;                // Index for the newly extended/branched segment

    // Update neighbor information
    hy[new_idx].from_who = from_idx;            // Keep track of where new segment extended/branched from
    if ( growth_type == "extend" ) {
        hy[from_idx].extend_to = new_idx;       // Keep track of where originating segment extended to
        hy[from_idx].can_extend = false;        // A segment cannot extend more than once
    } else if ( growth_type == "branch" ) {
        hy[from_idx].branch_to = new_idx;       // Keep track of where originating segment branched to
        hy[from_idx].can_branch = false;        // A segment cannot branch more than once
    }
    hy[new_idx].can_branch = true;              // The new segment can branch
    hy[new_idx].can_extend = true;              // The new segment can extend

    // Upate position info for new segment
    hy[new_idx].x1 = hy[from_idx].x2;                                   // Starting endpoint of new segment (x-coorindate)
    hy[new_idx].y1 = hy[from_idx].y2;                                   // Starting endpoint of new segment (y-coorindate)
    hy[new_idx].angle = hy[from_idx].angle + angle_new;                 // Angle of the new segment
    hy[new_idx].x2 = hy[new_idx].x1 + hy[new_idx].hyphaeLength * cos( hy[new_idx].angle );    // Ending endpoint of new segment (x-coorindate)
    hy[new_idx].y2 = hy[new_idx].y1 + hy[new_idx].hyphaeLength * sin( hy[new_idx].angle );    // Ending endpoint of new segment (y-coorindate)

    // Print out info:
    cout << "Segment " << from_idx << " " << growth_type << " to " << new_idx << endl;
    cout << "   Endpoints of new segment:       (" << hy[new_idx].x1 << ", " << hy[new_idx].y1 << ") to (" << hy[new_idx].x2 << ", " << hy[new_idx].y2 << ")" << endl;
    cout << "   Angle of new segment (radians): " << hy[new_idx].angle << endl;
}
