#ifndef GROWTHFUNCTIONS_H_
#define GROWTHFUNCTIONS_H_

#include <tuple>
#include <vector>
#include "hyphaeSegment.h"
#include "constants.h"
using namespace std;

class HyphaeSegment;

void add_new_hyphae( vector<HyphaeSegment> &hy, int from_idx, double angle_new, string growth_type );

#endif
