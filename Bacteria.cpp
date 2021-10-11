
#include "Bacteria.hpp"

bacterium::bacterium ()
{
    nodes.resize(nnode) ;
    ljnodes.resize(nnode-1) ;
    allnodes.resize(2*nnode-1) ;
    duplicate.resize(nnode) ;
    connection.resize(nbacteria) ;
    oldLoc.resize(2) ;
    
    pili.resize(nPili) ;
}
bool bacterium::SourceRegion()
{
    double sourceRadius = 10.0 ;
    double r2 = sqrt( pow( (nodes[(nnode-1)/2].x - domainx/2.0 ) ,2 ) + pow( (nodes[(nnode-1)/2].y - domainy/2.0 ) ,2 ) );
    if (r2< sourceRadius )
    {
        sourceWithin = true ;
        
    }
    //else is not written yet. once it is in source sourceWithin remains true forever
    
    return sourceWithin ;
}

