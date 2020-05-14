
#include "Bacteria.hpp"

bacterium::bacterium ()
{
    nodes.resize(nnode) ;
    ljnodes.resize(nnode-1) ;
    allnodes.resize(2*nnode-1) ;
    duplicate.resize(nnode) ;
    connection.resize(nbacteria) ;
    
    pili.resize(nPili) ;
}
