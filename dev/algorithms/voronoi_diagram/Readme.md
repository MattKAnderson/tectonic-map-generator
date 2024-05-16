# Voronoi Diagram

This is the directory for development of code to produce voronoi diagrams.

# status

Fortunes Algorithm implemented for producing a voronoi diagram in O(log N) time. Runs in about 2ms for 1000 site locations.

# TODOs
## Features / bug fixes
1. handle special case where new site.y equals the current intersection of two parabolae (correctness improvement... is this even possible? I'm not seeing it at tests of 10M+ sites)
1. clip to a smaller window

## optimization opportunities
1. finish implementing the beach line tree as a red-black tree 
1. Implement memory managers for 
   a. half-edges 
   a. vertices 
   a. arcs 
1. investigate why creating the vertex/region graphs is taking as long as it is