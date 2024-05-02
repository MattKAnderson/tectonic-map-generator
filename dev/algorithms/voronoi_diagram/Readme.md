# Voronoi Diagram

This is the directory for development of code to produce voronoi diagrams.

# status

Fortunes Algorithm implemented for producing a voronoi diagram in O(log N) time. Runs in about 5ms for 1000 seeds locations. Stable implementation into millions of seeds, but does not handle some special cases surrounding numerical precision that can cause some erros when number of seeds gets that large.

# TODOs
1. implement edges between vertices clipped to the BBOX bounds
1. implement RegionNode and VertexNode, RegionNode should define a polygon (ordered list of vertices) and include adjacencies
1. handle special case where new site.y equals the current intersection of two parabolae (correctness improvement, for large number of seeds)
1. implement the beach line tree as a red-black tree (runtime performance improvement)