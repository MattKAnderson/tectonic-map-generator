# Voronoi Diagram

This is the directory for development of code to produce voronoi diagrams.

# status

Fortunes Algorithm implemented for producing a voronoi diagram in O(log N) time. Runs in about 1ms for 1000 site locations.

# TODOs
## Features
1. re-organize the code (put implementation details in another header)
1. Implement Lloyd relaxation
1. Produce visualizations, add some results to this Readme
    1. image of bounded, image of cropped 
    1. movie of algorithm progression
    1. movie of Lloyd relaxation iterations
1. compare to the boost implementation on my machine

## Bugs
None observed

## optimization opportunities
1. (low priority) finish implementing the beach line tree as a red-black tree. Low priority because after measurement, the find algorithm only represents ~15-18% of runtime.  
