# Triangle circumcenter computation

Calculating the circumcenter of a triangle given it's vertices is one of the core geometry operations in Fortune's algorithm. This calculation is used to identify where intersection events would occur between 3 prospective regions. 

For N seeds, randomly spread over a square region, there is approximately 5.9N number of triangle circumcenter computations (determined empirically). 

This subdirectory hosts code to test the performance of two methods of computing the circumcenter. The Perpendicular Bisector method and a method derived from solving the linear system of equations formed by the constraint that the intersect is equidistant from all 3 of the triangles vertices.

# Results

The linear system method was faster to compute. This is due to requiring only a single division operation when computing the denominator of the 2x2 matrix determinant.

Average compute time and clock cycles on an Intel 13th gen processor:

| method | time (ns) | clock cycles |
|--------|-----------|--------------|
| Linear System | 22.6 | 116 |
| Perpendicular Bisector | 28.6 | 146 |
