#include <iostream>
#include <vector>
#include <string>

#include "Vector2D.hpp"

class plate {
public:
    plate() : drift({0.0, 0.0}) {};
    plate(Vector2D<double> drift) : drift(drift) {};
    Vector2D<double> movement(double timeunits);
    void update_drift(Vector2D<double> drift);
    Vector2D<double> drift;
private:
    //Vector2D<double> drift;
};