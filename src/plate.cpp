#include "plate.hpp"

Vector2D<double> plate::movement(double timeunits) {
    return timeunits * drift;
}

void plate::update_drift(Vector2D<double> drift) {
    this->drift = drift;
}
