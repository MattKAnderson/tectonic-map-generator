#include <geometry.hpp>

class LineClipper {
public:
    LineClipper(double min, double max): min_(min), max_(max) {}
    bool CohenSutherlandClip(const RealCoordinate& a, const RealCoordinate& b);
    RealCoordinate get_clipped_a() { return a_; };
    RealCoordinate get_clipped_b() { return b_; };
private:
    double min_;
    double max_;
    RealCoordinate a_;
    RealCoordinate b_;
    int outcode(const RealCoordinate& c);
    const static int INSIDE = 0;
    const static int LEFT = 1;
    const static int RIGHT = 2;
    const static int BOTTOM = 4;
    const static int TOP = 8;
};