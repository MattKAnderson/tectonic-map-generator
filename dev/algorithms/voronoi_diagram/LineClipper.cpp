#include <LineClipper.hpp>

int LineClipper::outcode(const RealCoordinate& c) {
    int code = INSIDE;
    if (c.x < min_) {
        code |= LEFT;
    }
    else if (c.x > max_) {
        code |= RIGHT;
    }
    if (c.y < min_) {
        code |= BOTTOM;
    }
    else if (c.y > max_) {
        code |= TOP;
    }
    return code;
}

bool LineClipper::CohenSutherlandClip(
    const RealCoordinate& a, const RealCoordinate& b
) {
    a_ = a;
    b_ = b;
    int code_a = outcode(a_);
    int code_b = outcode(b_);
    while(true) {
        if (!(code_a | code_b)) { return true; }
        else if (code_a & code_b) { return false; }
        
        double x, y;
        int code_out = code_a > code_b ? code_a : code_b;
        if (code_out & TOP) {
            x = a_.x + (b_.x - a_.x) * (max_ - a_.y) / (b_.y - a_.y);
            y = max_;
        }
        else if (code_out & BOTTOM) {
            x = a_.x + (b_.x - a_.x) * (min_ - a_.y) / (b_.y - a_.y);
            y = min_;
        }
        else if (code_out & RIGHT) {
            y = a_.y + (b_.y - a_.y) * (max_ - a_.x) / (b_.x - a_.x);
            x = max_;
        }
        else {
            y = a_.y + (b_.y - a_.y) * (min_ - a_.x) / (b_.x - a_.x);
            x = min_;
        }

        if (code_out == code_a) {
            a_ = {x, y};
            code_a = outcode(a_);
        }
        else {
            b_ = {x, y};
            code_b = outcode(b_);
        }
    }
}



