#pragma once
#include <math.h>
#include "GL/glut.h"

// A 2-D vector
class Vec2 { 
public:
    Vec2(GLfloat x, GLfloat y) { // Constructor
        this->x = x;
        this->y = y;
    }

    Vec2(const Vec2& vec) {    // Copy constructor
        this->x = vec.x;
        this->y = vec.y;
    }

    GLfloat getx()
    {
        return this->x;
    }

    GLfloat gety()
    {
        return this->y;
    }

    float Length() {               // Return length of vector
        return sqrt(x * x + y * y);
    }

    void Normalize() {             // Transform into normal vector in same direction
        float len = this->Length();
        this->x /= len;
        this->y /= len;
    }

    void Perpendicu() {            // Transform into vector perpendicular to current orientation
        float _x = this->x;
        this->x = -this->y;
        this->y = _x;
    }

    Vec2 operator+(const Vec2& b) {     // Overload + Operator
        return Vec2(this->x + b.x, this->y + b.y);

    }

    Vec2 operator-(const Vec2& b) {     // Overload - Operator
        return Vec2(this->x - b.x, this->y - b.y);

    }

    Vec2 operator*(const GLint c) {     // Overload * Operator for * by Constants
        float constant = float(c);
        return Vec2(this->x * c, this->y * c);
    }

private:
    //Define vector
    GLfloat x = 0.0;
    GLfloat y = 0.0;
};