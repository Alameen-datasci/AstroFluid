#pragma once

struct vec3 {
    double x, y, z;

    // Default constructor (automatically initializes to zero)
    vec3() : x(0.0), y(0.0), z(0.0) {}
    
    // Parameterized constructor
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Static constant for semantic readability
    static const vec3 zero;
};

vec3 operator+(const vec3&, const vec3&);
vec3 operator-(const vec3&, const vec3&);
vec3 operator*(const vec3&, double);
vec3 operator*(double, const vec3&);
vec3 operator/(const vec3&, double);

bool operator==(const vec3&, const vec3&);
bool operator!=(const vec3&, const vec3&);

double norm(const vec3&);
double dot(const vec3&, const vec3&);
vec3 cross(const vec3&, const vec3&);