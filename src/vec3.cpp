#include "vec3.hpp"
#include <cmath>


/* Vector Addition */
vec3 operator+(const vec3& a, const vec3& b) {
    return {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    };
}

/* Vector Subtraction */
vec3 operator-(const vec3& a, const vec3& b) {
    return {
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    };
}

/* Vector Multiplication with scalar */
vec3 operator*(double k, const vec3& v) {
    return {k * v.x, k * v.y, k * v.z};
}

vec3 operator*(const vec3& v, double k) {
    return {v.x * k, v.y * k, v.z * k};
}

/* Vector Division with scalar */
vec3 operator/(const vec3& v, double k) {
    return {v.x / k, v.y / k, v.z / k};
}

double norm(const vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Dot Product of Vectors
double dot(const vec3& a, const vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross Product of Vectors
vec3 cross(const vec3& a, const vec3& b) {
    return {
        a.y * b.z - b.y * a.z,
        a.z * b.x - a.x * b.z,
        a.x * b.y - b.x * a.y
    };
}

// Equality Operators
bool operator==(const vec3& a, const vec3& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool operator!=(const vec3& a, const vec3& b) {
    return !(a == b);
}

const vec3 vec3::zero = {0.0, 0.0, 0.0};