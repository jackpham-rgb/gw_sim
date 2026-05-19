// vec3.hpp -- minimal 3D vector type for the N-body simulation.
// Header-only. Plain double precision.

#pragma once

#include <cmath>

namespace gw {

struct Vec3 {
    double x, y, z;

    constexpr Vec3() : x(0), y(0), z(0) {}
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s)      const { return {x / s, y / s, z / s}; }
    Vec3 operator-()              const { return {-x, -y, -z}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s)      { x *= s; y *= s; z *= s; return *this; }

    double dot(const Vec3& o)  const { return x * o.x + y * o.y + z * o.z; }
    Vec3   cross(const Vec3& o) const {
        return { y * o.z - z * o.y,
                 z * o.x - x * o.z,
                 x * o.y - y * o.x };
    }
    double norm()  const { return std::sqrt(x * x + y * y + z * z); }
    double norm2() const { return x * x + y * y + z * z; }
    Vec3   unit()  const {
        double n = norm();
        return n > 0.0 ? Vec3{ x / n, y / n, z / n } : Vec3{};
    }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

} // namespace gw
