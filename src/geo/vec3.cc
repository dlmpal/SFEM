#include "vec3.h"
#include <cmath>

namespace sfem::geo
{
    //=============================================================================
    Vec3::Vec3(Scalar _x, Scalar _y, Scalar _z)
    {
        x[0] = _x;
        x[1] = _y;
        x[2] = _z;
    }
    //=============================================================================
    Vec3::Vec3(Scalar x1, Scalar y1, Scalar z1, Scalar x2, Scalar y2, Scalar z2)
    {
        x[0] = x2 - x1;
        x[1] = y2 - y1;
        x[2] = z2 - z1;
    }
    //=============================================================================
    Vec3 Vec3::operator+(const Vec3 &other) const
    {
        return Vec3(this->x[0] + other.x[0], this->x[1] + other.x[1], this->x[2] + other.x[2]);
    }
    //=============================================================================
    Vec3 Vec3::operator-(const Vec3 &other) const
    {
        return Vec3(this->x[0] - other.x[0], this->x[1] - other.x[1], this->x[2] - other.x[2]);
    }
    //=============================================================================
    Scalar Vec3::operator*(const Vec3 &other) const
    {
        return this->x[0] * other.x[0] + this->x[1] * other.x[1] + this->x[2] * other.x[2];
    }
    //=============================================================================
    Vec3 Vec3::operator*(Scalar a) const
    {
        return Vec3(this->x[0] * a, this->x[1] * a, this->x[2] * a);
    }
    //=============================================================================
    Scalar Vec3::mag() const
    {
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }
    //=============================================================================
    Vec3 Vec3::cross_prod(const Vec3 &other) const
    {
        return Vec3(this->x[1] * other.x[2] - this->x[2] * other.x[1],
                    this->x[2] * other.x[0] - this->x[0] * other.x[2],
                    this->x[0] * other.x[1] - this->x[1] * other.x[0]);
    }
    //=============================================================================
    Vec3 Vec3::normalize() const
    {
        Scalar _mag = mag();
        return Vec3(x[0] / _mag, x[1] / _mag, x[2] / _mag);
    }
    //=============================================================================
    Vec3 Vec3::unit_tan() const
    {
        return normalize();
    }
    //=============================================================================
    Vec3 Vec3::unit_norm() const
    {
        Scalar _mag = mag();

        Scalar xn = x[1] / _mag;
        Scalar yn = -x[0] / _mag;
        Scalar zn = x[2] / _mag;

        return Vec3(xn, yn, zn);
    }
}