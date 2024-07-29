#pragma once

#include "../common/config.h"
#include <array>

namespace sfem::geo
{
    /// @brief 3D vector
    struct Vec3
    {
        // Default constructor
        Vec3() = default;

        /// @brief Create a Vec3 equal to [x, y, z]
        Vec3(Scalar x, Scalar y, Scalar z);

        /// @brief Create a Vec3 from (x1, y1, z1) to (x2, y2, z2)
        Vec3(Scalar x1, Scalar y1, Scalar z1, Scalar x2, Scalar y2, Scalar z2);

        /// @brief Add two vectors
        Vec3 operator+(const Vec3 &other) const;

        /// @brief Subtract two vectors
        Vec3 operator-(const Vec3 &other) const;

        /// @brief Inner product between two vectors
        Scalar operator*(const Vec3 &other) const;

        /// @brief Return a vector scaled by a
        Vec3 operator*(Scalar a) const;

        /// @brief Vector magnitude
        Scalar mag() const;

        /// @brief 3D vector cross product
        Vec3 cross_prod(const Vec3 &other) const;

        /// @brief Return the normalized vector
        Vec3 normalize() const;

        /// @brief Return a unit vector in the direction of Vec3
        Vec3 unit_tan() const;

        /// @brief Returns a unit vector, normal to the original, in the 2D sense
        Vec3 unit_norm() const;

        /// @brief Vector xyz components
        std::array<Scalar, 3> x = {0};
    };
}