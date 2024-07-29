#pragma once

#include "logger.h"
#include <cmath>
#include <vector>

namespace sfem::math
{
    /// @brief Tranpose a matrix of size (r, c)
    inline void transpose(int r, int c, const Scalar m[], Scalar mt[])
    {
        for (auto i = 0; i < r; i++)
        {
            for (auto j = 0; j < c; j++)
            {
                mt[j * r + i] = m[i * c + j];
            }
        }
    }

    /// @brief Multiply two matrices of size (r, rc) and (rc, c)
    inline void matmult(int r, int c, int rc, const Scalar m1[], const Scalar m2[], Scalar m[])
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m[i * c + j] = 0.0;

                for (int k = 0; k < rc; k++)
                {
                    m[i * c + j] = std::fma(m1[i * rc + k], m2[k * c + j], m[i * c + j]);
                }
            }
        }
    }

    /// @brief Add two matrices of size (r, c)
    /// @note the matrix are scaled by a1 and a2, i.e. m = a1 * m1 + a2 * m2
    inline void matadd(int r, int c, const Scalar m1[], Scalar a1,
                       const Scalar m2[], Scalar a2, Scalar m[])
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m[i * c + j] = a1 * m1[i * c + j] + a2 * m2[i * c + j];
            }
        }
    }

    /// @brief Compute y += a*x
    inline void axpy(int size, Scalar a, const Scalar x[], Scalar y[])
    {
        for (int i = 0; i < size; i++)
        {
            y[i] = std::fma(a, x[i], y[i]);
        }
    }

    /// @brief Determinant of a square 3x3, 2x2 or 1x1 matrix
    /// @note The array m[] should always have size 3x3, regardless of r
    inline Scalar det(int r, const Scalar m[])
    {
        Scalar d = -1;

        if (r == 3)
        {
            d = m[0] * (m[4] * m[8] - m[7] * m[5]) -
                m[1] * (m[3] * m[8] - m[6] * m[5]) +
                m[2] * (m[3] * m[7] - m[4] * m[6]);
        }
        else if (r == 2)
        {
            d = m[0] * m[4] - m[1] * m[3];
        }
        else if (r == 1)
        {
            d = m[0];
        }
        else
        {
            Logger::instance().error("Determinant not defined for square matrix with dimension: " + std::to_string(r) + "\n", __FILE__, __LINE__);
        }

        return d;
    }

    /// @brief Inverse and determinant of a square 3x3, 2x2 or 1x1 matrix
    /// @note The arrays m[] and mi[] should always have size 3x3, regardless of r
    inline Scalar inv(int r, const Scalar m[], Scalar mi[])
    {
        // Determinant and inverse
        Scalar d = det(r, m);
        Scalar di = 1 / d;

        if (r == 3)
        {
            mi[0] = di * (m[4] * m[8] - m[5] * m[7]);
            mi[1] = -di * (m[1] * m[8] - m[2] * m[7]);
            mi[2] = di * (m[1] * m[5] - m[2] * m[4]);
            mi[3] = -di * (m[3] * m[8] - m[5] * m[6]);
            mi[4] = di * (m[0] * m[8] - m[2] * m[6]);
            mi[5] = -di * (m[0] * m[5] - m[2] * m[3]);
            mi[6] = di * (m[3] * m[7] - m[4] * m[6]);
            mi[7] = -di * (m[0] * m[7] - m[1] * m[6]);
            mi[8] = di * (m[0] * m[4] - m[1] * m[3]);
        }
        else if (r == 2)
        {
            mi[0] = di * m[4];
            mi[1] = -di * m[1];
            mi[3] = -di * m[3];
            mi[4] = di * m[0];
        }
        else if (r == 1)
        {
            mi[0] = di;
        }
        else
        {
            Logger::instance().error("Inverse not defined for square matrix with dimension: " + std::to_string(r) + "\n", __FILE__, __LINE__);
        }

        return d;
    }

    /// @brief Moore-Penrose pseudo-inverse of a 3x2, 3x1 or 2x1 matrix
    /// @note The arrays m[] and mi[] should always have size 3x3, regardless of r and c
    inline Scalar pinv(int r, int c, const Scalar m[], Scalar mi[])
    {
        // Transpose
        Scalar mt[3 * 3] = {0};
        transpose(3, 3, m, mt);

        // Intermediate product
        Scalar mtm[3 * 3] = {0};
        matmult(3, 3, 3, mt, m, mtm);

        // Invert intermediate product
        Scalar mtmi[3 * 3] = {0};
        inv(c, mtm, mtmi);

        // Multiply by tranpose
        matmult(3, 3, 3, mtmi, mt, mi);

        // Compute the determinant
        Scalar d = sqrt(det(c, mtm));

        return d;
    }
}
