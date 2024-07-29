#pragma once

#include "../common/config.h"
#include <vector>

namespace sfem::la
{
    class DenseMatrix
    {
    public:
        /// @brief Create a DenseMatrix
        /// @param n_rows Number of rows
        /// @param n_cols Number of columns
        /// @param value  Value to initialize all entries
        DenseMatrix(int n_rows, int n_cols, Scalar value = 0.0);

        /// @brief Create a DenseMatrix
        /// @param n_rows Number of rows
        /// @param n_cols Number of columns
        DenseMatrix(int n_rows, int n_cols, const std::vector<Scalar> &entries);

        /// @brief Get the number of rows
        int n_rows() const;

        /// @brief Get the number of columns
        int n_cols() const;

        /// @brief Get a reference to the matrix entries, e.g.
        /// the underlying std::vector
        const std::vector<Scalar> &entries() const;

        /// @brief Set all the entries uniformly
        void set_all(Scalar value);

        /// @brief Reurn the entry at the i, j position
        Scalar at(int i, int j) const;

        /// @brief Get the values of the i-th row
        std::vector<Scalar> get_row_values(int i) const;

        /// @brief Get the values of the j-th column
        std::vector<Scalar> get_col_values(int j) const;

        /// @brief Insert a value at the i, j position
        void insert(int i, int j, Scalar value);

        /// @brief Add a value to the i, j position
        void add(int i, int j, Scalar value);

        /// @brief Return the transpose
        DenseMatrix T() const;

        /// @brief Matrix multiplication
        DenseMatrix operator*(const DenseMatrix &other) const;

        /// @brief Scale the matrix by a constant alpha
        DenseMatrix operator*(Scalar alpha) const;

        /// @brief Matrix addition
        DenseMatrix operator+(const DenseMatrix &other) const;

        /// @brief Matrix addition (inplace)
        void operator+=(const DenseMatrix &other);

        /// @brief Add a constant alpha to the matrix
        DenseMatrix operator+(Scalar alpha) const;

        /// @brief Add a constant alpha to the matrix (inplace)
        void operator+=(Scalar alpha);

        /// @brief Matrix subtraction
        DenseMatrix operator-(const DenseMatrix &other) const;

        /// @brief Matrix subtraction (inplace)
        void operator-=(const DenseMatrix &other);

        /// @brief Subtract a constant alpha from the matrix
        DenseMatrix operator-(Scalar alpha) const;

        /// @brief Subtract a constant alpha from the matrix (inplace)
        void operator-=(Scalar alpha);

    private:
        /// @brief Number of rows
        int _n_rows;

        /// @brief Number of columns
        int _n_cols;

        /// @brief Entries
        std::vector<Scalar> _entries;
    };
}