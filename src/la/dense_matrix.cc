#include "dense_matrix.h"
#include "../common/logger.h"
#include "../common/error.h"
#include "../common/math.h"

namespace sfem::la
{
    //=============================================================================
    void CheckDenseMatrixSize(int n_rows, int n_cols)
    {
        if (n_rows <= 0 || n_cols <= 0)
        {
            Logger::instance().error("Invalid DenseMatrix size: " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + "\n",
                                     __FILE__, __LINE__);
        }
    }
    //=============================================================================
    void CheckDenseMatrixBounds(int i, int j, int n_rows, int n_cols)
    {
        if (i < 0 || i >= n_rows || j < 0 || j >= n_cols)
        {
            Logger::instance().error("Invalid indices (i, j)=(" + std::to_string(i) + ", " + std::to_string(j) + ") for DenseMatrix of size : " + std::to_string(n_rows) + " x " + std::to_string(n_cols) + "\n",
                                     __FILE__, __LINE__);
        }
    }
    //=============================================================================
    void CheckDenseMatrixSizesMult(int n_rows_1, int n_cols_1, int n_rows_2, int n_cols_2)
    {
        if (n_cols_1 != n_rows_2)
        {
            std::string message = "Cannot multiply matrices of size ";
            message += std::to_string(n_rows_1) + " x " + std::to_string(n_cols_1);
            message += " and ";
            message += std::to_string(n_rows_2) + " x " + std::to_string(n_cols_2) + "\n";
            Logger::instance().error(message, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    void CheckDenseMatrixSizesAdd(int n_rows_1, int n_cols_1, int n_rows_2, int n_cols_2)
    {
        if (n_rows_1 != n_rows_2 || n_cols_1 != n_cols_2)
        {
            std::string message = "Cannot add matrices of size ";
            message += std::to_string(n_rows_1) + " x " + std::to_string(n_cols_1);
            message += " and ";
            message += std::to_string(n_rows_2) + " x " + std::to_string(n_cols_2) + "\n";
            Logger::instance().error(message, __FILE__, __LINE__);
        }
    }
    //=============================================================================
    DenseMatrix::DenseMatrix(int n_rows, int n_cols, Scalar value)
        : _n_rows(n_rows), _n_cols(n_cols)
    {
        CheckDenseMatrixSize(n_rows, n_cols);
        _entries.resize(n_rows * n_cols, value);
    }
    //=============================================================================
    DenseMatrix::DenseMatrix(int n_rows, int n_cols, const std::vector<Scalar> &entries)
        : _n_rows(n_rows), _n_cols(n_cols), _entries(entries)
    {
        CheckDenseMatrixSize(n_rows, n_cols);
        if (static_cast<std::size_t>(n_rows * n_cols) != entries.size())
        {
            error::invalid_size_error(n_rows * n_cols, entries.size(), __FILE__, __LINE__);
        }
    }
    //=============================================================================
    int DenseMatrix::n_rows() const
    {
        return _n_rows;
    }
    //=============================================================================
    int DenseMatrix::n_cols() const
    {
        return _n_cols;
    }
    //=============================================================================
    const std::vector<Scalar> &DenseMatrix::entries() const
    {
        return _entries;
    }
    //=============================================================================
    void DenseMatrix::set_all(Scalar value)
    {
        std::fill(_entries.begin(), _entries.end(), value);
    }
    //=============================================================================
    Scalar DenseMatrix::at(int i, int j) const
    {
        CheckDenseMatrixBounds(i, j, _n_rows, _n_cols);
        return _entries[i * _n_cols + j];
    }
    //=============================================================================
    std::vector<Scalar> DenseMatrix::get_row_values(int i) const
    {
        std::vector<Scalar> row_values(_n_cols);
        for (int j = 0; j < _n_cols; j++)
        {
            row_values[j] = at(i, j);
        }
        return row_values;
    }
    //=============================================================================
    std::vector<Scalar> DenseMatrix::get_col_values(int j) const
    {
        std::vector<Scalar> col_values(_n_rows);
        for (int i = 0; i < _n_rows; i++)
        {
            col_values[i] = at(i, j);
        }
        return col_values;
    }
    //=============================================================================
    void DenseMatrix::insert(int i, int j, Scalar value)
    {
        CheckDenseMatrixBounds(i, j, _n_rows, _n_cols);
        _entries[i * _n_cols + j] = value;
    }
    //=============================================================================
    void DenseMatrix::add(int i, int j, Scalar value)
    {
        CheckDenseMatrixBounds(i, j, _n_rows, _n_cols);
        _entries[i * _n_cols + j] += value;
    }
    //=============================================================================
    DenseMatrix DenseMatrix::T() const
    {
        DenseMatrix trans(_n_cols, _n_rows);
        math::transpose(_n_rows, _n_cols, _entries.data(), trans._entries.data());
        return trans;
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator*(const DenseMatrix &other) const
    {
        CheckDenseMatrixSizesMult(_n_rows, _n_cols, other._n_rows, other._n_cols);
        DenseMatrix result(_n_rows, other._n_cols);
        math::matmult(_n_rows, other._n_cols, _n_cols, _entries.data(), other._entries.data(), result._entries.data());
        return result;
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator*(Scalar alpha) const
    {
        DenseMatrix result(_n_rows, _n_cols);
        math::axpy(_n_rows * _n_cols, alpha, _entries.data(), result._entries.data());
        return result;
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator+(const DenseMatrix &other) const
    {
        CheckDenseMatrixSizesAdd(_n_rows, _n_cols, other._n_rows, other._n_cols);
        DenseMatrix result(_n_rows, _n_cols);
        math::matadd(_n_rows, _n_cols, _entries.data(), 1, other._entries.data(), 1, result._entries.data());
        return result;
    }
    //=============================================================================
    void DenseMatrix::operator+=(const DenseMatrix &other)
    {
        CheckDenseMatrixSizesAdd(_n_rows, _n_cols, other._n_rows, other._n_cols);
        math::matadd(_n_rows, _n_cols, _entries.data(), 1, other._entries.data(), 1, _entries.data());
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator+(Scalar alpha) const
    {
        DenseMatrix result(_n_rows, _n_cols);
        for (std::size_t i = 0; i < _entries.size(); i++)
        {
            result._entries[i] = _entries[i] + alpha;
        }
        return result;
    }
    //=============================================================================
    void DenseMatrix::operator+=(Scalar alpha)
    {
        for (std::size_t i = 0; i < _entries.size(); i++)
        {
            _entries[i] += alpha;
        }
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator-(const DenseMatrix &other) const
    {
        CheckDenseMatrixSizesAdd(_n_rows, _n_cols, other._n_rows, other._n_cols);
        DenseMatrix result(_n_rows, _n_cols);
        math::matadd(_n_rows, _n_cols, _entries.data(), 1, other._entries.data(), -1, result._entries.data());
        return result;
    }
    //=============================================================================
    void DenseMatrix::operator-=(const DenseMatrix &other)
    {
        CheckDenseMatrixSizesAdd(_n_rows, _n_cols, other._n_rows, other._n_cols);
        math::matadd(_n_rows, _n_cols, _entries.data(), 1, other._entries.data(), -1, _entries.data());
    }
    //=============================================================================
    DenseMatrix DenseMatrix::operator-(Scalar alpha) const
    {
        DenseMatrix result(_n_rows, _n_cols);
        for (std::size_t i = 0; i < _entries.size(); i++)
        {
            result._entries[i] = _entries[i] - alpha;
        }
        return result;
    }
    //=============================================================================
    void DenseMatrix::operator-=(Scalar alpha)
    {
        for (std::size_t i = 0; i < _entries.size(); i++)
        {
            _entries[i] -= alpha;
        }
    }
}