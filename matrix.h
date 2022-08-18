#ifndef MATRIX2_H
#define MATRIX2_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

template <class T>
class matrix2
{
public:
    // the various constructors
    matrix2();
    matrix2(int nRows, int nCols);
    matrix2(int nRows, int nCols, const T *inputData);
    matrix2(const matrix2<T> &inputMatrix);
    matrix2(int nRows, int nCols, const std::vector<T> *inputData);

    // destructor
    ~matrix2();

    // configuration methods
    bool resize(int numRows, int numCols);

    // element access methods
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // overload the == operator
    bool operator==(const matrix2<T> &rhs);

    // overload +, -, and * operators (friends)
    template <class U>
    friend matrix2<U> operator+(const matrix2<U> &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator+(const U &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator+(const matrix2<U> &lhs, const U &rhs);

    template <class U>
    friend matrix2<U> operator-(const matrix2<U> &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator-(const U &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator-(const matrix2<U> &lhs, const U &rhs);

    template <class U>
    friend matrix2<U> operator*(const matrix2<U> &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator*(const U &lhs, const matrix2<U> &rhs);
    template <class U>
    friend matrix2<U> operator*(const matrix2<U> &lhs, const U &rhs);

private:
    int Sub2Ind(int row, int col);
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

/* ************************************************************************************************ */
// constructors and destructors

// the default constructor
template <class T>
matrix2<T>::matrix2()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

// construct an empty matrix (all elements 0)
template <class T>
matrix2<T>::matrix2(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = 0.0;
}

// construct from const linear array
template <class T>
matrix2<T>::matrix2(int nRows, int nCols, const T *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData[i];
}

// the copy constructor
template <class T>
matrix2<T>::matrix2(const matrix2<T> &inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// construct from std::vector
template <class T>
matrix2<T>::matrix2(int nRows, int nCols, const std::vector<T> *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData->at(i);
}

// destructor
template <class T>
matrix2<T>::~matrix2()
{
    if (m_matrixData != nullptr)
        delete[] m_matrixData;
}

/* ************************************************************************************************ */
// configuration functions

template <class T>
bool matrix2<T>::resize(int numRows, int numCols)
{
    m_nRows = numRows;
    m_nCols = numCols;
    m_nElements = (m_nRows * m_nCols);
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr)
    {
        for (int i = 0; i < m_nElements; i++)
            m_matrixData[i] = 0.0;
        return true;
    }
    else
        return false;
}

/* ************************************************************************************************ */
// element functions

template <class T>
T matrix2<T>::GetElement(int row, int col)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
        return m_matrixData[linearIndex];
    else
        return 0;
}

template <class T>
bool matrix2<T>::SetElement(int row, int col, T elementValue)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        m_matrixData[linearIndex] = elementValue;
        return true;
    }
    else
        return false;
}

template <class T>
int matrix2<T>::GetNumRows()
{
    return m_nRows;
}

template <class T>
int matrix2<T>::GetNumCols()
{
    return m_nCols;
}

/* ************************************************************************************************ */
// overloaded operator functions

/* the + operator ********************************************************************************* */
// matrix + matrix
template <class T>
matrix2<T> operator+(const matrix2<T> &lhs, const matrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar + matrix
template <class T>
matrix2<T> operator+(const T &lhs, const matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs + rhs.m_matrixData[i];

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix + scalar
template <class T>
matrix2<T> operator+(const matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] + rhs;

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

/* the - operator ********************************************************************************* */
// matrix - matrix
template <class T>
matrix2<T> operator-(const matrix2<T> &lhs, const matrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar - matrix
template <class T>
matrix2<T> operator-(const T &lhs, const matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs - rhs.m_matrixData[i];

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix - scalar
template <class T>
matrix2<T> operator-(const matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] - rhs;

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

/* the * operator ********************************************************************************* */
// scalar * matrix
template <class T>
matrix2<T> operator*(const T &lhs, const matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs * rhs.m_matrixData[i];

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * scalar
template <class T>
matrix2<T> operator*(const matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] * rhs;

    matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * matrix
template <class T>
matrix2<T> operator*(const matrix2<T> &lhs, const matrix2<T> &rhs)
{
    int l_numRows = lhs.m_nRows;
    int l_numCols = lhs.m_nCols;
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;

    // standard matrix multiplication condition
    if (l_numCols == r_numRows)
    {
        T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];

        // loop through each row of lhs
        for (int lhsRow = 0; lhsRow < l_numRows; lhsRow++)
        {
            // loop through each column of rhs
            for (int rhsCol = 0; rhsCol < r_numCols; rhsCol++)
            {
                T elementResult = 0.0;

                // loop through each element of this lhs row
                for (int lhsCol = 0; lhsCol < l_numCols; lhsCol++)
                {
                    // compute the lhs linear index
                    int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
                    // compute the rhs linear index
                    int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
                    // perform the calculation on these elements
                    elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
                }

                // store the result
                int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
                tempResult[resultLinearIndex] = elementResult;
            }
        }

        matrix2<T> result(l_numRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    }
    else
    {
        matrix2<T> result(1, 1);
        return result;
    }
}

/* the == operator ******************************************************************************** */

template <class T>
bool matrix2<T>::operator==(const matrix2<T> &rhs)
{
    // check if matrices are the same size, if not return false
    if ((this->m_nRows != rhs.m_nRows) || (this->m_nCols != rhs.m_nCols))
        return false;

    // check if the elements are equal
    bool flag = true;
    for (int i = 0; i < this->m_nElements; i++)
    {
        if (this->m_matrixData[i] != rhs.m_matrixData[i])
            flag = false;
    }
    return flag;
}

/* ************************************************************************************************ */
// private functions

template <class T>
int matrix2<T>::Sub2Ind(int row, int col)
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return (row * m_nCols) + col;
    else
        return -1;
}

#endif