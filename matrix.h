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
    bool Resize(int numRows, int numCols);
    void SetToIdentity();

    // element access methods
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // manipulation methods
    // compute matrix inverse
    bool Inverse();

    // overload the == operator
    bool operator==(const matrix2<T> &rhs);
    bool Compare(const matrix2<T> &matrix1, double tolerance);

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

    bool Separate(matrix2<T> *matrix1, matrix2<T> *matrix2, int colNum);

private:
    // public:
    bool Join(const matrix2<T> &matrix2);
    int Sub2Ind(int row, int col);
    bool IsSquare();
    bool CloseEnough(T f1, T f2);
    void SwapRow(int i, int j);
    void MultAdd(int i, int j, T multFactor);
    void MultRow(int i, T multFactor);
    int FindRowWithMaxElement(int colNumber, int startingRow);
    void PrintMatrix();

private:
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
bool matrix2<T>::Resize(int numRows, int numCols)
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

// function to convert the existing matrix into an identity matrix
template <class T>
void matrix2<T>::SetToIdentity()
{
    if (!IsSquare())
        throw std::invalid_argument("Cannot form an identity matrix that is not square!");

    for (int row = 0; row < m_nRows; row++)
    {
        for (int col = 0; col < m_nCols; col++)
        {
            if (col == row)
                m_matrixData[Sub2Ind(row, col)] = 1.0;
            else
                m_matrixData[Sub2Ind(row, col)] = 0.0;
        }
    }
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

template <class T>
bool matrix2<T>::Compare(const matrix2<T> &matrix1, double tolerance)
{
    // check that the matrices have the same dimensions
    int numRows1 = matrix1.m_nRows;
    int numCols1 = matrix1.m_nCols;
    if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
        return false;

    // loop over all the elements and compute the sum of differences
    double cumulativeSum = 0.0;
    for (int i = 0; i < m_nElements; i++)
    {
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += ((element1 - element2) * (element1 - element2));
    }

    // calculate the root mean square
    double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1) - 1));
    if (final < tolerance)
        return true;
    else
        return false;
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
        // if (this->m_matrixData[i] != rhs.m_matrixData[i])
        if (!CloseEnough(this->m_matrixData[i], rhs.m_matrixData[i]))
            flag = false;
    }
    return flag;
}

// separate the matrix into two parts, around the column number provided
// note that the output is returned into the two matrix2<T> pointers in the input argument list
template <class T>
bool matrix2<T>::Separate(matrix2<T> *matrix1, matrix2<T> *matrix2, int colNum)
{
    // compute the sizes of the new matrices
    int numRows = m_nRows;
    int numCols1 = colNum;
    int numCols2 = m_nCols - colNum;

    // resize the two matrices to the proper dimensions
    matrix1->Resize(numRows, numCols1);
    matrix2->Resize(numRows, numCols2);

    // loop over the original matrix and store data into
    // the appropriate elements of the two output matrices
    for (int row = 0; row < m_nRows; row++)
    {
        for (int col = 0; col < m_nCols; col++)
        {
            if (col < colNum)
                matrix1->SetElement(row, col, this->GetElement(row, col));
            else
                matrix2->SetElement(row, col - colNum, this->GetElement(row, col));
        }
    }
}

/* ************************************************************************************************ */
// private functions

// join two matrices
template <class T>
bool matrix2<T>::Join(const matrix2<T> &matrix2)
{
    // extract the information that we need from both matrices
    int numRows1 = m_nRows;
    int numRows2 = matrix2.m_nRows;
    int numCols1 = m_nCols;
    int numCols2 = matrix2.m_nCols;

    // throw error if the matrices have different number of rows
    if (numRows1 != numRows2)
        throw std::invalid_argument("Attempt to join matrices with different numbers of rows is invalid!");

    // allocate memory for the result
    // note that only the number of columns increases
    T *newMatrixData = new T[numRows1 * (numCols1 + numCols2)];

    // copy the two matrices into the new one
    int linearIndex, resultLinearIndex;
    for (int i = 0; i < numRows1; i++)
    {
        for (int j = 0; j < (numCols1 + numCols2); j++)
        {
            resultLinearIndex = (i * (numCols1 + numcol2)) + j;

            // if j is in the left hand matrix, we get data from there
            if (j < numCols1)
            {
                linearIndex = (i * numCols1) + j;
                newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
            }
            // otherwise, j must be in the right hand matrix, so we get data from there
            else
            {
                linearIndex = (i * numCols2) + (j - numCols1);
                newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
            }
        }
    }

    // update the stored data
    m_nCols = numCols1 + numCols2;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = newMatrixData[i];

    return true;
}

template <class T>
int matrix2<T>::Sub2Ind(int row, int col)
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return (row * m_nCols) + col;
    else
        return -1;
}

// test if values are close enough to a fixed tolerance
template <class T>
bool matrix2<T>::CloseEnough(T f1, T f2)
{
    return fabs(f1 - f2) < 1e-9;
}

#endif