#ifndef QBMATRIX2_H
#define QBMATRIX2_H

template <class T>
class qbMatrix2
{
public:
    // the various constructors
    qbMatrix2();
    qbMatrix2(int nRows, int nCols);
    qbMatrix2(int nRows, int nCols, const T *inputData);
    qbMatrix2(const qbMatrix2<T> &inputMatrix);

    // destructor
    ~qbMatrix2();

    // configuration methods
    bool resize(int numRows, int numCols);

    // element access methods
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // overload the == operator
    bool operator==(const qbMatrix2<T> &rhs);

    // overload +, -, and * operators (friends)
    template <class U>
    friend qbMatrix2<U> operator+(const qbMatrix2<U> &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator+(const U &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator+(const qbMatrix2<U> &lhs, const U &rhs);

    template <class U>
    friend qbMatrix2<U> operator-(const qbMatrix2<U> &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator-(const U &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator-(const qbMatrix2<U> &lhs, const U &rhs);

    template <class U>
    friend qbMatrix2<U> operator*(const qbMatrix2<U> &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator*(const U &lhs, const qbMatrix2<U> &rhs);
    template <class U>
    friend qbMatrix2<U> operator*(const qbMatrix2<U> &lhs, const U &rhs);

private:
    int Sub2Ind(int row, int col);
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

/* ************************************************************************************************ */
// constructors and destructors

// the default constructor
template <class T>
qbMatrix2<T>::qbMatrix2()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

// construct an empty matrix (all elements 0)
template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols)
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
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const T *inputData)
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
qbMatrix2<T>::qbMatrix2(const qbMatrix2<T> &inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// destructor
template <class T>
qbMatrix2<T>::~qbMatrix2()
{
    if (m_matrixData != nullptr)
        delete[] m_matrixData;
}

/* ************************************************************************************************ */
// configuration functions

template <class T>
bool qbMatrix2<T>::resize(int numRows, int numCols)
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
T qbMatrix2<T>::GetElement(int row, int col)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
        return m_matrixData[linearIndex];
    else
        return -1;
}

template <class T>
bool qbMatrix2<T>::SetElement(int row, int col, T elementValue)
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
int qbMatrix2<T>::GetNumRows()
{
    return m_nRows;
}

template <class T>
int qbMatrix2<T>::GetNumCols()
{
    return m_nCols;
}

/* ************************************************************************************************ */
// overloaded operator functions

/* the + operator ********************************************************************************* */
// matrix + matrix
template <class T>
qbMatrix2<T> operator+(const qbMatrix2<T> &lhs, const qbMatrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRow * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar + matrix
template <class T>
qbMatrix2<T> operator+(const T &lhs, const qbMatrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs + rhs.m_matrixData[i];

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix + scalar
template <class T>
qbMatrix2<T> operator+(const qbMatrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] + rhs;

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

/* the - operator ********************************************************************************* */
// matrix - matrix
template <class T>
qbMatrix2<T> operator-(const qbMatrix2<T> &lhs, const qbMatrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRow * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar - matrix
template <class T>
qbMatrix2<T> operator-(const T &lhs, const qbMatrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs - rhs.m_matrixData[i];

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix - scalar
template <class T>
qbMatrix2<T> operator-(const qbMatrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] - rhs;

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

/* the * operator ********************************************************************************* */
// scalar * matrix
template <class T>
qbMatrix2<T> operator*(const T &lhs, const qbMatrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs * rhs.m_matrixData[i];

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * scalar
template <class T>
qbMatrix2<T> operator*(const qbMatrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = (numRows * numCols);
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
        tempResult[i] = lhs.m_matrixData[i] * rhs;

    qbMatrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

#endif