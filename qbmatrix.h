#ifndef QBMATRIX2_H
#define QBMATRIX2_H

template <class T>
class qbMatrix2
{
public:
    // the various constructors
    qbMatrix2();
    qbMatrix2(int nRows, int nCols);
    qbMatrix2(int nRows, int nCols, const T *inputData); // pointer to a linear array of input data
    qbMatrix2(const qbMatrix2<T> &inputMatrix);          // copy constructor using another instance

    // destructor
    ~qbMatrix2();

    // configuration methods
    bool resize(int numRows, int numCols); // allows changing the size of the matrix

    // element access methods
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // overload the == operator
    bool operator==(const qbMatrix2<T> &rhs); // test equality of two matrices

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

#endif