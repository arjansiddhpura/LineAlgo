#ifndef QBMATRIX2_H
#define QBMATRIX2_H

template <class T>
class qbMatrix2
{
    public:
    // the various constructors
    qbMatrix2();
    qbMatrix2(int nRows, int nCols);
    qbMatrix2(int nRows, int nCols, const T *inputData);    // pointer to a linear array of input data
    qbMatrix2(const qbMatrix2<T>& inputMatrix);             // copy constructor using another instance

    // destructor
    ~qbMatrix2();

    // configuration methods
    
};

#endif