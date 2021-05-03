#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <fstream>
#include <list>


class Matrix
{
    /*
   wyglad tablic:
      [x][y]
              y[0] y[1]
                |   |
       x[0] - |0,0|0,1|
       x[1] - |1,0|1,1|
       x[2] - |2,0|2,1|
   */

private:
    double** _matrix;
    int _xSize, _ySize;

    double** InstanciateMatrix(int, int);

public:
    Matrix();
    Matrix(int x, int y);
    Matrix(int x, int y, double z);
    Matrix(int);
    Matrix(std::string path);
    ~Matrix();

    Matrix(const Matrix& m1);

    void store(std::string filename, std::string path) const;

    void operator = (const Matrix&);
    Matrix operator+(const Matrix& m2) const;
    Matrix operator-(const Matrix& m2) const;
    Matrix operator*(const Matrix& m2) const;


    Matrix addReturnCopy(const Matrix&) const;
    Matrix multiReturnCopy(const Matrix&) const;
    Matrix subbReturnCopy(const Matrix&) const;

    Matrix* multiReturnPointer(const Matrix&) const;
    Matrix* addReturnPointer(const Matrix&) const;
    Matrix* subbReturnPointer(const Matrix&) const;

    //Matrix* duplicate() const;
    void set(int y, int x, double val);
    double get(int x, int y) const;
    void PrintMatrix()const;



    int cols() const;
    int rows() const;

};

