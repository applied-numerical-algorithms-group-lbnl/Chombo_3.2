#ifndef _MATRIXTOCSC_H_
#define _MATRIXTOCSC_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>

using std::cerr;
using std::endl;

#define FORMAT_CSC (1)
#define FORMAT_COO (0)
// #define PTR_TO_ELEMENT (const MatrixElement *)

typedef struct {
  unsigned int row, col;
  double value;
} MatrixElement;

// function used for stl vector sort
bool elementsCmp(const MatrixElement &e1, const MatrixElement &e2);

class SparseMatrix {
  private:
    std::vector<MatrixElement> elements;
    unsigned int nnz, format, nRows, nCols;
    void updateNRowsCols();

  public:
    void addElement(const unsigned int r, const unsigned int c, const double v);
    void deleteElement(int r, int c);
    void removeLastElement();
    void sortElements();
    void transformToCSC();
    // void transformToCOO(); not needed for SuperLU
    void printArrays();
    unsigned int getNNZElements();
    unsigned int getNRows();
    unsigned int getNCols();
    void extractSuperLUData(int *rows, int *colPtrs, double *values);


    // Constructor
    SparseMatrix() {
      nnz = 0;
      nRows = 0;
      nCols = 0;
      format = FORMAT_COO;
    }
};

#endif // _MATRIXTOCSC_H_
