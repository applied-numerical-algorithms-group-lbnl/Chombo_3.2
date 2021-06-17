#include "MatrixToCSC.H"
// #define CMP_COL_ROW(r,c,i,j) (c[_x]!=c[_p] ? c[i]-c[j] : r[i]-r[j])

void SparseMatrix::addElement(const unsigned int r, const unsigned int c,
  const double v) {

  // we assume that the current format is COO, we cannot insert in CSC
  assert(format == FORMAT_COO);

  /* elements[nnz].row   = r;
  elements[nnz].col   = c;
  elements[nnz].value = v; */
  MatrixElement e;

  e.row   = r;
  e.col   = c;
  e.value = v;

  elements.push_back(e);

  // update the nRows in the SparseMatrix object
  if(nRows <= r)
    nRows = r+1; // we have to add 1 because we start from 0

  // same logic for nCols
  if(nCols <= c)
    nCols = c+1;

  nnz++;
}

void SparseMatrix::deleteElement(int r, int c) {
  // modify only if we are in COO
  assert(format == FORMAT_COO);

  // std::vector<int>::const_iterator itv;
  auto itv = elements.begin();

  while(itv != elements.end() && !(itv->row == r && itv->col == c) ) {
    itv++;
  }

  if(itv != elements.end()) { // if we have found the element, delete it!
    elements.erase(itv);
    updateNRowsCols();
    nnz--;
  }
}

void SparseMatrix::removeLastElement() {
  // modify only if we are in COO
  assert(format == FORMAT_COO);

  if(!elements.empty()) { // we can only delete if there is an element...
    elements.pop_back();
    nnz--;
    updateNRowsCols();
  }
}

// Based on the elements of the object of SparseMatrix, it updates the
// matrix dimensions
void SparseMatrix::updateNRowsCols() {
  nRows = 0, nCols = 0;

  for (auto& it : elements) {
    if(it.row+1 > nRows)
      nRows = it.row+1;
    if(it.col+1 > nCols)
      nCols = it.col+1;
  }
}

// we sort first by col and then by row
void SparseMatrix::sortElements() {
  // we cannot sort this array if the format is CSC
  assert(format == FORMAT_COO);

  /* ideally, we want nnz == maxElements (i.e., all nnz elements specified
  before sorting the SparseMatrix arrays) */
  /* if(nnz < maxElements) {
    std::cerr << "WARNING: not all nonzero elements were specified"
      "before sorting the SparseMatrix arrays" << '\n';
  } */

  // using qsort from glibc
  // qsort(elements, nnz, sizeof(MatrixElement), elementsCmp);

  sort(elements.begin(), elements.end(), &elementsCmp);
}

// function used by qsort to compare between two elements, first by col and
// then by row
bool elementsCmp(const MatrixElement &e1, const MatrixElement &e2) {
  if ( e1.col != e2.col )
    return e1.col < e2.col;
  else
    return e1.row  < e2.row;
}

void SparseMatrix::transformToCSC() {
  // we can only transform to CSC format if we are in COO
  assert(format == FORMAT_COO);

  int i, t1, t, p, k;

  // first thing we need to do is to sort the elements, first by col and
  // then by row
  sortElements();

  k = 0;
  /*  Remove duplicates */
  for(i = 1; i < nnz; i++) {
      if( (elements.at(k).row != elements.at(i).row) ||
          (elements.at(k).col != elements.at(i).col) ) {
        k++;

        elements.at(k).row   = elements.at(i).row;
        elements.at(k).col   = elements.at(i).col;
        elements.at(k).value = elements.at(i).value;
      } else {
        // if the rows and columns are the same between two elements,
        // add the values.
        elements.at(k).value += elements.at(i).value;
      }
  }

  nnz = k+1; // this is the nnz values without duplicates
  elements.at(0).col = 0; // the first value of the CSC col array is always 0
  t = 0;
  p = 1;

  /*  Build CSC column pointers in the col array */
  for(i = 1; i < nnz; i++) {
      t1 = elements.at(i).col;
      if( elements.at(i).col != t )
          elements.at(p++).col = i;
      t = t1;
  }

  elements.at(p).col = i; // this marks the last value of the compressed
                       // column vector Matrix->Col (i == nnz)

  format = FORMAT_CSC; // updating the format to CSC.
}

unsigned int SparseMatrix::getNNZElements() {
  return(nnz);
}

//DELETE THIS ASAP
#include<stdio.h>
void SparseMatrix::printArrays() {
  printf("Printing in ");

  if(format == FORMAT_COO)
    printf("COO format:\n\n");
  else
    printf("CSC format:\n\n");

  int i;

  printf("Row:\n");

  for(i = 0; i < nnz; i++) {
      printf("%d\n", elements.at(i).row);
  }

  printf("Col:\n");
  for(i = 0; i < nnz; i++) {
      printf("%d\n", elements.at(i).col);
  }

  printf("Elements:\n");
  for(i = 0; i < nnz; i++) {
      printf("%f\n", elements.at(i).value);
  }
}

void SparseMatrix::extractSuperLUData(int *rows, int *colPtrs, double *values) {

  if(format == FORMAT_COO)
    transformToCSC();

  // dumping data into the arrays
  for( int i = 0; i < nnz; i++ ) {
    values[i]  = elements.at(i).value;
    rows[i]    = elements.at(i).row;
    colPtrs[i] = elements.at(i).col;
  }
}

unsigned int SparseMatrix::getNRows() {
  return nRows;
}

unsigned int SparseMatrix::getNCols() {
  return nCols;
}

// void SparseMatrix::transformToCOO(); not needed for SuperLU
