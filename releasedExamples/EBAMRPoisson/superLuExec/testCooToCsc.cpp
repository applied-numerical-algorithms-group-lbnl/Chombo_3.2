#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Code adapted from the SuperLU example: superlu.c
#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "memusage.H"
#include "memtrack.H"

#ifdef CH_USE_SUPERLU
#include "slu_ddefs.h"
#include <stdio.h>
#endif // CH_USE_SUPERLU

#include "MatrixToCSC.H"

int main(int argc, char* argv[])
{
#ifdef CH_USE_SUPERLU
  SuperMatrix A, L, U, B;
  double   *rhs, *sol;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      nrhs, info, i;
  superlu_options_t options;
  SuperLUStat_t stat;
#endif // CH_USE_SUPERLU

  SparseMatrix sm;
  int *rows, *cols, nnz, m, n;
  double *vals;

  sm.addElement(1,1,21);
  sm.addElement(3,3,5);
  sm.addElement(2,1,12);
  sm.addElement(3,4,21);
  sm.addElement(1,0,12);
  sm.addElement(0,2,21);
  sm.addElement(4,4,18);
  sm.addElement(4,1,12);
  sm.addElement(0,0,19);
  sm.addElement(2,2,16);
  sm.addElement(4,0,12);
  sm.addElement(0,3,21);

  /* With the following instructions we can verify the correct
  behaviour of the SparseMatrix methods.

  cout << sm.getNRows() << " " << sm.getNCols() << endl;
  cout << sm.getNNZElements() << endl;
  sm.deleteElement(4,4);
  sm.deleteElement(3,4);
  sm.removeLastElement();
  cout << sm.getNNZElements() << endl;
  cout << sm.getNRows() << " " << sm.getNCols() << endl;
  */

  /* For verification purposes, we can print the arrays
     before and after transforming to CSC. */

  // sm.printArrays();

  sm.transformToCSC();

  // sm.printArrays();

  m = sm.getNRows();
  n = sm.getNCols();

  cout << "SM(" << m << "," << n << ")" << endl;

  nnz = sm.getNNZElements();

  // allocate space for the SuperLU arrays
  vals = (double *) malloc(nnz * sizeof(double));
  rows = (int    *) malloc(nnz * sizeof(int   ));
  cols = (int    *) malloc(nnz * sizeof(int   ));

  sm.extractSuperLUData(rows, cols, vals);

  cout << "Printing data in SuperLU format..." << endl;

  for (int i = 0; i < nnz; i++) {
    cout << rows[i] << ", " << cols[i] << ": " << vals[i] <<endl;
  }

#ifdef CH_USE_SUPERLU
  // Create matrix A in the format expected by SuperLU.
  dCreate_CompCol_Matrix(&A, m, n, nnz, vals, rows, cols,
    SLU_NC, SLU_D, SLU_GE);

  nrhs = 1;

  if ( !(rhs = doubleMalloc(m * nrhs)) )
    ABORT("Malloc fails for rhs[].");

  for (i = 0; i < m; ++i)
    rhs[i] = 1.0;

  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
  if ( !(perm_r = intMalloc(m)) )
    ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) )
    ABORT("Malloc fails for perm_c[].");

  // Set the default input options.
  set_default_options(&options);
  options.ColPerm = NATURAL;

  // Initialize the statistics variables.
  StatInit(&stat);

  // Solve the linear system.
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  sol = (double*) ((DNformat*) B.Store)->nzval;

  printf("dgssv() returns INFO = %d\n", info);

  printf("\nOk, let's print the solution...\n");

  for (i = 0; i < m; i++)
    printf("%.8E\n", sol[i]);

  dPrint_CompCol_Matrix((char*) "A", &A);

  /* We can also print the factoring matrices L, U: */
  dPrint_CompCol_Matrix((char*) "U", &U);
  dPrint_SuperNode_Matrix((char*) "L", &L);
  print_int_vec((char *)"\nperm_r", m, perm_r);



  // De-allocate storage
  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
#endif // CH_USE_SUPERLU

  return(0);
}
