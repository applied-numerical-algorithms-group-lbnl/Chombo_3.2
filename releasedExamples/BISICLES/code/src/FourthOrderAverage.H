#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERAVERAGE_H_
#define _FOURTHORDERAVERAGE_H_

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "MayDay.H"

#include "NamespaceHeader.H"

/// utility functions for 4th-order averaging

///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
*/
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        int a_sgn = 1);


///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
     using one-sided corrections at domain boundaries
*/
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        const ProblemDomain&  a_domain,
                        int a_sgn = 1);


///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
*/
void fourthOrderAverageCell(FArrayBox& a_phi,
                            int a_sgn = 1);

///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
     using one-sided corrections at domain boundaries
*/
void fourthOrderAverageCell(FArrayBox& a_phi,
                            const ProblemDomain&  a_domain,
                            const Box&            a_bx,
                            int a_sgn = 1);

#include "NamespaceFooter.H"

#endif
