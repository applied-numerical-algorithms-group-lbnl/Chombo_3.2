#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _LEVELMAPPEDDERIVATIVES_H_
#define _LEVELMAPPEDDERIVATIVES_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelSigmaCS.H"

// for now, include "MappedDerivatives.H for derivComponent
//#include "MappedDerivatives.H"

#include "NamespaceHeader.H"

/** A set of utility functions for computing derivatives in the mapped
    sigma coordinate system used by ice-sheet models.
    if we represent (x,y,z) as (sigma, xtilde, ytilde) in 3d, where
    xtilde = x, ytilde = y, and sigma = (z_surface - z)/H, where H =
    ice-sheet thickness, then we have the following geometric
    relations: 

d/(dx) = d/(d xtilde) + (1/H)Delta_x*(d/d sigma)
d/(dy) = d/(d ytilde) + (1/H)Delta_y*(d/d sigma)
d/(d sigma) = -(1/H)d/(d sigma)

where Delta_x and Delta_y are geometric factors computed in the SigmaCS.

These utility functions implement this differencing so that it's
convenient to use. 

*/

/// convenient component lookup -- componentDir*SpaceDim + derivDir
int inline derivComponent(int derivDir, int componentDir)
{return componentDir*SpaceDim + derivDir;}




/// compute cell-centered derivatives
/** it's assumed that any ghost-cell boundary conditions have already
    been set 
    -- derivGhost is the number of ghost cells in ccDeriv in which to 
       compute the derivative. A check is made to ensure that ccData 
       has enough ghost cells to make this possible. Note also that 
       derivGhost may be negative.
*/
void 
computeCCDerivatives(LevelData<FArrayBox>& a_ccDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_derivGhost = IntVect::Zero,
		     const bool a_maskOneSide = false);

/// compute face-centered derivatives of cell-centered data on all faces
void
computeFCDerivatives(LevelData<FluxBox>& a_fcDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_derivGhost = IntVect::Zero);



/// compute face-centered derivatives of cell-centered data
/** This version takes cell-centered derivatives in the mapped
    coordinates and averages them to faces and applies the correction
    required due to the mapping. This is essentially designed to
    support the horizontal ice-velocity elliptic operator.

    fcDeriv should be a face-centered FArrayBox

    ccData is the cell-centered data (needed on ghost cells normal 
           to the face in faceDir)

    ccDeriv contains the cell-centered derivatives used to compute
            tangential derivatives. These should be indexed in the
            same way as the cell-centered derivatives returned by the
            computeCCDerivatives function, and need to be defined on  
            one ghost cell in the faceDir direction.

    faceDerivBox should be a face-centered box corresponding to the 
                 region in fcDeriv where derivatives are needed.     
     
*/
void
computeFCDerivatives(LevelData<FluxBox>& a_fcDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelData<FArrayBox>& a_ccDeriv,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_ghostVect = IntVect::Zero);


#include "NamespaceFooter.H"

#endif
