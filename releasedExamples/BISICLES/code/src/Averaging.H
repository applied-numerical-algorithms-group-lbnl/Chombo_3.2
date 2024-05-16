#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AVERAGING_H_
#define _AVERAGING_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

/** A set of utility functions for averaging data from a fine grid to
    a coarse grid
*/

void horizontalAverage(LevelData<FArrayBox>& a_crseData,
                       const LevelData<FArrayBox>& a_fineData,
                       int a_nRef);

void horizontalAverageFace(LevelData<FluxBox>& a_crseData,
                           const LevelData<FluxBox>& a_fineData,
                           int a_nRef);

void averageAllDim(LevelData<FArrayBox>& a_crseData, 
                   const LevelData<FArrayBox>& a_fineData,
                   int a_nRef);


void averageAllDimFace(LevelData<FluxBox>& a_crseData, 
                       const LevelData<FluxBox>& a_fineData,
                       int a_nRef);

/// single-box versions

void horizontalAverage(FArrayBox& a_crseData,
                       const FArrayBox& a_fineData,
                       const Box& a_crseBox,
                       int a_nRef);

void horizontalAverageFace(FluxBox& a_crseData,
                           const FluxBox& a_fineData,
                           const Box& a_crseCellBox,
                           int a_nRef);

void averageAllDim(FArrayBox& a_crseData, 
                   const FArrayBox& a_fineData,
                   const Box& a_crseBox,
                   int a_nRef);


void averageAllDimFace(FluxBox& a_crseData, 
                       const FluxBox& a_fineData,
                       const Box& a_crseBox,
                       int a_nRef);


#include "NamespaceFooter.H"

#endif
