#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  FillFromReference.H
// ============
//
// function to 
//

#ifndef _FILLFROMREFERENCE_H_
#define _FILLFROMREFERENCE_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "CoarseAverage.H"
#include "RealVect.H"


#include "NamespaceHeader.H"

/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference FArrayBox
/**
   
 */

void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const FArrayBox& a_srcData,
		       const RealVect&  a_destDx,
		       const RealVect&  a_srcDx,
		       const IntVect& a_srcGhost,
		       bool a_verbose, 
		       CoarseAverage::averageType a_avgType = CoarseAverage::arithmetic);



/// function to fill a single level's worth of state data 
/// (in a LevelData<FArrayBox>) from a reference  LevelData<FArrayBox>
/**
   
 */
void FillFromReference(LevelData<FArrayBox>& a_destData,
		       const LevelData<FArrayBox>&  a_srcData,
		       const RealVect&  a_destDx,
		       const RealVect&  a_srcDx,
		       bool a_verbose, 
		       CoarseAverage::averageType a_avgType = CoarseAverage::arithmetic);


// fill a single LevelData wth a flattened AMR hierarchy's worth of data
/** This version does not perserve const-ness of the src data -- data in 
    covered regions is replaced by averaged-down data.
 */
void
flattenCellData(LevelData<FArrayBox>& a_destData,
                const RealVect& a_destDx,
                Vector<LevelData<FArrayBox>* >& a_srcData,
                const Vector<RealVect>& a_srcDx,
                bool a_verbose, 
		CoarseAverage::averageType a_avgType = CoarseAverage::arithmetic);


// fill a single LevelData wth a flattened AMR hierarchy's worth of data
/** This version does not perserve const-ness of the src data -- data in 
    covered regions is replaced by averaged-down data. Version to use 
    refcountedPtrs
 */
void
flattenCellData(LevelData<FArrayBox>& a_destData,
                const RealVect& a_destDx,
                Vector<RefCountedPtr<LevelData<FArrayBox> >  >& a_srcData,
                const Vector<RealVect>& a_srcDx,
                bool a_verbose, 
		CoarseAverage::averageType a_avgType = CoarseAverage::arithmetic);


// fill a single LevelData wth a flattened AMR hierarchy's worth of data
/** This version perserves const-ness of the src data, but requires some 
    constraints on the fine meshes (grids on all levels must be coarsenable
    to the destination resolution.
 */
void
flattenCellDataConst(LevelData<FArrayBox>& a_destData,
                     const RealVect& a_destDx,
                     const Vector<LevelData<FArrayBox>* >& a_srcData,
                     const Vector<RealVect>& a_srcDx,
                     bool a_verbose, 
		     CoarseAverage::averageType a_avgType = CoarseAverage::arithmetic);

#include "NamespaceFooter.H"

#endif
