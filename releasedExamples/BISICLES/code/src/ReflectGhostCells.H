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
//  ReflectGhostCells.H
// ============
//
// Simple function to fill domain ghost cells with reflected values
//

#ifndef _REFLECTGHOSTCELLS_H_
#define _REFLECTGHOSTCELLS_H_

#include "ProblemDomain.H"
#include "LevelData.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

void 
ReflectGhostCells(LevelData<FArrayBox>& a_phi,
		  const ProblemDomain& a_domain,
		  const int a_dir,
		  const Side::LoHiSide a_side);

void 
ReflectGhostCells(FArrayBox& a_phi,
		  const ProblemDomain& a_domain,
		  const IntVect& a_ghostVect,
		  const int dir, 
		  const Side::LoHiSide a_side);

#include "NamespaceFooter.H"

#endif
