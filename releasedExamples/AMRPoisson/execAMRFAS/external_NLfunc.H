#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTERNALNLFUNC_H_
#define _EXTERNALNLFUNC_H_

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelDataOps.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

void setNL_Level(LevelData<FArrayBox>&  a_NL,
                 LevelData<FArrayBox>&  a_dNL,
                 const LevelData<FArrayBox>&  a_u);

void setNL_piece(Vector<LevelData<FArrayBox>* > a_NL,
                 Vector<LevelData<FArrayBox>* > a_dNL,
                 Vector<LevelData<FArrayBox>* > a_u,
                 int a_finestLevel);

#include "NamespaceFooter.H"
#endif
