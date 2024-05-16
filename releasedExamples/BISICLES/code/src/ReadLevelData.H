#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// ReadLevelData.H
//===========================================================================

#ifndef _READLEVELDATA_H_
#define _READLEVELDATA_H_
#include "AMRIO.H"
#include "NamespaceHeader.H"

//fill a LevelData<FArrayBox> from data stored in an AMR file 
void readLevelData(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_data,
		   Real& a_dx,
		   const std::string a_file,
		   const Vector<std::string>& a_names,
		   int a_nComp);


///create one or more AMR Hieracrhies from data stored in an AMR file
/**
   
   \param Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > a_data : each 
   Vector<RefCountedPtr<LevelData<FArrayBox> > > is an AMR hierarchy.
   \param a_crseDx : mesh spacing on level 0
   \param a_refRatio : refinement ratio between meshes.
   \param a_file : file name to read from
   \param a_name : for each name in a_names look for either 
   (a) a single component LevelData named "name"
   (b) a multi-component LevelData, stored consecutively and beginning with name
   \param a_nComp : the number of components in each heirachy
*/
void readMultiLevelData
(Vector< Vector <RefCountedPtr<LevelData<FArrayBox> > > >& a_data,
 Real& a_crseDx,
 Vector<int>& a_refRatio,
 const std::string a_file,
 const Vector<std::string>& a_names, 
 int a_nComp);

#include "NamespaceFooter.H"

#endif
