#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BISICLES_VERSION_H_
#define _BISICLES_VERSION_H_
#include "CHOMBO_VERSION.H"
#include "BaseNamespaceHeader.H"
/**
   BISICLES VERSION NUMBER
 */
#define BISICLES_VERSION_MAJOR 1
#define BISICLES_VERSION_MINOR 1
#define BISICLES_PATCH_NUMBER 0
/**
   The Chombo version that this BISICLES version is intended for
 */
#define BISICLES_CHOMBO_VERSION 3.2.5
#define BISICLES_CHOMBO_SVN_REVISION 23113

/** signal that this is the svn trunk version of BISICLES, rather than 
    the release version. This should be changed in the release version 
    of this file */
#define BISICLES_TRUNK TRUE

/**
   \mainpage BISICLES: An ice flow model based on the Chombo 
   Block-Structured Adaptive Mesh Refinement Library

   <a href="http://seesar.lbl.gov/anag/bisicles">BISICLES webpage</a>

*/


#include "BaseNamespaceFooter.H"
#endif
