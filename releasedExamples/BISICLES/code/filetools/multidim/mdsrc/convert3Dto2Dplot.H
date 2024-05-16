#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONVERT3DTO2DPLOT_H_
#define _CONVERT3DTO2DPLOT_H_


//#define SAVE_SPACEDIM CH_SPACEDIM
#undef CH_SPACEDIM
#define CH_SPACEDIM  2
#include <string>
#include "IntVect.H"
//#undef CH_SPACEDIM
//#define CH_SPACEDIM SAVE_SPACEDIM
//#define CH_SPACEDIM 2
//#undef SAVE_SPACEDIM

//#include "UsingNamespace.H"

using namespace Chombo;

using std::string;

int convert3DNCTo2DPlot(const string& fname1d,
                        const string& fname2d);



#endif
