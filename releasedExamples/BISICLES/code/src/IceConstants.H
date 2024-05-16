#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _ICECONSTANTS_H_
#define _ICECONSTANTS_H_
#include "BaseNamespaceHeader.H"


#define GROUNDEDMASKVAL (1)
#define groundedmaskval (1)
#define FLOATINGMASKVAL (2)
#define floatingmaskval (2)
#define OPENSEAMASKVAL (4)
#define openseamaskval (4)
#define OPENLANDMASKVAL (8)
#define openlandmaskval (8)
#define OCEANMASKVAL OPENSEAMASKVAL
#define oceanmaskval openseamaskval

#if defined(CH_USE_FLOAT) || defined(CH_LANG_CC)

#define seconds_per_tropical_year     (31556926)
#define SECONDS_PER_TROPICAL_YEAR     (31556926)


#define HUGE_NORM (1.0e10)
#define TINY_NORM (1.0e-10)
#define tiny_thickness (1.0e-2)
#define TINY_THICKNESS (1.0e-2)
#define huge_thickness (1.0e+4)
#define HUGE_THICKNESS (1.0e+4)
#define HUGE_VEL (1.0e+6)
#define huge_vel (1.0e+6)
#define TINY_VEL (1.0e-3)
#define tiny_vel (1.0e-3)
#define ICE_DENSITY 918.0
#define SEA_WATER_DENSITY 1028.0
#define GRAVITY 9.81

#elif defined(CH_USE_DOUBLE) && !defined(CH_LANG_CC)

#define seconds_per_tropical_year     (31556926d0)
#define SECONDS_PER_TROPICAL_YEAR     (31556926d0)

#define HUGE_NORM (1.0d10)
#define TINY_NORM (1.0d-10)
#define tiny_thickness (1.0d-2)
#define TINY_THICKNESS (1.0d-2)
#define huge_thickness (1.0e+4)
#define HUGE_THICKNESS (1.0e+4)
#define HUGE_VEL (1.0d+6)
#define huge_vel (1.0d+6)
#define TINY_VEL (1.0d-2)
#define tiny_vel (1.0d-2)

#else

#error One of CH_USE_{FLOAT,DOUBLE} must be defined

#endif

#include "BaseNamespaceFooter.H"
#endif /*CH_ICECONSTANTS_H*/
