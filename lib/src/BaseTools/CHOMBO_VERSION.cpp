#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CHOMBO_VERSION.H"

#include "BaseNamespaceHeader.H"

//
// The definition of our version string.
//
#define bl_str(s) # s
#define bl_xstr(s) bl_str(s)

void getChomboVersion(string& a_version,
                      string& a_flags,
                      string& a_date,
                      string& a_time)
{
  a_version = bl_xstr(CHOMBO_VERSION);
  a_flags   = bl_xstr(CHOMBO_FLAGS);
  a_date    = __DATE__;
  a_time    = __TIME__;
}

#undef bl_str
#undef bl_xstr

#include "BaseNamespaceFooter.H"
