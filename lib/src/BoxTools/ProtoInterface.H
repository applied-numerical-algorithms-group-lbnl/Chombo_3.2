#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PROTOINTERFACE_H_
#define _PROTOINTERFACE_H_

#ifdef USE_PROTO
#include "Proto.H"
#endif
#include "IntVect.H"
#include "Box.H"
#include "BaseFab.H"



#ifdef USE_PROTO


///
/**
   These functions are meant to interface with the Proto language.
*/
namespace ProtoCh
{
  typedef Proto::Point    pr_pt;
  typedef Chombo::IntVect ch_iv;
  typedef Chombo::Box     ch_box;
  typedef Proto::Box      pr_box;

  ///get point from intvect
  extern pr_pt  getPoint( const ch_iv& a_iv)
  {
    pr_pt retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = a_iv[idir];
    }
    return retval;

  }

  /// gets proto box from chombo box
  inline pr_box getProtoBox(const ch_box& a_box);
  {
    pr_pt ptlo = getPoint(a_box.smallEnd());
    pr_pt pthi = getPoint(a_box.bigEnd());
    return pr_box(ptlo, pthi);
  }

  ///get intvect from point
  inline ch_iv  getIntVect(const  pr_pt  & a_pt);
  {
    ch_iv retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = a_pt[idir];
    }
    return retval;
  }

  ///get chombo box from proto box
  inline Box getBox(const pr_ & a_bx);
  {
    IntVect ivlo = getIntVect(a_bx.low());
    IntVect ivhi = getIntVect(a_bx.high());
    return Box(ivlo, ivhi);
  }

} 

#endif

#endif
