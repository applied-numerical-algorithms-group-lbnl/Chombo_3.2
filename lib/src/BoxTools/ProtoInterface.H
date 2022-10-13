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


using Proto::Point;
using Proto::BoxData;
using Proto::Var;
///these functions are meant to interface with the Proto language.
namespace ProtoCh
{


///get point from intvect
  extern Point      getPoint( const CH_XD::IntVect& a_iv);

/// gets proto box from chombo box
  extern Proto::Box  getProtoBox( const CH_XD::Box& a_box);

///get intvect from point
  extern CH_XD::IntVect getIntVect(const  Point  & a_pt);

///get chombo box from proto box
  extern CH_XD::Box getBox(const Proto::Box & a_bx);

///
  template <typename T, int ncomp>
  void  aliasBoxData(BoxData<T, ncomp> & a_boxdata, const CH_XD::BaseFab<T>& a_bfr, const int a_startcomp = 0)
    {
      CH_assert( a_bfr.nComp() == (ncomp + a_startcomp));
      Proto::Box region = ProtoCh::getProtoBox(a_bfr.box());
      a_boxdata.define(a_bfr.dataPtr(a_startcomp), region, ncomp);
    }


} 

#endif

#endif
