#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ProtoInterface.H"
#include "SPACE.H"
#include "UsingNamespace.H"

#ifdef USE_PROTO

///get point from intvect
Point
ProtoCh::getPoint( const CH_XD::IntVect& a_iv)
{
  Point retval;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    retval[idir] = a_iv[idir];
  }
  return retval;

}

/// gets proto box from chombo box
Proto::Box  
ProtoCh::getProtoBox( const CH_XD::Box& a_box)
{
  Point ptlo = getPoint(a_box.smallEnd());
  Point pthi = getPoint(a_box.bigEnd());
  return Proto::Box(ptlo, pthi);
}

///get intvect from point
CH_XD::IntVect 
ProtoCh::getIntVect(const  Point  & a_pt)
{
  CH_XD::IntVect retval;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    retval[idir] = a_pt[idir];
  }
  return retval;
}

///get chombo box from proto box
CH_XD::Box 
ProtoCh::getBox(const Proto::Box & a_bx)
{
  CH_XD::IntVect ivlo = getIntVect(a_bx.low());
  CH_XD::IntVect ivhi = getIntVect(a_bx.high());
  return CH_XD::Box(ivlo, ivhi);
}



#endif


