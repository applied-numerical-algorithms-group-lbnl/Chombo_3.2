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
#include "ProblemDomain.H"



#ifdef USE_PROTO


///
/**
   These functions are meant to interface with the Proto language.
*/
namespace ProtoCh
{
  typedef Proto::Point           pr_pt;
  typedef Chombo::IntVect        ch_iv;
  typedef Chombo::Box            ch_box;
  typedef Proto::Box             pr_box;
  typedef Proto::ProblemDomain   pr_dom;
  typedef Chombo::ProblemDomain  ch_dom;
  ///get point from intvect
  extern pr_pt  getPoint( const ch_iv& a_iv)
  {
    pr_pt retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_iv[idir];
    }
    return retval;
  }
  
  ///get intvect from point
  inline ch_iv  getIntVect(const  pr_pt  & a_pt)
  {
    ch_iv retval;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval[idir] = a_pt[idir];
    }
    return retval;
  }


  ///get chombo box from proto box
  inline ch_box getBox(const pr_box & a_bx)
  {
    ch_iv ivlo = getIntVect(a_bx.low());
    ch_iv ivhi = getIntVect(a_bx.high());
    return ch_box(ivlo, ivhi);
  }
  
  /// gets chombo problem domain from proto problem domain
  inline ch_dom getChomboDomain(const pr_dom& a_dom)
  {
    pr_box                box_pr = a_dom.box();
    std::array<bool, DIM> per_pr = a_dom.periodicity();
    ch_box                box_ch = getBox(box_pr);
    bool                  per_ch[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      per_ch[idir] = per_pr[idir];
    }
    return ch_dom(box_ch, per_ch);
  }
  
  /// more consistent with other names here
  inline ch_dom getDomain(const pr_dom& a_dom)
  {
    return getChomboDomain(a_dom);
  }
  
  /// gets proto box from chombo box
  inline pr_box getProtoBox(const ch_box& a_box)
  {
    pr_pt ptlo = getPoint(a_box.smallEnd());
    pr_pt pthi = getPoint(a_box.bigEnd());
    return pr_box(ptlo, pthi);
  }
  
  /// gets proto problem domain from Chombo problem domain
  inline pr_dom getProtoDomain(const ch_dom& a_dom)
  {
    ch_box                box_ch = a_dom.domainBox();
    pr_box                box_pr = getProtoBox(box_ch);
    const bool*           per_ch = a_dom.isPeriodicVect(); //bool[DIM]
    std::array<bool, DIM> per_pr;
    for(int idir = 0; idir < DIM; idir++)
    {
      per_pr[idir] = per_ch[idir];
    }
    return pr_dom(box_pr, per_pr);
  }




} 

#endif

#endif
