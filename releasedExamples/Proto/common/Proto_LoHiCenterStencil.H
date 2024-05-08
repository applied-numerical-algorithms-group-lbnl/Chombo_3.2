
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Vector.H"

#ifndef __Lo_Hi_Center_Stencil__
#define __Lo_Hi_Center_Stencil__
/// 
/**
   Tools to manage center/lo/hi stencil paradigm everywhere in the vector operators.
   For now, it seems simplest to just leave everything just hanging out there as 
   static functions and public data.
   This will at least standardize apply and computing the one-sided boxes.
**/
namespace PrCh_Lo_Hi_Center
{
  typedef Proto::Stencil<double>        pr_sten;
  typedef Proto::Point                  pr_pt;
  typedef Proto::Box                    pr_box;
  typedef Proto::BoxData<double, 1  >   pr_box_data_sca;
  typedef Proto::Shift                  pr_shift;
  typedef Proto::ProblemDomain          pr_dom;
  typedef Chombo::MayDay                ch_mayday;
  ///
  class LHC_MetaData
  {
  public:
    pr_sten m_lo_sten;
    pr_sten m_hi_sten;
    pr_sten m_ce_sten;
    int     m_has_lo ;
    int     m_has_hi ;
    pr_box  m_lo_box ;
    pr_box  m_hi_box ;
    pr_box  m_ce_box ;

    /// 
    void setBoxes(int     a_has_lo ,
                  int     a_has_hi ,
                  pr_box  a_lo_box ,
                  pr_box  a_hi_box ,
                  pr_box  a_ce_box )
    {
      m_has_lo = a_has_lo;
      m_has_hi = a_has_hi;
      m_lo_box = a_lo_box;
      m_hi_box = a_hi_box;
      m_ce_box = a_ce_box;
    }
    
  };

  /// compute boxes for face-centered destination
  int
  loHiCenterFace(pr_box        &  a_lo_box,
                 pr_box        &  a_hi_box,
                 pr_box        &  a_ce_box,
                 int           &  a_has_lo,
                 int           &  a_has_hi,
                 const pr_box  &  a_flux_box, //called inBox in VTO.cpp
                 const pr_box  &  a_valid,
                 const int     &  a_face_dir,
                 const pr_dom  &  a_domain)
  {
    CH_TIME("loHiCenterFace");
    /**
       This is writteen following ViscousTensorOp::loHiCenterFace.
       The box calculus here is a bit fussy so I am trying
       to keep this as close to the original as possible.
       Proto Box syntax is different in places, however.  They
       like to return stuff instead of change in place.
    **/

    //default to that blessed case where we are not near the boundary
    a_ce_box = a_flux_box;
    a_has_lo = 0;
    a_has_hi = 0;
    ///I am keeping this comment as written because the date is hilarious (dtg 12/21/2023)
    /**
     // if we're periodic in a_dir, then all boxes are
     // contained in the domain (DFM 2/5/10)
     **/
    bool periodic = a_domain.isPeriodic(a_face_dir);
    if(!periodic)
    {
      pr_box domain_face_box = a_domain.box();
      domain_face_box.m_high[a_face_dir] += 1;
      domain_face_box.recomputeSize();  // required to keep box data consistent after above shenanigans.

      //have to catch return values because proto has a const semantic for grow
      a_ce_box =  a_flux_box;
      a_ce_box =  a_ce_box.grow(a_face_dir, 1);
      a_ce_box &= domain_face_box;
      a_ce_box =  a_ce_box.grow(a_face_dir,-1);

      //flux box intersected with domain
      pr_box in_box = a_flux_box;
      in_box  &= domain_face_box;
        
      //Shift/intersect/shift back for high and low sides
      /// to see if we are at the domain boundary
      //Have to catch return values because shift is constant in proto
      pr_box tmp_hi = in_box;
      pr_box tmp_lo = in_box;
      tmp_hi = tmp_hi.shift(a_face_dir, 1); 
      tmp_lo = tmp_lo.shift(a_face_dir,-1);
        
      tmp_hi &= domain_face_box;
      tmp_lo &= domain_face_box;
        
      tmp_hi = tmp_hi.shift(a_face_dir,-1);
      tmp_lo = tmp_lo.shift(a_face_dir, 1);
        
      if(tmp_hi != in_box)
      {
        a_has_hi = 1;
        a_hi_box = pr_box(in_box.m_low, in_box.m_high);
        a_hi_box.m_low[a_face_dir] = in_box.m_high[a_face_dir]; //single plane in diffdir
        a_hi_box.recomputeSize();//required to keep box data consistent after above shenanigans.
      }
      if(tmp_lo != in_box)
      {
        a_has_lo = 1;
        a_lo_box = pr_box(in_box.m_low, in_box.m_high);
        a_lo_box.m_high[a_face_dir] = in_box.m_low[a_face_dir]; //single plane in diffdir
        a_lo_box.recomputeSize();//required to keep box data consistent after above shenanigans.
      }
    }
    return 0;
  }
  ///
  /**
     Apply the composite stencil on the appropriate boxes.
     1. apply the centered stencil on the centered box
     2. if(hasLo) (Lo stencil gets applied on the lo box )
     3. if(hasHi) (Hi stencil gets applied on the hi box )
     This uses the standard chombo semantic of dest, src.

     Note 1: The source for the lo and hi are not necessarily the same as the center
     stencil's source so we send the low and high source as  pointers.    

     Note 2: Lo and Hi stencil execution does not always
     happen so those pointers can be null. 

     Note 3: In the viscous tensor operator,
     velocity stencil applies rely upon ghost cells; 
     grad stencil applies rely upon one-sided-ness.    

     Note 4: Yes it does really work that way.
     Make fun if you like, but this operator has been widely used. 
     Good engineering is not always elegant.

     Note 5: In that vein, I am unsatisfied with the LoHiCenter API 
     but I cannot think of anything better. 
  **/
  int  loHiCenterApply(pr_box_data_sca         & a_dst,
                       const pr_box_data_sca   & a_src_ce,
                       const pr_box_data_sca   * a_src_lo_ptr,
                       const pr_box_data_sca   * a_src_hi_ptr,
                       const LHC_MetaData      & a_sten) 
  {
    CH_TIME("loHiCenterApply");
    ///always.   
    bool initToZero = true; double applyScale = 1;
    //this is src, dest so a_dst is the output
    {
      CH_TIME("actual_cen_apply");
      a_sten.m_ce_sten.apply(  a_src_ce, a_dst, a_sten.m_ce_box, initToZero, applyScale);
    }

    if(a_sten.m_has_lo == 1)
    {
      CH_TIME("actual_lo_apply");
      if(a_src_lo_ptr == NULL)
      {
        ch_mayday::Error("loHiCenterApply: Error: has_lo == 1 and src_lo_ptr == NULL");
      }
      //this is src, dest so a_dst is the output
      a_sten.m_lo_sten.apply(*a_src_lo_ptr, a_dst, a_sten.m_lo_box, initToZero, applyScale);
    }
    if(a_sten.m_has_hi == 1)
    {
      CH_TIME("actual_hi_apply");
      if(a_src_hi_ptr == NULL)
      {
        ch_mayday::Error("loHiCenterApply: Error: has_hi == 1 and src_hi_ptr == NULL");
      }
      //this is src, dest so a_dst is the output
      a_sten.m_hi_sten.apply(*a_src_hi_ptr, a_dst, a_sten.m_hi_box, initToZero, applyScale);
    }
    return 0;
  }//end function loHiCenterApply

} //end namespace Pr_Ch_AMR_AMRElliptic
#endif
