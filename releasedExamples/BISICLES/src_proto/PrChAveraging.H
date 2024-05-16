#ifndef __PR_CH_Averaging_
#define __PR_CH_Averaging_
#include "Proto.H"

namespace PrChTools
{
  ///
  class InterpolationUtilities
  {
  public:
    typedef  Proto::ProblemDomain                  pr_dom;
    typedef  Proto::Point                          pr_pt;
    typedef  Proto::Shift                          pr_shift;
    typedef  Proto::Box                            pr_box;
    typedef  Proto::BoxData<double, 1    >         pr_box_data;
    typedef  Proto::Stencil<double>                pr_sten;
    
    ///zeroth order interpolation
    static vector<pr_sten>
    getZerothInterpStencils(int    a_ref_rat)
    {
      pr_box ref_box(pr_pt::Zeros(), pr_pt::Ones(a_ref_rat -1));
      int numPts = ref_box.size();
      vector<pr_sten> retval(numPts);

      pr_pt dst_ref_ratio = pr_pt::Ones(a_ref_rat);
      pr_pt src_ref_ratio = pr_pt::Ones(1);
      pr_shift   home_shift(pr_pt::Zeros());
      pr_sten    home_sten(home_shift, 1.);
      int ipt = 0;
      for(auto boxit = ref_box.begin(); boxit != ref_box.end(); ++boxit)
      {
        pr_pt box_loc = *boxit;
        pr_sten pt_sten  =  home_sten;

        pt_sten.srcRatio()  = src_ref_ratio;
        pt_sten.destRatio() = dst_ref_ratio;
        pt_sten.destShift() = box_loc;
        retval[ipt]= pt_sten;
        ipt++;
      }
    
      return retval;
    }
  };
  class AveragingUtilities
  {
  public:
    ///no data. all static functions
    AveragingUtilities(){;}
    ~AveragingUtilities(){;}

    typedef  Proto::ProblemDomain                  pr_dom;
    typedef  Proto::Point                          pr_pt;
    typedef  Proto::Shift                          pr_shift;
    typedef  Proto::Box                            pr_box;
    typedef  Proto::BoxData<double, 1    >         pr_box_data;
    typedef  Proto::Stencil<double>                pr_sten;

    ///
    static pr_sten
    getStencilFromRefBox(pr_box a_ref_box,
                         int    a_ref_rat)
    {
      int numPts = a_ref_box.size();
      double weight = 1.0/double(numPts);
      pr_pt src_ref_ratio = pr_pt::Ones(a_ref_rat);
      pr_pt dst_ref_ratio = pr_pt::Ones(1);
      pr_pt dst_shift     = pr_pt::Zeros();
      /**
         Yes, I know I am setting things on an  empty stencil.
         Right now, adding a stencil with non unity refinement ratio 
         to an empty stencil without doing this causes a runtime error.
         dtg 11-28-2023
      **/
      pr_sten whole_sten;
      whole_sten.srcRatio()  = src_ref_ratio;
      whole_sten.destRatio() = dst_ref_ratio;
      for(auto boxit = a_ref_box.begin(); boxit != a_ref_box.end(); ++boxit)
      {
        pr_pt ref_pt = *boxit;
        pr_sten pt_sten(pr_shift(ref_pt), weight, dst_ref_ratio, dst_shift, src_ref_ratio);
        whole_sten += pt_sten;
      }
    
      return whole_sten;
    }
    /// for averaging stuff based in faces 
    /**
       set coarse value to arithmetic average of all fine face values.   Only faces in the same facedir plane are used.
       The valid box is cell centered.
    **/
    static void
    averageDownFace(pr_box_data       & a_coar,
                    const pr_box_data & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat,
                    int                 a_face_dir)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Zeros());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      // only include faces on the same facedir plane
      ref_box.m_high[a_face_dir] = 0;
      ref_box.recomputeSize();  // required to keep box data consistent after above shenanigans.
      
      pr_sten whole_sten = getStencilFromRefBox(ref_box, a_ref_rat);

      pr_box_data& fine_cast = (pr_box_data&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      pr_box coar_face_valid = a_coar_valid;
      coar_face_valid.m_high[a_face_dir] += 1;
      coar_face_valid.recomputeSize();  // required to keep box data consistent after above shenannigans.
      
      whole_sten.apply(fine_cast, a_coar, coar_face_valid, initToZero, applyScale);

    }
    /// for averaging stuff based in cells
    /**
       set coarse value to arithmetic average of all fine cell values.
    **/
    static void
    averageDownCell(pr_box_data       & a_coar,
                    const pr_box_data & a_fine,
                    const pr_box      & a_coar_valid,
                    int                 a_ref_rat)
    {
      pr_box unit_box(pr_pt::Zeros(), pr_pt::Zeros());
      pr_pt pt_ref_rat = pr_pt::Ones( a_ref_rat);
      pr_box ref_box = unit_box.refine(pt_ref_rat);
      pr_sten whole_sten = getStencilFromRefBox(ref_box, a_ref_rat);

      
      pr_box_data& fine_cast = (pr_box_data&) a_fine;
      bool initToZero = true; double applyScale = 1.;
      whole_sten.apply(fine_cast, a_coar, a_coar_valid, initToZero, applyScale);
    }
  };  //end class AveragingUtilities
}     //end namespace Chombo

#endif //include guard
