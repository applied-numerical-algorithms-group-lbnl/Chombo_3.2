#include "AMRMultiGrid.H"

#ifndef ___Base_Simplified_AMRLevelOp___
#define ___Base_Simplified_AMRLevelOp___

#include "NamespaceHeader.H"

/// Base_Simplified_AMRLevelOp provides reasonable defaults for many of AMRLevelOps functions 
/**
   class  Base_Simplified_AMRLevelOp

   This class includes all the virtual functions that AMRLevelOp, MGLevelOp and LevelOp require.
   I funnel them all down into a few pure virtual functions that each separate operator implements.
   Everything is implemented in the simplest way I can figure out so lots of niggling optimizations
   have been ignored (axby is used for all arithmetic, for example).    If someone uses this interface
   enough to want them back, I will consider that a victory.

   The number of purely virtual functions here is still large but this has no references to proto. 
   The next layer deals with proto<---->chombo conversions  and will deal only with the functions here
   which are purely virtual in this class.    It will use proto tools to reduce further the number of 
   functions required by the operator.
  
   Dad joke: 
   A purely virtual function's = 0 means that the function has not been sullied with an implementation.
   A purely virtuous function is a completely different thing.  
**/
template<typename level_data_type>
class Base_Simplified_AMRLevelOp : public AMRLevelOp<level_data_type>
{
public:

  typedef Box                            ch_box;
  typedef MayDay                         ch_mayday;
  typedef ProblemDomain                  ch_dom;
  typedef AMRLevelOp<level_data_type>    ch_amrlevelop;
  ///This class has no data.  
  Base_Simplified_AMRLevelOp():ch_amrlevelop()
  {}

  ///This class has no data 
  virtual ~Base_Simplified_AMRLevelOp(){}

    
         
  ///a_phiThisLevel += I[2h->h](a_correctCoarse)
  virtual void prolongIncrement(level_data_type      & a_phiThisLevel,
                                const level_data_type& a_correctCoarse) = 0;


  ///Default implementations independent of Proto go here.
    
  /// a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  virtual void restrictResidual(level_data_type      & a_resCoar,
                                level_data_type      & a_phiMedi,
                                const level_data_type& a_rhsMedi)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::restrictResidual");
    ///this makes a fine version of the residual
    level_data_type resMedi;
    this->create(resMedi, a_rhsMedi);
    residual(resMedi, a_phiMedi, a_rhsMedi);
    restriction(a_resCoar, resMedi, 2);
  }

  
    
  ///lhs = rhs - l(phi)
  virtual void
  residual(level_data_type       & a_lhs,
           const level_data_type & a_phi,
           const level_data_type & a_rhs,
           bool a_homogeneous = false)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::residual");
    ///here lhs = lphi
    a_lhs.setVal(0.);
    level_data_type lphi;
    this->create(lphi, a_rhs);
      
    applyOp(lphi, a_phi, a_homogeneous);
    ///lhs = rhs - lphi
    this->axby(a_lhs, a_rhs, lphi, 1.0, -1.0);
  }


  /// precond is just relax here.   This can be overwritten by the derived class.
  virtual void
  preCond(level_data_type        & a_cor,
          const level_data_type  & a_res)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::preCond");
    /**
       I have not seen it to be  all that important what 
       this number is (as long as it is not too high).
       I did once get grief for setting it to 40.
       So I guess we should keep this under 40. --dtg
    **/
    int num_relax = 27;

    this->relax(a_cor, a_res, num_relax);
  }

  ///
  virtual void
  applyOp(level_data_type      & a_lhs,
          const level_data_type& a_phi,
          bool a_homogeneous = false)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::applyOp");
    ch_dom domain = this->getDomainFromData(a_phi);
    level_data_type& castphi = (level_data_type&)(a_phi);
    fillGhostDataOutsideDomain(castphi, domain);
    fillGhostDataOverFineFineInterface(castphi);
    this->applyOperator(a_lhs, a_phi);
  }



  ///
  virtual void
  incr  (level_data_type       & a_lhs,
         const level_data_type & a_x,
         double a_scale)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::incr");
    ///a_lhs = 1. a_lhs + a_scale*a_x + a_lhs
    this->axby(a_lhs, a_lhs , a_x, 1.0,  a_scale);
  }

  ///
  virtual void
  scale(level_data_type      & a_lhs,
        const Real           & a_scale)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::scale");
    ///a_lhs = 0. *a_lhs + a_scale*a_lhs
    this->axby(a_lhs, a_lhs , a_lhs, 0.0,  a_scale);
  }


  ///Apply the AMR operator, including flux matching.
  virtual void
  AMROperator(level_data_type      & a_lphMedi,
              const level_data_type& a_phiFine,
              const level_data_type& a_phiMedi,
              const level_data_type& a_phiCoar,
              bool a_homogeneousDomBC,
              ch_amrlevelop*  a_finerOp)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMROperator");
    ch_dom domain = getDomainFromData(a_phiMedi);
    level_data_type& castphi = (level_data_type&)(a_phiMedi);
    this->fillGhostDataOutsideDomain(castphi, domain);
    this->fillGhostDataOverCoarseFineInterface(castphi, a_phiCoar);
    this->fillGhostDataOverFineFineInterface(castphi);
    applyOperator(a_lphMedi, a_phiMedi);
    this->refluxFineCorrection(a_lphMedi, a_phiFine, a_phiMedi, a_finerOp);
  }
    
  ///
  virtual void
  AMRResidual(level_data_type      & a_resMedi,
              const level_data_type& a_phiFine,
              const level_data_type& a_phiMedi,
              const level_data_type& a_phiCoar,
              const level_data_type& a_rhsMedi,
              bool a_homogeneousDomBC,
              ch_amrlevelop*  a_finerOp)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMRResidual");
    /// residual here holds lph
    this->AMROperator(a_resMedi, a_phiFine, a_phiMedi, a_phiCoar,
                      a_homogeneousDomBC, a_finerOp);
    /// residual = a_rhs - lphi
    this->axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);
      
  }

  ///assume no finer AMR level
  virtual void
  AMRResidualNF(level_data_type      & a_resMedi,
                const level_data_type& a_phiMedi,
                const level_data_type& a_phiCoar,
                const level_data_type& a_rhsMedi,
                bool a_homogeneousBC)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMRResidualNF");
    /// here resMedi = L(phiMedi)
    AMROperatorNF(a_resMedi, a_phiMedi, a_phiCoar, a_homogeneousBC);
    /// residual = a_rhs - lphi
    this->axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);
  }

  ///assume no coarser AMR level
  virtual void AMRResidualNC(level_data_type      & a_resMedi,
                             const level_data_type& a_phiFine,
                             const level_data_type& a_phiMedi,
                             const level_data_type& a_rhsMedi,
                             bool a_homogeneousBC,
                             ch_amrlevelop* a_finerOp)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMRResidualNC");
    /// after this, a_resMedi holds lphi
    AMROperatorNC(a_resMedi, a_phiFine, a_phiMedi, a_homogeneousBC, a_finerOp);
    /// residual = a_rhs - lphi
    this->axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);

  }


  ///assume no finer level
  virtual void AMROperatorNF(level_data_type      & a_lphMedi,
                             const level_data_type& a_phiMedi,
                             const level_data_type& a_phiCoar,
                             bool a_homogeneousBC)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMROperatorNF");
    ch_dom domain = this->getDomainFromData(a_phiMedi);
    level_data_type& castphi_medi = (level_data_type&)(a_phiMedi);
    level_data_type& castphi_coar = (level_data_type&)(a_phiCoar);
    fillGhostDataOutsideDomain(castphi_medi, domain);
    fillGhostDataOverCoarseFineInterface(castphi_medi, castphi_coar);
    applyOperator(a_lphMedi, a_phiMedi);
      
  }

  ///assume no coarser AMR level
  virtual void
  AMROperatorNC(level_data_type      & a_lphMedi,
                const level_data_type& a_phiFine,
                const level_data_type& a_phiMedi,
                bool a_homogeneousBC,
                ch_amrlevelop* a_finerOp)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMROperatorNC");
    ch_dom domain = getDomainFromData(a_phiMedi);
    level_data_type& castphi_medi = (level_data_type&)(a_phiMedi);
    level_data_type& castphi_fine = (level_data_type&)(a_phiFine);
    fillGhostDataOutsideDomain(        castphi_medi, domain);
    fillGhostDataOverFineFineInterface(castphi_medi);
    applyOperator(a_lphMedi, a_phiMedi);
    //finer level should always be defined if you got here.
    this->refluxFineCorrection(a_lphMedi, castphi_fine, castphi_medi, a_finerOp);
  }

 

  ///a_resMedi = a_resMedi - L(a_phiMedi, a_phiCoar)
  virtual void
  AMRUpdateResidual(level_data_type      & a_resMedi,
                    const level_data_type& a_phiMedi,
                    const level_data_type& a_phiCoar)
  {
    CH_TIME("Base_Simplified_AMRLevelOp::AMRUpdateResidual");
    level_data_type lphMedi;
    this->create(lphMedi, a_resMedi);
    AMROperatorNF(lphMedi, a_phiMedi, a_phiCoar, false);
    this->axby(a_resMedi, a_resMedi, lphMedi, 1.0, -1);
  }

  ///The following functions are new (not from AMRLevelOp).
    
  ///
  virtual ch_dom
  getDomainFromData(const level_data_type & a_data) = 0;

  ///
  virtual void
  fillGhostDataOutsideDomain(level_data_type & a_phi,
                             const ch_dom    & a_domain) = 0;

  ///
  virtual void
  fillGhostDataOverFineFineInterface(level_data_type &  a_phiFine) = 0;

  ///
  virtual void
  copyTo(level_data_type       &  a_dst,
         const level_data_type &  a_src) = 0;
    
  ///For this one, the derived class may assume all ghost cells have been set correctly
  virtual void 
  applyOperator(level_data_type        & a_lph,
                const level_data_type  & a_phi) = 0;

  /// Restrict fine data to coarse data using volume-weighted averaging.   
  /**
     Restrict is a c++ keyword, so we are restricted from naming a function restrict.
     To avoid this restriction, I named the function restriction instead of restrict.
  **/
  virtual void
  restriction(level_data_type & a_resCoar,
              level_data_type & a_resFine,
              int ref_rat) = 0;

  /// amr thing
  virtual void
  fillGhostDataOverCoarseFineInterface(level_data_type       &  a_phiFine, 
                                       const level_data_type &  a_phiCoar) = 0;
  /// amr thing
  virtual void
  refluxFineCorrection(level_data_type       &  a_lphCoar,
                       const level_data_type &  a_phiFine, 
                       const level_data_type &  a_phiCoar,
                       ch_amrlevelop*   a_finerOp)=0;
//  {
//    ch_mayday::Error("Doh! not implemented");
//  }
//
}; // end class Base_Simplified_AMRLevelOp

#include "NamespaceFooter.H"
  


#endif
