#include "AMRMultiGrid.H"
#include "Proto.H"

#ifndef AMRMultigrid_Proto_Interface

///
/**
   Infrastructurue for using AMRMultiGrid as a solver front end and a
   proto back end.  I intend to rewrite in this fashion all of the
   some of the non-EB elliptic operators in Chombo3.2, including
   Poisson (constant coefficient and variable), ViscousTensor and
   Resistivity (variable coefficient for both).  I ignore all the
   use_homogeneous flags in the base classes because all domain bcs
   here should be homogeneous.  If you have inhomogeneous boundary
   conditions, slip them into the RHS.  Your data holder needs to be
   able to give its domain.  

   Basic Inheritance Structure will look like this:
   OperatorClass: Base_Proto_AMRLevelOp: Base_Simplified_AMRLevelOp

   Base_Simplified_AMRLevelOp funnels down the  AMRLevelOp interface 
   to a more managable number of functions.

   Base_Proto_AMRLevelOp uses proto functions to narrow that number of functions down further.

   OperatorClass should only have to apply its operator, fill in ghost cell boundary 
   conditions outside the domain, and 
   to a more managable number of functions.
**/
namespace PrCh_AMR_Elliptic
{
  typedef Chombo::Box  ch_box;

  ///
  /**
     This class includes all the virtual functions that AMRLevelOp, MGLevelOp and LevelOp require.
     I funnel them all down into a few pure virtual functions that each separate operator implements.
     Everything is implemented in the simplest way I can figure out so lots of niggling optimizations
     have been ignored (axby is used for all arithmetic, for example).    If someone uses this interface
     enough to want them back, I will consider that a victory.

     The number of purely virtual functions here is still large but this has no references to proto. 
     The next layer deals with proto<---->chombo conversions  and will deal only with the functions here
     which are purely virtual in this class.    It will use proto tools to reduce further the number of 
     functions required by operator.
  
     
     A purely virtual function's = 0 means that the function has not been sullied with an implementation.
     A purely virtuous function is a completely different thing.  
   **/
  template<typename level_data_type>
  class Base_Simplified_AMRLevelOp : public Chombo::AMRLevelOp<level_data_type>
  {
  public:

    typedef Chombo::Box           ch_box;
    typedef Chombo::ProblemDomain ch_dom;

    ///This class has no data.  
    Base_Simplified_AMRLevelOp():AMRLevelOp<>() {}

    ///This class has no data
    virtual ~Base_Simplified_AMRLevelOp(){}

    ///functions that all this will funnel into.
    virtual void
    relax(level_data_type       & a_correction,
          const level_data_type & a_residual, int a_iterations) = 0 ;

    ///
    virtual ch_dom
    getDomainFromData(const level_data_type & a_data) = 0;
    
    /// 
    virtual void
    fillGhostDataOutsideDomain(level_data_type & a_phi,
                               const ch_dom    & a_domain) = 0;

    ///
    virtual void
    fillGhostDataOverCoarseFineInterface(level_data_type       &  a_phiFine, 
                                         const level_data_type &  a_phiCoar) = 0;
    ///
    virtual void
    refluxFineCorrection(level_data_type       &  a_lphCoar,
                         const level_data_type &  a_phiFine, 
                         const level_data_type &  a_phiCoar,
                         Base_Proto_AMRLevelOp<level_data_type>*  a_finerOp) = 0;
    
    ///exchange
    virtual void
    fillGhostDataOverFineFineInterface(level_data_type       &  a_phiFine) = 0;

    ///these two can be on different layouts
    virtual void
    copyTo(level_data_type       &  a_dst,
           const level_data_type &  a_src) = 0;
    
    
    ///
    /**
       Create a coarsened  (by two) version of the input data.  This does not include averaging
       the data.   So if a_fine is over a Box of (0, 0, 0) (63, 63, 63), a_fine should
       be over a Box (0, 0, 0) (31, 31, 31).
    */
    virtual void createCoarser(T& a_coarse, const T& a_fine, bool ghosted) = 0;
    
    ///
    /**
       Return the norm of  a_rhs.
       a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
    */
    virtual Real norm(const T& a_rhs, int a_ord) = 0;
    
    ///
    /**
       Compute and return the dot product of a_1 and a_2.   In most contexts, this
       means return the sum over all data points of a_1*a_2.
    */
    virtual double
    dotProduct(const level_data_type & a_fir,
               const level_data_type & a_sec)  = 0;
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(level_data_type         & a_lph,
                  const level_data_type   & a_phi) = 0;


    ///
    /**
       Create data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
       (just  make a holder the same size).
    **/
    virtual void
    create(level_data_type       & a_lhs,
           const level_data_type & a_rhs) = 0;
    
    ///
    /**
       This is another one where you just have to create the data (call level_data_type::define, for example)
       And every box of Coar is a coarsening of the same box on fine by refrat.
    */
    virtual void createCoarsened(level_data_type&       a_phiCoar,
                                 const level_data_type& a_phiFine,
                                 const int&     a_refRat) = 0;
    
    ///
    /**
       calculate restricted residual
       a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
    */
    virtual void restrictResidual(T& a_resCoarse, T& a_phiFine, const T& a_rhsFine) = 0;

    ///
    /**
       correct the fine solution based on coarse correction
       a_phiThisLevel += I[2h->h](a_correctCoarse)
    */
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;
    
    ///a_resCoarse = average(a_resfine) (cannot call this restrict)
    virtual void restriction(T& a_resCoarse, T& a_resFine) = 0;
    
    ///
    /**
       return the refinement ratio to next coarser level.
       return 1 when there are no coarser Base_Proto_AMRLevelOp objects.
    */
    virtual int refToCoarser() = 0;
    
    ///to make arithmetic sane
    /**
       a_output = a_aco*a_xval + a_bco*a_yval.
       This is not fortran so different arguments can point to the same data.
    **/
    virtual void 
    axby(level_data_type         & a_output,
         const level_data_type   & a_xval,
         const level_data_type   & a_yval,
         const double            & a_bco,
         const double            & a_aco) = 0;
         
    ///
    virtual void setToZero(level_data_type & a_lhs) = 0;

    ///
    /**
       correct the fine solution based on coarse correction
       a_phiThisLevel += I[2h->h](a_correctCoarse)
    */
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;
  
    /// The stuff below here should reduce to the stuff above
    /// LinearOp  functions

    ///
    /**
       calculate restricted residual
       a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
    */
    virtual void restrictResidual(level_data_type      & a_resCoar,
                                  level_data_type      & a_phiMedi,
                                  const level_data_type& a_rhsMedi)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::restrictResidual");
      ///this makes a fine version of the residual
      level_data_type resMedi;
      create(resMedi, a_rhsMedi);
      residual(resMedi, a_phiMedi, a_rhsMedi);
      restriction(resCoar, resMedi);
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
      applyOperator(a_lhs, a_phi);
      ///lhs = rhs - lphi
      axby(a_lhs, a_rhs, a_lhs, 1.0, -1.0);
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

      relax(a_cor, a_res, num_relax);
    }

    ///
    virtual void
    applyOp(level_data_type      & a_lhs,
            const level_data_type& a_phi,
            bool a_homogeneous = false)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::applyOp");
      ch_dom domain = getDomainFromData(a_phi);
      level_data_type& castphi = (level_data_type&)(a_phi);
      fillGhostDataOutsideDomain(castphi, domain);
      fillGhostDataOverFineFineInterface(castphi);
      applyOperator(a_lhs, a_phi);
    }


    ///
    virtual void assign(    level_data_type& a_lhs, const level_data_type& a_rhs)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::assign");
      //a_lhs = 0. * a_lhs + 1.* a_rhs
      axby(a_lhs, a_lhs, a_rhs, 0.0, 1.0, );
    }

    ///
    virtual void
    incr  (level_data_type       & a_lhs,
           const level_data_type & a_x,
           double a_scale)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::incr");
      ///a_lhs = 1. a_lhs + a_scale*a_x + a_lhs
      axby(a_lhs, a_lhs , a_x, 1.0,  a_scale);
    }

    ///
    virtual void
    scale(level_data_type      & a_lhs,
          const Real           & a_scale)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::scale");
      ///a_lhs = 0. *a_lhs + a_scale*a_lhs
      axby(a_lhs, a_lhs , a_lhs, 0.0,  a_scale);
    }


    ///Apply the AMR operator, including flux matching.
    virtual void AMROperator(level_data_type      & a_lphMedi,
                             const level_data_type& a_phiFine,
                             const level_data_type& a_phiMedi,
                             const level_data_type& a_phiCoar,
                             bool a_homogeneousDomBC,
                             Base_Proto_AMRLevelOp<level_data_type>*  a_finerOp)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMROperator");
      ch_dom domain = getDomainFromData(a_phiMedi);
      level_data_type& castphi = (level_data_type&)(a_phiMedi);
      fillGhostDataOutsideDomain(castphi, domain);
      fillGhostDataOverCoarseFineInterface(castphi, a_phiCoar);
      fillGhostDataOverFineFineInterface(castphi);
      applyOperator(a_lphMedi, a_phiMedi);
      refluxFineCorrection(a_lphMedi, a_phiFine, a_phiMedi, a_finerOp);
    }
    
    ///
    virtual void
    AMRResidual(level_data_type      & a_residual,
                const level_data_type& a_phiFine,
                const level_data_type& a_phiMedi,
                const level_data_type& a_phiCoar,
                const level_data_type& a_rhsMedi,
                bool a_homogeneousDomBC,
                Base_Proto_AMRLevelOp<level_data_type>*  a_finerOp)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRResidual");
      /// residual here holds lph
      AMROperator(a_residual, a_phiFine, a_phiMedi, a_phiCoar, a_finerOp);
      /// residual = a_rhs - lphi
      axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);
      
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
      axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);
    }

    ///assume no coarser AMR level
    virtual void AMRResidualNC(level_data_type      & a_resMedi,
                               const level_data_type& a_phiFine,
                               const level_data_type& a_phiMedi,
                               const level_data_type& a_rhsMedi,
                               bool a_homogeneousBC,
                               Base_Proto_AMRLevelOp<level_data_type>* a_finerOp)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRResidualNC");
      /// after this, a_resMedi holds lphi
      AMROperatorNC(a_resMedi, a_phiFine, a_phiMedi, a_homogeneousBC, a_finerOp);
      /// residual = a_rhs - lphi
      axby(a_resMedi, a_rhsMedi, a_resMedi, 1.0, -1.0);

    }


    ///assume no finer level
    virtual void AMROperatorNF(level_data_type      & a_lphMedi,
                               const level_data_type& a_phiMedi,
                               const level_data_type& a_phiCoar,
                               bool a_homogeneousBC)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMROperatorNF");
      ch_dom domain = getDomainFromData(a_phiMedi);
      level_data_type& castphi = (level_data_type&)(a_phiMedi);
      fillGhostDataOutsideDomain(castphi, domain);
      fillGhostDataOverCoarseFineInterface(castphi, a_phiCoar);
      applyOperator(a_lphMedi, a_phiMedi);
      
    }

    ///assume no coarser AMR level
    virtual void
    AMROperatorNC(level_data_type      & a_lphMedi,
                  const level_data_type& a_phiFine,
                  const level_data_type& a_phiMedi,
                  bool a_homogeneousBC,
                  Base_Proto_AMRLevelOp<level_data_type>* a_finerOp)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMROperatorNC");
      ch_dom domain = getDomainFromData(a_phiMedi);
      level_data_type& castphi = (level_data_type&)(a_phiMedi);
      fillGhostDataOutsideDomain(        castphi, domain);
      fillGhostDataOverFineFineInterface(castphi);
      applyOperator(a_lphMedi, a_phiMedi);
    }

 
    ///  a_resCoar = I[h-2h]( a_rhsMedi - L(phiMedi, phiCoar)  )
    virtual void
    AMRRestrict(level_data_type      & a_resCoar,
                const level_data_type& a_rhsMedi,
                const level_data_type& a_phiMedi,
                const level_data_type& a_phiCoar
                bool a_skip_res)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRRestrict");
      level_data_type resMedi;
      create(resMedi, a_rhsMedi);
      AMROperatorNF(resMedi, a_phiMedi. a_phiCoar);
      restricton(a_resCoar, resMedi);
    }

    ///a_phiFine += I[2h->h](a_phiCoar) 
    virtual void
    AMRProlong(level_data_type      & a_phiFine,
               const level_data_type& a_phiCoar)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRProlong");
      level_data_type phiCoFi;
      createCoarsened(phiCoFi, a_phiFine);
      copyTo(phiCoFi, a_phiCoar);
      prolongIncrement(a_phiFine, phiCoFi);
    }

    ///a_resMedi = a_resMedi - L(a_phiMedi, a_phiCoar)
    virtual void
    AMRUpdateResidual(level_data_type      & a_resMedi,
                      const level_data_type& a_phiMedi,
                      const level_data_type& a_phiCoar)
    {
      CH_TIME("Base_Simplified_AMRLevelOp::AMRUpdateResidual");
      level_data_type lphMedi;
      create(lphMedi, a_resMedi);
      AMROperatorNF(lphMedi, a_phiMedi, a_phiCoar);
      axby(a_resMedi, a_resMedi, lphMedi, 1.0, -1);
    }


  };


  ///
  /**
     level_data_typehis layer deals with proto<---->chombo conversions  and will deal only with
     the purely virtual functions of Base_Simplified_AMRLevelOp.
     This layer does heavily use the proto interface.   Only operator-dependent stuff
     is left out.   I think this means that we can have the application operator only 
     have to make a flux or apply the operator.
  **/
  template<int ncomp>
  class Base_Proto_AMRLevelOp: public Base_Simplified_AMRLevelOp<Proto::LevelBoxData<ncomp> >
  {
  public:
    typedef Proto::LevelBoxData<ncomp> level_data_type;
    
    ///this holds data that is common to all of these things.
    typedef Chombo::Box           ch_box;
    typedef Chombo::ProblemDomain ch_dom;

    ///This class has no data.  
    Base_Simplified_AMRLevelOp():AMRLevelOp<>() {}

    ///This class has no data
    virtual ~Base_Simplified_AMRLevelOp(){}

    ///functions that all this will funnel into.
    virtual void
    relax(level_data_type       & a_correction,
          const level_data_type & a_residual, int a_iterations) = 0 ;

    ///
    virtual ch_dom
    getDomainFromData(const level_data_type & a_data) = 0;
    
    /// 
    virtual void
    fillGhostDataOutsideDomain(level_data_type & a_phi,
                               const ch_dom    & a_domain) = 0;

    ///
    virtual void
    fillGhostDataOverCoarseFineInterface(level_data_type       &  a_phiFine, 
                                         const level_data_type &  a_phiCoar) = 0;
    ///
    virtual void
    refluxFineCorrection(level_data_type       &  a_lphCoar,
                         const level_data_type &  a_phiFine, 
                         const level_data_type &  a_phiCoar,
                         Base_Proto_AMRLevelOp<level_data_type>*  a_finerOp) = 0;
    
    ///exchange
    virtual void
    fillGhostDataOverFineFineInterface(level_data_type       &  a_phiFine) = 0;

    ///these two can be on different layouts
    virtual void
    copyTo(level_data_type       &  a_dst,
           const level_data_type &  a_src) = 0;
    
    
    ///
    /**
       Create a coarsened  (by two) version of the input data.  This does not include averaging
       the data.   So if a_fine is over a Box of (0, 0, 0) (63, 63, 63), a_fine should
       be over a Box (0, 0, 0) (31, 31, 31).
    */
    virtual void createCoarser(T& a_coarse, const T& a_fine, bool ghosted) = 0;
    
    ///
    /**
       Return the norm of  a_rhs.
       a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
    */
    virtual Real norm(const T& a_rhs, int a_ord) = 0;
    
    ///
    /**
       Compute and return the dot product of a_1 and a_2.   In most contexts, this
       means return the sum over all data points of a_1*a_2.
    */
    virtual double
    dotProduct(const level_data_type & a_fir,
               const level_data_type & a_sec)  = 0;
    


    ///
    /**
       Create data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
       (just  make a holder the same size).
    **/
    virtual void
    create(level_data_type       & a_lhs,
           const level_data_type & a_rhs) = 0;
    
    ///
    /**
       This is another one where you just have to create the data (call level_data_type::define, for example)
       And every box of Coar is a coarsening of the same box on fine by refrat.
    */
    virtual void createCoarsened(level_data_type&       a_phiCoar,
                                 const level_data_type& a_phiFine,
                                 const int&     a_refRat) = 0;
    
    ///
    /**
       calculate restricted residual
       a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
    */
    virtual void restrictResidual(T& a_resCoarse, T& a_phiFine, const T& a_rhsFine) = 0;

    ///
    /**
       correct the fine solution based on coarse correction
       a_phiThisLevel += I[2h->h](a_correctCoarse)
    */
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;
    
    ///a_resCoarse = average(a_resfine) (cannot call this restrict)
    virtual void restriction(T& a_resCoarse, T& a_resFine) = 0;
    
    ///
    /**
       return the refinement ratio to next coarser level.
       return 1 when there are no coarser Base_Proto_AMRLevelOp objects.
    */
    virtual int refToCoarser() = 0;

  };

  /// our friend mr. 7 pt laplacian
  class Proto_Helmholtz_Op: public   Base_Simplified_AMRLevelOp<1>
  {
  public:
    typedef Proto::LevelBoxData<1> level_data_type;
    typedef Chombo::MayDay         ch_mayday;
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(level_data_type         & a_lph,
                  const level_data_type   & a_phi)
    {
      ch_mayday::Error("Doh! not implemented");
    }
  }


#endif
