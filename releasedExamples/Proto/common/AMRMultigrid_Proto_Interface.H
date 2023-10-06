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

    ///This class has no data
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
    refluxFineCorrection(level_data_type       &  a_phiFine, 
                         const level_data_type &  a_phiCoar,
                         const ch_dom          &  a_domFine,
                         const ch_dom          &  a_domCoar,
                         int                      a_refRatio) = 0;
    
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

    /// The stuff below here should reduce to the stuff above
    /// LinearOp  functions

    ///lhs = rhs - l(phi)
    virtual void
    residual(level_data_type       & a_lhs,
             const level_data_type & a_phi,
             const level_data_type & a_rhs,
             bool a_homogeneous = false)
    {
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
    virtual void applyOp(   T& a_lhs, const T& a_phi, bool a_homogeneous = false)
    {
      ch_dom domain = getDomainFromData(a_phi);
      level_data_type& castphi = (level_data_type&)(a_phi);
      fillGhostCellsOutsideDomain(castphi);
      applyOperator(a_lhs, a_phi);
    }


    ///
    virtual void assign(    T& a_lhs, const T& a_rhs)
    {
      //a_lhs = 0. * a_lhs + 1.* a_rhs
      axby(a_lhs, a_lhs, a_rhs, 0.0, 1.0, );
    }

    ///
    /**
       Increment by scaled amount (a_lhs += a_scale*a_x).
    */
    virtual void
    incr  (level_data_type       & a_lhs,
           const level_data_type & a_x,
           double a_scale)
    {
      ///a_lhs = 1. a_lhs + a_scale*a_x + a_lhs
      axby(a_lhs, a_lhs , a_x, 1.0,  a_scale);
    }

    ///
     virtual void scale(     T& a_lhs, const Real& a_scale)
    {
      ///a_lhs = 0. *a_lhs + a_scale*a_lhs
      axby(a_lhs, a_lhs , a_lhs, 0.0,  a_scale);
    }


    ///MGLevelOp interface

    ///
    /**
       Create a coarsened  (by two) version of the input data.  This does not include averaging
       the data.   So if a_fine is over a Box of (0, 0, 0) (63, 63, 63), a_fine should
       be over a Box (0, 0, 0) (31, 31, 31).
    */
    virtual void createCoarser(T& a_coarse, const T& a_fine, bool ghosted) = 0;

    ///
    /**
       Use your relaxation operator to remove the high frequency wave numbers from
       the correction so that it may be averaged to a coarser refinement.
       A point relaxtion scheme, for example takes the form
       a_correction -= lambda*(L(a_correction) - a_residual).
    */

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



////AMRLevelOp interface
    ///
    /**
       return the refinement ratio to next coarser level.
       return 1 when there are no coarser Base_Proto_AMRLevelOp objects.
    */
    virtual int refToCoarser() = 0;

    ///
    /**
       a_residual = a_rhs - L(a_phiFine, a_phi, a_phiCoarse)
    */
    virtual void AMRResidual(T& a_residual, const T& a_phiFine, const T& a_phi,
                             const T& a_phiCoarse, const T& a_rhs,
                             bool a_homogeneousDomBC,
                             Base_Proto_AMRLevelOp<T>*  a_finerOp) = 0;

    ///
    /**
       a_residual = a_rhs - L^nf(a_phi, a_phiCoarse)
       assume no finer AMR level
    */
    virtual void AMRResidualNF(T& a_residual, const T& a_phi, const T& a_phiCoarse,
                               const T& a_rhs, bool a_homogeneousBC) = 0;

    ///
    /**
       a_residual = a_rhs - L(a_phiFine, a_phi)
       assume no coarser AMR level
    */
    virtual void AMRResidualNC(T& a_residual, const T& a_phiFine, const T& a_phi,
                               const T& a_rhs, bool a_homogeneousBC,
                               Base_Proto_AMRLevelOp<T>* a_finerOp) = 0;

    ///
    /**
       Apply the AMR operator, including coarse-fine matching
    */
    virtual void AMROperator(T& a_LofPhi,
                             const T& a_phiFine, const T& a_phi,
                             const T& a_phiCoarse,
                             bool a_homogeneousDomBC,
                             Base_Proto_AMRLevelOp<T>*  a_finerOp) = 0;

    ///
    /**
       Apply the AMR operator, including coarse-fine matching.
       assume no finer AMR level
    */
    virtual void AMROperatorNF(T& a_LofPhi,
                               const T& a_phi,
                               const T& a_phiCoarse,
                               bool a_homogeneousBC) = 0;

    ///
    /**
       Apply the AMR operator, including coarse-fine matching
       assume no coarser AMR level
    */
    virtual void AMROperatorNC(T& a_LofPhi,
                               const T& a_phiFine,
                               const T& a_phi,
                               bool a_homogeneousBC,
                               Base_Proto_AMRLevelOp<T>* a_finerOp) = 0;

    ///
    /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection)) */
    virtual void AMRRestrict(T& a_resCoarse, const T& a_residual, const T& a_correction,
                             const T& a_coarseCorrection, bool a_skip_res) = 0;

    ///
    /** a_correction += I[2h->h](a_coarseCorrection) */
    virtual void AMRProlong(T& a_correction, const T& a_coarseCorrection) = 0;

    ///
    /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
    virtual void AMRUpdateResidual(T& a_residual, const T& a_correction,
                                   const T& a_coarseCorrection) = 0;


  };

///
/**
   Factory to create AMRLevelOps
*/
  template <class T>
  class Base_Proto_AMRLevelOpFactory : public AMRLevelOpFactory<T>
  {
  public:
    virtual ~Base_Proto_AMRLevelOpFactory()
    {
    }

    /**
       Create an operator at an  index space = coarsen(a_fineIndexSpace, 2^a_depth)
       Return NULL if no such Multigrid level can be created at this a_depth.
       If a_homoOnly = true, then only homogeneous boundary conditions will be needed.
    */
    virtual MGLevelOp<T>* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                  int                  a_depth,
                                  bool                 a_homoOnly = true) = 0;

    ///
    /**
       return a new operator.  this is done with a new call.
       caller is responsible for deletion
    */
    virtual AMRLevelOp<T>* AMRnewOp(const ProblemDomain& a_indexSpace)=0;

    ///
    /**
       return refinement ratio  to next finer level.
    */
    virtual int refToFiner(const ProblemDomain& a_indexSpace) const =0;

  };
  
  ///
  /**
     This layer deals with proto<---->chombo conversions  and will deal only with
     the purely virtual functions of Base_Simplified_AMRLevelOp.
     This layer does heavily use the proto interface.   Only operator-dependent stuff
     is left out.   I think this means that we can have the application operator only 
     have to make a flux or apply the operator.
  **/
  template<int ncomp>
  class Base_Proto_AMRLevelOp: Base_Simplified_AMRLevelOp<Proto::LevelBoxData<ncomp> >
  {
  public:
    ///this holds data that is common to all of these things.

  };


#endif
