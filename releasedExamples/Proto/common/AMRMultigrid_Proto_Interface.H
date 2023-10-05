
#ifndef AMRMultigrid_Proto_Interface

///
/**
   Infrastructurue for using AMRMultigrid as a solver front end 
   and a proto back end.  These amount to a reimplementation of the
   non-EB elliptic operators in Chombo.   
 */
template<typename level_data_type>
class Base_Proto_AMRLevelOp : public Chombo::AMRLevelOp<level_data_type>
{
public:


  
  ///
  Base_Proto_AMRLevelOp()
    :MGLevelOp<T>()
  {
  }

  ///
  virtual ~Base_Proto_AMRLevelOp()
  {
  }
  /// LinearOp  functions

  ///
  /**
     Say you are  solving L(phi) = rhs.   Make a_lhs = L(a_phi) - a_rhs.   If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void residual(  T& a_lhs, const T& a_phi, const T& a_rhs, bool a_homogeneous = false) = 0;

  ///
  /**
     Given the current state of the residual the correction, apply your preconditioner to a_cor.
   */
  virtual void preCond(   T& a_cor, const T& a_residual)       = 0;

  ///
  /**
     In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).  If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void applyOp(   T& a_lhs, const T& a_phi, bool a_homogeneous = false) = 0;

  ///
  /**
     Creat data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
     just  make a holder the same size.
   */
  virtual void create(    T& a_lhs, const T& a_rhs) = 0;

  ///
  /**
     Set a_lhs  equal to a_rhs.
   */
  virtual void assign(    T& a_lhs, const T& a_rhs)       = 0;

  ///
  /**
     Compute and return the dot product of a_1 and a_2.   In most contexts, this
     means return the sum over all data points of a_1*a_2.
  */
  virtual Real dotProduct(const T& a_1, const T& a_2)     = 0;
  
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(const T& a_1, const int a_sz, const T a_2[], Real a_mdots[])
  {
    for (int j=0; j<a_sz; j++)
      {
        a_mdots[j] = dotProduct(a_1, a_2[j]);
      }
  }

  ///
  /**
     Increment by scaled amount (a_lhs += a_scale*a_x).
  */
  virtual void incr  (    T& a_lhs, const T& a_x, Real a_scale) = 0;

  ///
  /**
     Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
  */
  virtual void axby(      T& a_lhs, const T& a_x, const T& a_y, Real a_a, Real a_b) = 0;

  ///
  /**
     Multiply the input by a given scale (a_lhs *= a_scale).
  */
  virtual void scale(     T& a_lhs, const Real& a_scale)  = 0;

  ///
  /**
     Return the norm of  a_rhs.
     a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
  */
  virtual Real norm(const T& a_rhs, int a_ord) = 0;

  ///
  /**
     Set a_lhs to zero.
   */
  virtual void setToZero(T& a_lhs) = 0;

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
  virtual void relax(T& a_correction, const T& a_residual, int a_iterations) = 0 ;

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

  ///
  /**
   */
  virtual void createCoarsened(T&       a_lhs,
                               const T& a_rhs,
                               const int&     a_refRat) = 0;

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


#endif
