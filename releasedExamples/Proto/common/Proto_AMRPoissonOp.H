
#ifndef _AMRPOISSONOP_H_
#define _AMRPOISSONOP_H_


#include "AMRMultiGrid.H"
#include "Proto.H"

///
/**
   Operator for solving (alpha I + beta*Laplacian)(phi) = rho
   over an AMR hierarchy.
*/
class ProtoAMRPoissonOp : public AMRLevelOp<Proto::LevelBoxData<Real> >
{
public:

/** 
AMRMultgrid defines an operator interface.

These are the virtual functions I have to deal with:
LinearOp (from LinearSolver.H)":

  virtual void residual(  T& a_lhs, const T& a_phi, const T& a_rhs, bool a_homogeneous = false) = 0;
  virtual void preCond(   T& a_cor, const T& a_residual)       = 0;
  virtual void applyOp(   T& a_lhs, const T& a_phi, bool a_homogeneous = false) = 0;
  virtual void create(    T& a_lhs, const T& a_rhs) = 0;
  virtual void assign(    T& a_lhs, const T& a_rhs)       = 0;
  virtual Real dotProduct(const T& a_1, const T& a_2)     = 0;
  virtual void incr  (    T& a_lhs, const T& a_x, Real a_scale) = 0;
  virtual void axby(      T& a_lhs, const T& a_x, const T& a_y, Real a_a, Real a_b) = 0;
  virtual void scale(     T& a_lhs, const Real& a_scale)  = 0;
  virtual Real norm(const T& a_rhs, int a_ord) = 0;
  virtual void setToZero(T& a_lhs) = 0;

MGLevelOp (from MultiGrid.H):

  virtual void createCoarser(T& a_coarse, const T& a_fine, bool ghosted) = 0;
  virtual void relax(T& a_correction, const T& a_residual, int a_iterations) = 0 ;
  virtual void restrictResidual(T& a_resCoarse, T& a_phiFine, const T& a_rhsFine) = 0;
  virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;

  AMRLevelOp (from AMRMultigrid.H":

  virtual int refToCoarser() = 0;
  virtual void AMRResidual(T& a_residual, const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse, const T& a_rhs,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  virtual void AMRResidualNF(T& a_residual, const T& a_phi, const T& a_phiCoarse,
                             const T& a_rhs, bool a_homogeneousBC) = 0;

  virtual void AMRResidualNC(T& a_residual, const T& a_phiFine, const T& a_phi,
                             const T& a_rhs, bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;

  virtual void AMROperator(T& a_LofPhi,
                           const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  virtual void AMROperatorNF(T& a_LofPhi,
                             const T& a_phi,
                             const T& a_phiCoarse,
                             bool a_homogeneousBC) = 0;

  virtual void AMROperatorNC(T& a_LofPhi,
                             const T& a_phiFine,
                             const T& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;

  virtual void AMRRestrict(T& a_resCoarse, const T& a_residual, const T& a_correction,
                           const T& a_coarseCorrection, bool a_skip_res) = 0;

  virtual void AMRProlong(T& a_correction, const T& a_coarseCorrection) = 0;

  virtual void AMRUpdateResidual(T& a_residual, const T& a_correction,
                                 const T& a_coarseCorrection) = 0;

  virtual void createCoarsened(T&       a_lhs,
                               const T& a_rhs,
                               const int&     a_refRat) = 0;
**/


  typedef Proto::LevelBoxData<Real> pr_lbd;

  ///
  virtual void residual(  pr_lbd& a_lhs, const pr_lbd& a_phi, const pr_lbd& a_rhs, bool a_homogeneous = false)

  {
  }

  ///
  virtual void preCond(   pr_lbd& a_cor, const pr_lbd& a_residual)

  {
  }
  
  ///
  virtual void applyOp(   pr_lbd& a_lhs, const pr_lbd& a_phi, bool a_homogeneous = false)
  {
  }

  ///
  virtual void create(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
  {
  }

  ///
  virtual void assign(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
  {
  }

  ///
  virtual Real dotProduct(const pr_lbd& a_1, const pr_lbd& a_2)
  {
  }

  ///
  virtual void incr  (    pr_lbd& a_lhs, const pr_lbd& a_x, Real a_scale)
  {
  }

  ///
  virtual void axby(      pr_lbd& a_lhs, const pr_lbd& a_x, const pr_lbd& a_y, Real a_a, Real a_b)
  {
  }
  
  virtual void scale(     pr_lbd& a_lhs,   const Real& a_scale)
  {
  }
  
  virtual Real norm(const pr_lbd& a_rhs, int a_ord)
  {
  }
  
  virtual void setpr_lbdoZero(pr_lbd& a_lhs)
  {
  }

///MGLevelOp (from MultiGrid.H):

  ///
  virtual void createCoarser(pr_lbd& a_coarse, const pr_lbd& a_fine, bool ghosted)
  {
  }

  ///
  virtual void relax(pr_lbd& a_correction, const pr_lbd& a_residual, int a_iterations)
  {
  }
  
  ///
  virtual void restrictResidual(pr_lbd& a_resCoarse, pr_lbd& a_phiFine, const pr_lbd& a_rhsFine)
  {
  }

  
  virtual void prolongIncrement(pr_lbd& a_phiThisLevel, const pr_lbd& a_correctCoarse)
  {
  }

  ///
  virtual int refToCoarser()
  {
  }
  
  virtual void AMRResidual(pr_lbd& a_residual, const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                           const pr_lbd& a_phiCoarse, const pr_lbd& a_rhs,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<pr_lbd>*  a_finerOp) = 0;

  virtual void AMRResidualNF(pr_lbd& a_residual, const pr_lbd& a_phi, const pr_lbd& a_phiCoarse,
                             const pr_lbd& a_rhs, bool a_homogeneousBC) = 0;

  virtual void AMRResidualNC(pr_lbd& a_residual, const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                             const pr_lbd& a_rhs, bool a_homogeneousBC,
                             AMRLevelOp<pr_lbd>* a_finerOp) = 0;

  virtual void AMROperator(pr_lbd& a_LofPhi,
                           const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                           const pr_lbd& a_phiCoarse,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<pr_lbd>*  a_finerOp) = 0;

  virtual void AMROperatorNF(pr_lbd& a_LofPhi,
                             const pr_lbd& a_phi,
                             const pr_lbd& a_phiCoarse,
                             bool a_homogeneousBC) = 0;

  virtual void AMROperatorNC(pr_lbd& a_LofPhi,
                             const pr_lbd& a_phiFine,
                             const pr_lbd& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<pr_lbd>* a_finerOp) = 0;

  virtual void AMRRestrict(pr_lbd& a_resCoarse, const pr_lbd& a_residual, const pr_lbd& a_correction,
                           const pr_lbd& a_coarseCorrection, bool a_skip_res) = 0;

  virtual void AMRProlong(pr_lbd& a_correction, const pr_lbd& a_coarseCorrection) = 0;

  virtual void AMRUpdateResidual(pr_lbd& a_residual, const pr_lbd& a_correction,
                                 const pr_lbd& a_coarseCorrection) = 0;

  virtual void createCoarsened(pr_lbd&       a_lhs,
                               const pr_lbd& a_rhs,
                               const int&     a_refRat) = 0;


  /** going to try to do this with strong construction --placeholder
   */
  ProtoAMRPoissonOp(const ProblemDomain& a_indexSpace);

  ///
  /**
   */
  virtual ~ProtoAMRPoissonOp()
  {
  }
private:
  /** going to try to do this with strong construction
   */
  ProtoAMRPoissonOp();

};


#endif
