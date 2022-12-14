
#ifndef _AMRPOISSONOP_H_
#define _AMRPOISSONOP_H_


#include "AMRMultiGrid.H"
#include "Proto.H"
#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

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
private:
  /** going to try to do this with strong construction
   */
  ProtoAMRPoissonOp();

  
  unsigned int m_refToCoar;
  unsigned int m_refToFine;
  pr_ldb       m_resid;
  pr_ldb       m_residC;
  pr_ldb       m_deltaC; 
  pr_ldb       m_phic;

  bool   m_hasCoarser; //not strictly necessary--could check the pointer
  double m_alpha;     
  double m_beta;       
  double m_dx;        
  Box     m_domain;       
  double m_lambda; //relaxation coeff
  
  pr_lbd m_scratch;

  array<Point,MG_NUM_COLORS>  s_colors;
  Stencil<double> m_negoperator; //always need -Lap(phi)
  Stencil<double> m_restrict;

  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_updateOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpRhs;
  array<Stencil<double>, MG_NUM_COLORS> m_prolong;
public:

  typedef Proto::LevelBoxData<Real> pr_lbd;
  virtual void applyDomainBoundaryConditions (pr_ldb& a_phi)
  {
    
  }
  virtual void homogeneousCFInterp(pr_lbd& a_phi)
  {
    m_phiC.setVal(0);
    Proto::InterpStencil<Real> interpSten = Proto::InterpStencil<Real>::Build(2, m_refToCoar);
    a_phi.interpBoundaries(m_phiC, interpSten); 
  }

  ///
  virtual void residual(  pr_lbd& a_lhs, const pr_lbd& a_phi, const pr_lbd& a_rhs, bool a_homogeneous = false)
  {
    a_lhs.setVal(0.);
    applyOp(a_lhs, a_phi, a_homogeneous);
    a_lhs.increment(a_rhs, -1.);
  }

  ///
  virtual void preCond(   pr_lbd& a_cor, const pr_lbd& a_residual)
  {
    
    relax(a_cor, a_residual, 27);
  }

  ///
  virtual void applyOp(   pr_lbd& a_lhs, const pr_lbd& a_phi, bool a_homogeneous = false)
  {
    //fill ghost cells
    homogeneousCFInterp(a_phi);
    applyDomainBoundaryConditions(a_phi);
    auto dbl = a_lhs.layout();
    //apply the actual operator
    for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
    {
      auto valid = dbl[dit];
      a_lph[dit] |= m_negoperator(phi[dit], valid);
      a_lph *= -1.0;
    }
  }

  ///
  virtual void create(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
  {
    a_lhs.define(a_rhs.disjointBoxLayotu(), a_rhs.ghostVect());
  }

  ///
  virtual void assign(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
  {
    a_lhs.setVal(0.);
    a_lhs.increment(a_rhs);
  }

  ///
  virtual Real dotProduct(const pr_lbd& a_1, const pr_lbd& a_2)
  {
    pr_lbd& dotprod = m_scratch;
    Real sumdotprod = 0;
    for(auto dit = a_1.begin(); dit != a_1.end(); ++dit)
    {
      auto valid = dbl[dit];
      dotprod[dit].setVal(0.);
      dotprod[dit] += a_1[dit];
      dotprod[dit] *= a_2[dit];
      sumdotprod += dotprod[dit].sum();
    }

#ifdef CH_MPI
    Real sendBuf = sumdotprod;
    int result = MPI_Allreduce(&sendBuf, &sumdotprod, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Communication error in Proto_AMRPoissonOp::dotProduct");
    }
#endif

  }

  ///
  virtual void incr  (    pr_lbd& a_lhs, const pr_lbd& a_x, Real a_scale)
  {
    pr_lbd& increm = m_scratch;
    for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
    {
      increm[dit].setVal(0);
      increm[dit] += a_x[dit];
      increm[dit] *= a_scale;
      a_lhs[dit]  += increm[dit];
    }
  }

  ///
  virtual void axby(      pr_lbd& a_lhs, const pr_lbd& a_x, const pr_lbd& a_y, Real a_a, Real a_b)
  {
    a_lhs.setVal(0);
    a_lhs.incr(a_x, a_a);
    a_lhs.incr(a_y, a_b);
  }
  
  virtual void scale(     pr_lbd& a_lhs,   const Real& a_scale)
  {
    a_lhs
  }
  
  virtual Real norm(const pr_lbd& a_rhs, int a_ord)
  {
  }
  
  virtual void setToZero(pr_lbd& a_lhs)
  {
    a_lhs.setVal(0.);
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

};


#endif
