#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SGMULTIGRID_H_
#define _SGMULTIGRID_H_

#include <cmath>
#include <memory>
#include "Proto.H"
#include "LevelData.H"
#include "Copier.H"
#include "BaseFab.H"

#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

using  Proto::Stencil;
using  Proto::Point;
using std::array;
/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.  periodic boundary conditions.
class SGMultigridLevel
{
public:

  
  /// called from the finest level
  SGMultigridLevel(const double & a_alpha,
                   const double & a_beta,
                   const double & a_dx,
                   const CH_XD::DisjointBoxLayout & a_grids,
                   bool useDenseStencil);
                   

  ///  
  ~SGMultigridLevel()
  {
  }

  /// 
  void residual(CH_XD::LevelData<CH_XD::FArrayBox >       & a_res,
                const CH_XD::LevelData<CH_XD::FArrayBox > & a_phi,
                const CH_XD::LevelData<CH_XD::FArrayBox > & a_rhs);


  /// 
  void relax(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phi,
             const CH_XD::LevelData<CH_XD::FArrayBox >& a_rhs);


  /// average down residual to next coarser level.
  void restrictResidual(CH_XD::LevelData<CH_XD::FArrayBox >      & a_resc,
                        const CH_XD::LevelData<CH_XD::FArrayBox >& a_res);

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phiFine,
                        const CH_XD::LevelData<CH_XD::FArrayBox >& a_deltaCoarse);


  /// 
  void applyOp(CH_XD::LevelData<CH_XD::FArrayBox >       & a_lph,
               const CH_XD::LevelData<CH_XD::FArrayBox > & a_phi);

  /// Multigrid v-cycle.
  void vCycle(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phi,
              const CH_XD::LevelData<CH_XD::FArrayBox >& a_rhs);

  void
  enforceBoundaryConditions(const CH_XD::LevelData<CH_XD::FArrayBox >& a_phi)
  {
    CH_TIME("enforce_bcs");
    CH_XD::LevelData<CH_XD::FArrayBox>& phi = const_cast<CH_XD::LevelData<CH_XD::FArrayBox >& >(a_phi);
    phi.exchange(m_exchangeCopier);
  }

  std::shared_ptr<SGMultigridLevel> m_coarser;

private:
  void defineStencils();
  void defineCoarserObjects();
  void getMultiColors();

  /// weak construction introduces unnecessary complications
  SGMultigridLevel()
  {
  }

  bool   m_useDenseStencil;
  bool   m_hasCoarser; //not strictly necessary--could check the pointer
  double m_alpha;     
  double m_beta;       
  double m_dx;        
  CH_XD::DisjointBoxLayout     m_grids;       
  CH_XD::Copier m_exchangeCopier;
  double m_lambda; //relaxation coeff
  array<Point,MG_NUM_COLORS>  s_colors;
  Stencil<double> m_negoperator; //always need -Lap(phi)
  Stencil<double> m_restrict;

  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_updateOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpRhs;
  array<Stencil<double>, MG_NUM_COLORS> m_prolong;

  CH_XD::LevelData<CH_XD::FArrayBox > m_resid;
  CH_XD::LevelData<CH_XD::FArrayBox > m_residC;
  CH_XD::LevelData<CH_XD::FArrayBox > m_deltaC; 
  // for flop counts
  unsigned long long int numPointsThisProc(const CH_XD::DisjointBoxLayout& a_dbl) const
  {
    unsigned long long int retval = a_dbl.numPointsThisProc();
    return retval;
  }
  unsigned long long int m_numPointsThisProc; // for flop counts
};


///class that outsiders actually call
class SGMultigrid
{
public:
  static int  s_numSmoothDown  ;
  static int  s_numSmoothUp    ;
  static int  s_numSmoothBottom;

  /// 
  SGMultigrid(const double & a_alpha,
              const double & a_beta,
              const double & a_dx,
              const CH_XD::DisjointBoxLayout & a_grids,
              bool a_useDenseStencil);

  /// 
  void residual(CH_XD::LevelData<CH_XD::FArrayBox >       & a_res,
                const CH_XD::LevelData<CH_XD::FArrayBox > & a_phi,
                const CH_XD::LevelData<CH_XD::FArrayBox > & a_rhs);

  /// 
  void applyOp(CH_XD::LevelData<CH_XD::FArrayBox >       & a_lph,
               const CH_XD::LevelData<CH_XD::FArrayBox > & a_phi);

  /// Multigrid v-cycle.
  void vCycle(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phi,
              const CH_XD::LevelData<CH_XD::FArrayBox >& a_rhs);

  /// average down residual to next coarser level.
  void restrictResidual(CH_XD::LevelData<CH_XD::FArrayBox >      & a_resc,
                        const CH_XD::LevelData<CH_XD::FArrayBox >& a_res)
  {
    m_finest->m_coarser->restrictResidual(a_resc, a_res);
  }

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phiFine,
                        const CH_XD::LevelData<CH_XD::FArrayBox >& a_deltaCoarse)
  {
    m_finest->m_coarser->prolongIncrement(a_phiFine, a_deltaCoarse);
  }

  /// 
  void relax(CH_XD::LevelData<CH_XD::FArrayBox >      & a_phi,
             const CH_XD::LevelData<CH_XD::FArrayBox >& a_rhs)
  {
    m_finest->relax(a_phi, a_rhs);
  }
  ///
  void
  enforceBoundaryConditions(CH_XD::LevelData<CH_XD::FArrayBox >& a_phi)
  {
    m_finest->enforceBoundaryConditions(a_phi);
  }
private:
  /// weak construction introduces unnecessary complications
  SGMultigrid()
  {
  }

  std::shared_ptr<SGMultigridLevel> m_finest;
};
#endif
