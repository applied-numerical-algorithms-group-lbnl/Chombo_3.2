
#ifndef _AMRPOISSONOP_H_
#define _AMRPOISSONOP_H_


#include "AMRMultiGrid.H"
#include "Proto.H"
#include "Proto_AMRData.H"
#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

namespace Chombo
{
///
/**
   Operator for solving (alpha I + beta*Laplacian)(phi) = rho
   over an AMR hierarchy.
*/
  class PR_CH_AMRPoissonOp : public Chombo::AMRLevelOp<Proto::LevelBoxData<double> >
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
    virtual double dotProduct(const T& a_1, const T& a_2)     = 0;
    virtual void incr  (    T& a_lhs, const T& a_x, double a_scale) = 0;
    virtual void axby(      T& a_lhs, const T& a_x, const T& a_y, double a_a, double a_b) = 0;
    virtual void scale(     T& a_lhs, const double& a_scale)  = 0;
    virtual double norm(const T& a_rhs, int a_ord) = 0;
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

#include "PR_CH_TypeDefs.H"
    /// strong construction only
    PR_CH_AMRPoissonOp();


    //sent in directly or implied
    unsigned int m_refCoar;
    unsigned int m_refFine;

    bool   m_hasCoar; 
    bool   m_hasFine;
    
    double m_alpha;     
    double m_beta;       
    double m_dx;        

    pr_dbl m_dblMedi;
    pr_dbl m_dblFine;
    pr_dbl m_dblCoar;
    
    ///internal
    pr_lbd       m_resid;
    pr_lbd       m_residC;
    pr_lbd       m_deltaC; 
    pr_lbd       m_phiC;
    pr_pt        m_ghost;
    double m_lambda; //relaxation coeff
    pr_lbd m_scratch;
    array<Proto::Point,MG_NUM_COLORS>  s_colors;
    
    pr_sten m_negoperator; //always need -Lap(phi)
    pr_sten m_restrict;

    array<pr_sten, MG_NUM_COLORS> m_relaxOpPhi;
    array<pr_sten, MG_NUM_COLORS> m_updateOpPhi;
    array<pr_sten, MG_NUM_COLORS> m_relaxOpRhs;
    array<pr_sten, MG_NUM_COLORS> m_prolong;

    array<std::string, 2*DIM> m_bcNames;
    

  public:

    ///called by the factory class
    PR_CH_AMRPoissonOp(pr_dbl& a_dblMedi,  
                       double a_alpha,
                       double a_beta,
                       double a_dx,
                       pr_pt  a_ghost,
                       array<std::string, 2*DIM> a_domainBC,
                       std::pair<pr_dbl, int>* a_fine,   // = NULL if not there, int is the refinement ratio
                       std::pair<pr_dbl, int>* a_coar)   // = NULL if not there, int is the refinement ratio
    {
      m_alpha = a_alpha;
      m_beta  = a_beta;
      m_hasFine = (a_fine != NULL);
      m_hasCoar = (a_coar != NULL);

      if(m_hasFine)
      {
        m_dblFine = a_fine->first;
        m_refFine = a_fine->second;
      }
      if(m_hasCoar)
      {
        m_dblCoar = a_fine->first;
        m_refCoar = a_fine->second;
        m_residC.define(m_dblCoar,  a_ghost);
        m_deltaC.define(m_dblCoar,  a_ghost);
      }
      m_dx = a_dx;
      double diagvalu = -2.*DIM/(a_dx*a_dx);
      double realdiag = a_alpha + a_beta*diagvalu;
      m_lambda = 1./realdiag;
      m_ghost = a_ghost;
      
    }


    Box
    getBCBox(Bx a_valid,
             unsigned int idir,
             Side::LoHiSide a_sd)
    {
      Point loapp, hiapp
    }
    Stencil<double>
    getGhostFillStencil(unsigned int    a_idir,
                        Side::LoHiSide, a_sd,
                        string          a_bcname,
                        unsigned int a_ghostLayer) //which layer of ghost cells are we filling
    {
      int isign = sign(a_sd);
      int shiftSign = -isign;
      int shiftval  = (a_ghostLayer+1)*shiftSign;
      
      Real coeff = 0;
      if(a_bcname == string("Dirichlet"))
      {
        coeff = -1;//Dirichlet is reflect odd
      }
      else if(a_bcname == string("Neumann"))
      {
        coeff = 1; //Neumann is reflect even
      }
      else
      {
        MayDay::Error("PR_CH_AMRPoissonOp::getGhostFillStencil: unknown bc string");
      }
          
      Stencil<double> retval = (coeff)*Shift(Point::Basis(a_idir, shiftval));
      
      return retval;

             
    }
    //fill one ghost cell outside the domain
    inline void
    applyDomainBoundaryConditions(BoxData<Real, 1>& a_src,
                                  unsigned int      a_idir,
                                  Side::LoHiSide    a_sd,
                                  string            a_bcname,
                                  Proto::DataIndex  a_dit)
    {
      int isign = sign(a_sd);
      auto bx = a_src.box();
      auto probdom = m_dblMedi.problemDomain();
      
      if(!probdom.contains(bx))
      {
        Box valid = m_dblMedi[a_dit];

        Stencil<double> dombcsten = getGhostFillStencil(a_idir, a_sd, a_bcname, 0);
        Bx bcbox = getBCBox(valid, a_idir, a_sd);
        for(int ibc = 0; ibc < 2*DIM; ibc++)
        {
          for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
          {
            Box bcsrcbox = m_BCApplyBoxes[ibc][ivec];
            m_regStencilBCS[ibc][ivec].apply(castSrc, castSrc,
                                             bcsrcbox, initToZero, 1.0);
          }
        }
      }
    }
    
    virtual void applyDomainBoundaryConditions (pr_lbd& a_phi)
    {
      CH_TIME("applyDomainBoundaryConditions");
      for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          for(unsigned int idir = 0; idir < DIM; idir++)
          {
            unsigned int index = ebp_index(idir, sit());
            string bcname = m_bcNames[index];
            applyDomainBoundaryConditions(a_phi[*dit], idir, sit(), bcname, *dit);
          }
      }

    }
    virtual void homogeneousCFInterp(pr_lbd& a_phi) 
    {
      if(m_hasCoar)
      {
        MayDay::Error("cf interp not there yet");
        m_phiC.setVal(0);
//        Proto::InterpStencil<double> interpSten = Proto::InterpStencil<double>::Build(2, m_refToCoar);
//        a_phi.interpBoundaries(m_phiC, interpSten);
      }
      a_phi.exchange();
    }

    ///
    virtual void residual(  pr_lbd& a_lhs, const pr_lbd& a_phi, const pr_lbd& a_rhs, bool a_homogeneous = false)
    {
      a_lhs.setVal(0.);
      applyOp(a_lhs, a_phi, a_homogeneous);
      pr_lbd& castrhs = (pr_lbd &) a_rhs;
      a_lhs.increment(castrhs, -1.);
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
      pr_lbd& castphi = (pr_lbd &)(a_phi);
      homogeneousCFInterp(castphi);
      applyDomainBoundaryConditions(castphi);
      auto dbl = a_lhs.layout();
      //apply the actual operator
      for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
      {
        auto index = *dit;
        auto valid = dbl[index];
        a_lhs[index] |= m_negoperator(a_phi[index], valid);
        a_lhs[index] *= -1.0;
      }
    }

    ///
    virtual void create(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
    {
      a_lhs.define(a_rhs.disjointBoxLayout(), a_rhs.ghostVect());
    }

    ///
    virtual void assign(    pr_lbd& a_lhs, const pr_lbd& a_rhs)
    {
      a_lhs.setVal(0.);
      pr_lbd & castrhs = (pr_lbd& )(a_rhs);
      a_lhs.increment(castrhs);
    }

    ///
    virtual double dotProduct(const pr_lbd& a_1, const pr_lbd& a_2)
    {
      pr_lbd& dotprod = m_scratch;
      double sumdotprod = 0;
      auto dbl = a_1.disjointBoxLayout();
      for(auto dit = a_1.begin(); dit != a_1.end(); ++dit)
      {
        auto valid = dbl[*dit];
        dotprod[*dit].setVal(0.);
        dotprod[*dit] += a_1[*dit];
        dotprod[*dit] *= a_2[*dit];
        sumdotprod += dotprod[*dit].sum();
      }

#ifdef CH_MPI
      double sendBuf = sumdotprod;
      int result = MPI_Allreduce(&sendBuf, &sumdotprod, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

      if (result != MPI_SUCCESS)
      {
        MayDay::Error("Communication error in PR_CH_AMRPoissonOp::dotProduct");
      }
#endif
      return sumdotprod;
    }

    ///
    virtual void incr  (    pr_lbd& a_lhs, const pr_lbd& a_x, double a_scale)
    {
      pr_lbd& increm = m_scratch;
      for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
      {
        increm[*dit].setVal(0);
        increm[*dit] += a_x[*dit];
        increm[*dit] *= a_scale;
        a_lhs[ *dit]  += increm[*dit];
      }
    }

    ///
    virtual void axby(      pr_lbd& a_lhs, const pr_lbd& a_x, const pr_lbd& a_y, double a_a, double a_b)
    {
      a_lhs.setVal(0);
      incr(a_lhs,a_x, a_a);
      incr(a_lhs,a_y, a_b);
    }
  
    virtual void scale(     pr_lbd& a_lhs,   const double& a_scale)
    {
      a_lhs *= a_scale;
    }
  
    virtual double norm(const pr_lbd& a_rhs, int a_ord)
    {
      return a_rhs.absMax();
    }
  
    virtual void setToZero(pr_lbd& a_lhs)
    {
      a_lhs.setVal(0.);
    }

///MGLevelOp (from MultiGrid.H):
    vector<Stencil<double> >  m_regStencilBCS[2*DIM];

    ///
    virtual void createCoarser(pr_lbd& a_coarse, const pr_lbd& a_fine, bool ghosted)
    {
      auto dblFine = a_fine.disjointBoxLayout();
      auto ghost   = a_fine.ghostVect();
      auto dblCoar = dblFine.coarsen(Proto::Point::Ones(2));
      a_coarse.define(dblCoar, ghost);
    }

    ///
    virtual int refToCoarser()
    {
      return m_refCoar;
    }

    ///
    virtual void relax(pr_lbd& a_phi, const pr_lbd& a_rhs, int a_maxiter);
  
    ///
    virtual void restrictResidual(pr_lbd& a_resCoarse, const pr_lbd& a_phiFine, const pr_lbd& a_rhsFine);

  
    virtual void prolongIncrement(pr_lbd& a_phiThisLevel, const pr_lbd& a_correctCoarse);

  
    virtual void AMRResidual(pr_lbd& a_residual, const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                             const pr_lbd& a_phiCoarse, const pr_lbd& a_rhs,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<pr_lbd>*  a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRResidualNF(pr_lbd& a_residual, const pr_lbd& a_phi, const pr_lbd& a_phiCoarse,
                               const pr_lbd& a_rhs, bool a_homogeneousBC)
    {
      MayDay::Error("not doing AMR yet");
    }
    

    virtual void AMRResidualNC(pr_lbd& a_residual, const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                               const pr_lbd& a_rhs, bool a_homogeneousBC,
                               AMRLevelOp<pr_lbd>* a_finerOp) 
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMROperator(pr_lbd& a_LofPhi,
                             const pr_lbd& a_phiFine, const pr_lbd& a_phi,
                             const pr_lbd& a_phiCoarse,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<pr_lbd>*  a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }


    virtual void AMROperatorNF(pr_lbd& a_LofPhi,
                               const pr_lbd& a_phi,
                               const pr_lbd& a_phiCoarse,
                               bool a_homogeneousBC)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMROperatorNC(pr_lbd& a_LofPhi,
                               const pr_lbd& a_phiFine,
                               const pr_lbd& a_phi,
                               bool a_homogeneousBC,
                               AMRLevelOp<pr_lbd>* a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRRestrict(pr_lbd& a_resCoarse, const pr_lbd& a_residual, const pr_lbd& a_correction,
                             const pr_lbd& a_coarseCorrection, bool a_skip_res)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRProlong(pr_lbd& a_correction, const pr_lbd& a_coarseCorrection)
    {
      MayDay::Error("not doing AMR yet");
    }


    virtual void AMRUpdateResidual(pr_lbd& a_residual, const pr_lbd& a_correction,
                                   const pr_lbd& a_coarseCorrection)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void createCoarsened(pr_lbd&       a_lhs,
                                 const pr_lbd& a_rhs,
                                 const int&     a_refRat)

    {
      MayDay::Error("not doing AMR yet");
    }

    /** going to try to do this with strong construction --placeholder
     */
    PR_CH_AMRPoissonOp(const ProblemDomain& a_indexSpace)
    {
      MayDay::Error("Still need a constructor");
    }

    ///
    /**
     */
    virtual ~PR_CH_AMRPoissonOp()
    {
    }

  };

  class PR_CH_AMRPoissonOpFactory: public Chombo::AMRLevelOpFactory<Proto::LevelBoxData<double> > 
  {
//    PR_CH_AMRPoissonOp(pr_dbl& a_dblMedi,  
//                       double a_alpha,
//                       double a_beta,
//                       double a_dx,
//                       pr_pt  a_ghost,
//                       array<std::string, 2*DIM> a_domainBC,
//                       std::pair<pr_dbl, int>* a_fine,   // = NULL if not there, int is the refinement ratio
//                       std::pair<pr_dbl, int>* a_coar)   // = NULL if not there, int is the refinement ratio
  public:
#include "PR_CH_TypeDefs.H"
    typedef  Proto::AMRData<Real>   pr_amrdata;
    PR_CH_AMRPoissonOpFactory(pr_amrdata a_data, 
                              double     a_dxLevel0)
    {
      m_data     = a_data;
      m_dxLevel0 = a_dxLevel0;
    }


  public:
    virtual AMRLevelOp<Proto::LevelBoxData<double> >* AMRnewOp(const ProblemDomain& a_indexSpace)
    {
      MayDay::Error("not implemented");
      return NULL;
    }

    MGLevelOp<Proto::LevelBoxData<double> >* MGnewOp(const ProblemDomain& a_indexSpace,
                                                     int                  a_depth,
                                                     bool                 a_homoOnly)
    {
      MayDay::Error("not implemented");
      return NULL;
    }

  ///
  /**
     return refinement ratio  to next finer level.
   */
  virtual int refToFiner(const ProblemDomain& a_indexSpace) const =0;

  private:
    //I have decided that I am team strong construction
    PR_CH_AMRPoissonOpFactory();
    
    pr_amrdata m_data;
    double     m_dxLevel0;
  };
} //end namespace
#endif
