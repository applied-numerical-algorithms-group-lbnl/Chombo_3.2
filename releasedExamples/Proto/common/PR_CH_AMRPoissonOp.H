
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

namespace PrChCommon
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
    pr_lbd1       m_resid;
    pr_lbd1       m_residC;
    pr_lbd1       m_deltaC; 
    pr_lbd1       m_phiC;
    pr_lbd1       m_scratch;
    
    pr_pt         m_ghost;
    double m_lambda; //relaxation coeff
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

    struct
    bc_eval_t
    {
      pr_box  m_applyBox;
      pr_sten m_stencil;
    };

    ///pick stencil for boundary conditions
    pr_sten
    getGhostFillStencil(unsigned int    a_idir,
                        ch_sd::LoHiSide a_sd,
                        string          a_bcname,
                        unsigned int a_ghostLayer) const //which layer of ghost cells are we filling
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
        ch_mayday::Error("PR_CH_AMRPoissonOp::getGhostFillStencil: unknown bc string");
      }
          
      pr_sten retval = (coeff)*pr_shift(pr_pt::Basis(a_idir, shiftval));
      
      return retval;
    }    
    
    // 
    int ebp_index(int a_idir, ch_sd::LoHiSide a_side)  const
    {
      return a_idir + a_side*DIM;
    }
    
    ///returns empty vector if there are no bcs to apply
    vector<bc_eval_t>
    getBCMetaData(pr_box         a_valid) const
    {
      vector<bc_eval_t> retval;
      pr_box dombx = m_dblMedi.domain().domainBox();
      pr_pt lobox = a_valid.low();
      pr_pt hibox = a_valid.high();
      pr_pt lodom = dombx.low();
      pr_pt hidom = dombx.high();
      for(unsigned int idir = 0; idir < DIM; idir++)
      {
        pr_box srcbx = a_valid.grow(pr_pt::Basis(idir));
        
        if(lobox[idir] == lodom[idir])
        {
          unsigned int  index = ebp_index(idir, ch_sd::Lo);
          string       bcname = m_bcNames[index];
          pr_sten   dombcsten = getGhostFillStencil(idir, ch_sd::Lo, bcname, 0);

          pr_pt loapp  = srcbx.low();
          pr_pt hiapp  = srcbx.high();
          loapp[idir]  = lodom[idir]-1;
          hiapp[idir]  = lodom[idir]-1;

          //for numbers of ghost cells > 1, need to do more shifting.   See Chombo4/src/EBProto/implem/Proto_SecondOrderCellImplem.H for an example
          pr_box appBox0(loapp, hiapp);  
          bc_eval_t entry;
          entry.m_applyBox = appBox0;
          entry.m_stencil  = dombcsten;
          retval.push_back(entry);
        }
        if(hibox[idir] == hidom[idir])
        {
          unsigned int index = ebp_index(idir, ch_sd::Hi);
          string      bcname = m_bcNames[index];
          pr_sten  dombcsten = getGhostFillStencil(idir, ch_sd::Hi, bcname, 0);

          pr_pt loapp = srcbx.low();
          pr_pt hiapp = srcbx.high();
          loapp[idir] = hidom[idir]+1;
          hiapp[idir] = hidom[idir]+1;
          //for numbers of ghost cells > 1, need to do more shifting.   See Chombo4/src/EBProto/implem/Proto_SecondOrderCellImplem.H for an example
          pr_box appBox0(loapp, hiapp);
          bc_eval_t entry;
          entry.m_applyBox = appBox0;
          entry.m_stencil  = dombcsten;
          retval.push_back(entry);
        }
      }
      return retval;

    }
    //fill one ghost cell outside the domain
    inline void
    applyDomainBoundaryConditions(pr_bd1          & a_src,
                                  unsigned int      a_idir,
                                  ch_sd::LoHiSide   a_sd,
                                  string            a_bcname,
                                  pr_di             a_dit)
    {
      int isign = sign(a_sd);
      auto bx = a_src.box();
      auto probdom = m_dblMedi.domain();
      auto probbox = probdom.domainBox();
      if(!probbox.contains(bx))
      {
        pr_box valid = m_dblMedi[a_dit];

        std::vector<bc_eval_t> bcmeta = getBCMetaData(valid);
        for(int ibc = 0; ibc < bcmeta.size(); ibc++)
        {
          const pr_box &  bcsrcbox  = bcmeta[ibc].m_applyBox;
          const pr_sten&  bcstencil = bcmeta[ibc].m_stencil;

          bcstencil.apply(a_src. a_src, bcsrcbox, true, 1.0); //true is to make this copy not increment
        }
      }
    }
    
    virtual void applyDomainBoundaryConditions (pr_lbd1 & a_phi)
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
    }
    virtual void homogeneousCFInterp(pr_lbd1& a_phi) 
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
    virtual void residual(  pr_lbd1& a_lhs, const pr_lbd1& a_phi, const pr_lbd1& a_rhs, bool a_homogeneous = false)
    {
      a_lhs.setVal(0.);
      applyOp(a_lhs, a_phi, a_homogeneous);
      pr_lbd1& castrhs = (pr_lbd1 &) a_rhs;
      a_lhs.increment(castrhs, -1.);
    }

    ///
    virtual void preCond(   pr_lbd1& a_cor, const pr_lbd1& a_residual)
    {
    
      relax(a_cor, a_residual, 27);
    }

    ///
    virtual void applyOp(   pr_lbd1& a_lhs, const pr_lbd1& a_phi, bool a_homogeneous = false)
    {
      //fill ghost cells
      pr_lbd1& castphi = (pr_lbd1 &)(a_phi);
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
    virtual void create(    pr_lbd1& a_lhs, const pr_lbd1& a_rhs)
    {
      a_lhs.define(a_rhs.disjointBoxLayout(), a_rhs.ghostVect());
    }

    ///
    virtual void assign(    pr_lbd1& a_lhs, const pr_lbd1& a_rhs)
    {
      a_lhs.setVal(0.);
      pr_lbd1 & castrhs = (pr_lbd1& )(a_rhs);
      a_lhs.increment(castrhs);
    }

    ///
    virtual double dotProduct(const pr_lbd1& a_1, const pr_lbd1& a_2)
    {
      pr_lbd1& dotprod = m_scratch;
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
    virtual void incr  (    pr_lbd1& a_lhs, const pr_lbd1& a_x, double a_scale)
    {
      pr_lbd1& increm = m_scratch;
      for(auto dit = a_lhs.begin(); dit != a_lhs.end(); ++dit)
      {
        increm[*dit].setVal(0);
        increm[*dit] += a_x[*dit];
        increm[*dit] *= a_scale;
        a_lhs[ *dit]  += increm[*dit];
      }
    }

    ///
    virtual void axby(      pr_lbd1& a_lhs, const pr_lbd1& a_x, const pr_lbd1& a_y, double a_a, double a_b)
    {
      a_lhs.setVal(0);
      incr(a_lhs,a_x, a_a);
      incr(a_lhs,a_y, a_b);
    }
  
    virtual void scale(     pr_lbd1& a_lhs,   const double& a_scale)
    {
      a_lhs *= a_scale;
    }
  
    virtual double norm(const pr_lbd1& a_rhs, int a_ord)
    {
      return a_rhs.absMax();
    }
  
    virtual void setToZero(pr_lbd1& a_lhs)
    {
      a_lhs.setVal(0.);
    }

///MGLevelOp (from MultiGrid.H):
    std::vector<pr_sten>  m_regStencilBCS[2*DIM];

    ///
    virtual void createCoarser(pr_lbd1& a_coarse, const pr_lbd1& a_fine, bool ghosted)
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
    virtual void relax(pr_lbd1& a_phi, const pr_lbd1& a_rhs, int a_maxiter);
  
    ///
    virtual void restrictResidual(pr_lbd1& a_resCoarse, const pr_lbd1& a_phiFine, const pr_lbd1& a_rhsFine);

  
    virtual void prolongIncrement(pr_lbd1& a_phiThisLevel, const pr_lbd1& a_correctCoarse);

  
    virtual void AMRResidual(pr_lbd1& a_residual, const pr_lbd1& a_phiFine, const pr_lbd1& a_phi,
                             const pr_lbd1& a_phiCoarse, const pr_lbd1& a_rhs,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<pr_lbd1>*  a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRResidualNF(pr_lbd1& a_residual, const pr_lbd1& a_phi, const pr_lbd1& a_phiCoarse,
                               const pr_lbd1& a_rhs, bool a_homogeneousBC)
    {
      MayDay::Error("not doing AMR yet");
    }
    

    virtual void AMRResidualNC(pr_lbd1& a_residual, const pr_lbd1& a_phiFine, const pr_lbd1& a_phi,
                               const pr_lbd1& a_rhs, bool a_homogeneousBC,
                               AMRLevelOp<pr_lbd1>* a_finerOp) 
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMROperator(pr_lbd1& a_LofPhi,
                             const pr_lbd1& a_phiFine, const pr_lbd1& a_phi,
                             const pr_lbd1& a_phiCoarse,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<pr_lbd1>*  a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }


    virtual void AMROperatorNF(pr_lbd1& a_LofPhi,
                               const pr_lbd1& a_phi,
                               const pr_lbd1& a_phiCoarse,
                               bool a_homogeneousBC)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMROperatorNC(pr_lbd1& a_LofPhi,
                               const pr_lbd1& a_phiFine,
                               const pr_lbd1& a_phi,
                               bool a_homogeneousBC,
                               AMRLevelOp<pr_lbd1>* a_finerOp)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRRestrict(pr_lbd1& a_resCoarse, const pr_lbd1& a_residual, const pr_lbd1& a_correction,
                             const pr_lbd1& a_coarseCorrection, bool a_skip_res)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void AMRProlong(pr_lbd1& a_correction, const pr_lbd1& a_coarseCorrection)
    {
      MayDay::Error("not doing AMR yet");
    }


    virtual void AMRUpdateResidual(pr_lbd1& a_residual, const pr_lbd1& a_correction,
                                   const pr_lbd1& a_coarseCorrection)
    {
      MayDay::Error("not doing AMR yet");
    }

    virtual void createCoarsened(pr_lbd1&       a_lhs,
                                 const pr_lbd1    & a_rhs,
                                 const int&     a_refRat)

    {
      MayDay::Error("not doing AMR yet");
    }

    /** going to try to do this with strong construction --placeholder
     */
    PR_CH_AMRPoissonOp(const pr_probdom& a_indexSpace)
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
    virtual AMRLevelOp<Proto::LevelBoxData<double> >* AMRnewOp(const ch_probdom& a_indexSpace)
    {
      MayDay::Error("not implemented");
      return NULL;
    }

    MGLevelOp<Proto::LevelBoxData<double> >* MGnewOp(const ch_probdom & a_indexSpace,
                                                     int           a_depth,
                                                     bool          a_homoOnly)
    {
      MayDay::Error("not implemented");
      return NULL;
    }

    ///
    /**
       return refinement ratio  to next finer level.
    */
    virtual int refToFiner(const ch_probdom& a_indexSpace) const
    {
      MayDay::Error("not implemented");
      return 2;
    }

  private:
    //I have decided that I am team strong construction
    PR_CH_AMRPoissonOpFactory();
    
    pr_amrdata m_data;
    double     m_dxLevel0;
  };
}
#endif
