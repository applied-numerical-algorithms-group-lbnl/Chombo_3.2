#include "AMRMultiGrid.H"
#include "Proto.H"

#ifndef ___AMRMultiGrid_Proto_Interface___
#define ___AMRMultiGrid_Proto_Interface___

///
/**
   Namespace PrCh_AMR_Elliptic

   Infrastructurue for using AMRMultiGrid as a solver front end and a
   proto back end.  I intend to rewrite in this fashion all of the
   some of the non-EB elliptic operators in Chombo3.2, including
   Poisson (constant coefficient and variable), ViscousTensor and
   Resistivity (variable coefficient for both).  I ignore all the
   use_homogeneous flags in the base classes because all domain bcs
   here should be homogeneous.  If you have inhomogeneous boundary
   conditions, slip them into the RHS.  Your data holder needs to be
   able to give its domain.  

   Brian van Straalen's wildly successful original templated
   elliptic solver to do the Martin-Cartwright algorithm has become
   a very large interface.  The base class hierarchy
   LevelOp/MGLevelOp/AMRLevelOp has (as of October 2023) 83 virtual
   functions.  Production codes can get large as more and more
   aggressive steps are taken to improve performance.  The
   AMRMultiGrid_Proto_Interface classes are meant to

   (1) Simplify using AMRMultiGrid.  
   (2) Provide the Proto backend for people who want to run on the device.

   The class structure reflects those priorities. 
   (1) Base_Simplified_AMRLevelOp provides reasonable defaults for many of these functions (independent of the proto interface).
   (2) Base_Proto_AMRLevelOp uses the proto interface to further winnow away at 83.
     
   In the end, all operators have to implement 7 functions.
**/
namespace PrCh_AMR_Elliptic
{
  typedef Chombo::Box  ch_box;

  ///
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
  class Base_Simplified_AMRLevelOp : public Chombo::AMRLevelOp<level_data_type>
  {
  public:

    typedef Chombo::Box           ch_box;
    typedef Chombo::ProblemDomain ch_dom;

    ///This class has no data.  
    Base_Simplified_AMRLevelOp():AMRLevelOp<>() {}

    ///This class has no data
    virtual ~Base_Simplified_AMRLevelOp(){}

    
    /// output = aco*xval + bco*yval
    /**
       Standard function that I use for all arithmetic.  Because this is not fortran,
       it does not matter if any of the arguments point to overlapping things.
     **/
    virtual void 
    axby(level_data_type         & a_output,
         const level_data_type   & a_xval,
         const level_data_type   & a_yval,
         const double            & a_bco,
         const double            & a_aco) = 0;
         
    ///a_phiThisLevel += I[2h->h](a_correctCoarse)
    virtual void prolongIncrement(T& a_phiThisLevel, const T& a_correctCoarse) = 0;


    ///Default implementations independent of Proto go here.
    
    /// a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
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
    virtual void
    AMROperator(level_data_type      & a_lphMedi,
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
     Base_Proto_AMRLevelOp
     All the proto-dependent default implementations go here.
     The distinction is probably academic but it does make for  smaller classes.

     This layer deals with proto<---->chombo conversions and will deal
     only with the purely virtual functions of
     Base_Simplified_AMRLevelOp.  This layer leans heavily upon the
     proto interface.  Only operator-dependent stuff is left out.

     The operator class has to implement the following:
     (1) applyOperator
     (2) relax 
     (3) fillGhostCellsOutsideDomain
     (4) refluxFineCorrection
     (5) refToCoarser.

     Some of these could be reasonably hardwired (force gsrb, first
     order prolongation, etc.), If this gets reused enough, we can
     do that. But for now, five shall be the number (and 5 << 83).  
     (4) and (5) can be delayed because they are amr things.

     An important design point here is that both base classes have no
     data.  I prefer base classes to just be base classes with all
     data living in the derived classs.  There are many default
     funcion implementations here.  All will have timers.  If, for a
     particular operator, one of these functions become odious, that operator
     can supercede it because all these functions are virtual,
  **/
  template<int ncomp>
  class Base_Proto_AMRLevelOp: public Base_Simplified_AMRLevelOp<Proto::LevelBoxData<ncomp> >
  {
  public:
    typedef Proto::LevelBoxData<double, ncomp> pr_lbd;
    
    ///this holds data that is common to all of these things.
    typedef Chombo::Box           ch_box;
    typedef Chombo::ProblemDomain ch_dom;

    ///This class has no data.  
    Base_Proto_AMRLevelOp():AMRLevelOp<>() {}

    ///This class has no data
    virtual ~Base_Proto_AMRLevelOp(){}

    ///
    virtual ch_dom
    getDomainFromData(const pr_lbd & a_data)
    {
      CH_TIME("Base_Proto_AMRLevelOp::getDomainFromData");
      auto   dbl_pr = a_data.layout();
      pr_dom dom_pr = dbl.domain();
      ch_dom retval = ProtoCh::getProblemDomain(dom_pr);
    }      
    
    ///exchange
    virtual void
    fillGhostDataOverFineFineInterface(pr_lbd       &  a_phiFine)
    {
      CH_TIME("Base_Proto_AMRLevelOp::fillGhostDataOverFineFineInterface");
      a_phiFine.exchange();
    }

    ///these two can be on different layouts
    virtual void
    copyTo(pr_lbd       &  a_dst,
           const pr_lbd &  a_src)
    {
      CH_TIME("Base_Proto_AMRLevelOp::copyTo");
      a_src.copyTo(a_dst);
    }
    
    ///
    virtual
    void createCoarser(pr_lbd                   & a_coar,
                       const pr_lbd             & a_fine,
                       bool ghosted)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      auto   dbl_fine = a_data.layout();
      auto   dbl_coar = dbl_fine.coarsen(2);
      pr_pt ghost = a_fine.ghost();
      a_coar.define(dbl_coar, ghost);
    }
    
    /// Cheers.
    virtual double  norm(const pr_lbd & a_rhs, int a_ord)
    {
      CH_TIME("Base_Proto_AMRLevelOp::norm");
      return a_rhs.absMax();
    }
    
    ///From EBLevelDataOps (credit to NASA's Mike Barad)
    double parallelSum(const double& a_value)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      // Find the sum of all a_value's
      Real sum = a_value;
#ifdef CH_MPI
      Real sendBuf = a_value;
      int result = MPI_Allreduce(&sendBuf, &sum, 1, MPI_CH_REAL,MPI_SUM, Chombo_MPI::comm);

      if (result != MPI_SUCCESS)
      {
        MayDay::Error("Base_Proto_AMRLevelOp::parallelSum communication error");
      }
#endif
      return sum;
    }

    ///
    /**
       These two have to on the same layout for this to work.
       Fortunately, I think that is forced by the calling function.
    **/
    virtual double
    dotProduct(const pr_lbd & a_fir,
               const pr_lbd & a_sec)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarser");
      auto dbl = a_fir.layout();
      auto dit = a_fir.begin();
      double localSum = 0.;
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto& firfab = a_fir[dit[ibox]];
        const auto& secfab = a_sec[dit[ibox]];
        auto multfab = firfab*secfab;
        localSum += multfab.sum();
      }
      double mpisum = parallelSum(localSum);
      return mpisum;
    }
    
    ///
    virtual void
    create(pr_lbd       & a_lhs,
           const pr_lbd & a_rhs)
    {
      CH_TIME("Base_Proto_AMRLevelOp::create");
      auto grids = a_rhs.layout();
      auto ghost = a_rhs.ghost();
      a_lhs.define(grids, ghost);
    }
    
    ///
    virtual void createCoarsened(pr_lbd       & a_phiCoar,
                                 const pr_lbd & a_phiFine,
                                 const int    & a_refRat)
    {
      CH_TIME("Base_Proto_AMRLevelOp::createCoarsened (the AMR one)");
      auto gridsFine = a_phiFine.layout();
      auto ghostFine = a_gridsFine.ghost();
      auto gridsCoar = gridsFine.coarsen(    Point::Ones(a_refRat));
      a_phiCoar.define(gridsCoar, ghostFine);
    }
    
    ///a_phiThisLevel += I[2h->h](a_correctCoarse).  factor of 2.  indifferent to layout.
    virtual void
    prolongIncrement(pr_lbd              & a_phiMedi,
                     const pr_lbd        & a_corCoar)
    {
      //borrowed from Phil's levelmultigrid example in proto/example
      pr_lbd localCor;
      createCoarser(   localCor, a_phiMedi);
      a_corCoar.copyTo(localCoar);
      pr_interp_sten interp_stencil = pr_interp_sten::Constant(Point::Ones(2));
      auto grids =  a_phiMedi.layout();
      auto dit   =  a_phiMedi.begin();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto & corfab =   localCor[dit[ibox]];
        const auto & phifab =  a_phiMedi[dit[ibox]];
      
        phifab += interp_stencil(corfab);
      }
    }

    
    ///a_resCoarse = average(a_resfine) on processor. 
    virtual void
    restriction(pr_lbd & a_resCoar,
                pr_lbd & a_resFine)
    {
      ///taken from Phil's levelmultigrid example
      auto A = Stencil<double>::AvgDown(2);
      auto grids =  a_phiMedi.layout();
      auto dit   =  a_phiMedi.begin();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        auto & finefab =   a_resFine[dit[ibox]];
        auto & coarfab =   a_resCoar[dit[ibox]];
        coarfab.setVal(0.);
        coarfab |= A(finefab);
      }
    }

    ///
    virtual void setToZero(level_data_type & a_lhs)
    {
      a_lhs.setToZero();
    }

    /// amr thing
    virtual void
    fillGhostDataOverCoarseFineInterface(pr_lbd       &  a_phiFine, 
                                         const pr_lbd &  a_phiCoar)
    {
      CH_TIME("Base_Proto_AMRLevelOp::fillGhostDataOverCoarseFineInterface");
      ch_mayday::Error("I think we can leave this until we are doing AMR");
    }
    

  };
} //end namespace PrCh_AMR_Elliptic

#endif
