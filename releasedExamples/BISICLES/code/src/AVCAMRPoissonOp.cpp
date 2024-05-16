#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "AMRMultiGrid.H"
#include "Misc.H"
#include "AMRPoissonOpF_F.H"
#include "ViscousTensorOpF_F.H"
#include "ExtrapBCF_F.H"
#include "AVCAMRPoissonOp.H"
#include "AVCAMRPoissonOpF_F.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

void AVCAMRPoissonOp::residualI(LevelData<FArrayBox>&       a_lhs,
                               const LevelData<FArrayBox>& a_phi,
                               const LevelData<FArrayBox>& a_rhs,
                               bool                        a_homogeneous)
{
  CH_TIME("AVCAMRPoissonOp::residualI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("AVCAMRPoissonOp::residualIBC");

    for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }
  }

  phi.exchange(phi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit()];
      const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
      FORT_AVCCOMPUTERES1D
#elif CH_SPACEDIM == 2
      FORT_AVCCOMPUTERES2D
#elif CH_SPACEDIM == 3
      FORT_AVCCOMPUTERES3D
#else
      This_will_not_compile!
#endif
	                 (CHF_FRA1(a_lhs[dit],0),
                          CHF_CONST_FRA1(phi[dit],0),
                          CHF_CONST_FRA1(a_rhs[dit],0),
                          CHF_CONST_REAL(m_alpha),
                          CHF_CONST_FRA1((*m_aCoef)[dit],0),
                          CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                          CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                          CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                          CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                          This_will_not_compile!
#endif
                          CHF_BOX(region),
                          CHF_CONST_REAL(m_dx));

      CH_assert(a_lhs[dit].norm(0) < 1.0e+10);

    } // end loop over boxes

  int dbg=0;dbg++;

}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void AVCAMRPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::preCond");

  // diagonal term of this operator in:
  //
  //       alpha * a(i)
  //     + beta  * sum_over_dir (b(i-1/2*e_dir) + b(i+1/2*e_dir)) / (dx*dx)
  //
  // The inverse of this is our initial multiplier.

  int ncomp = a_phi.nComp();

  CH_assert(m_lambda.isDefined());
  CH_assert(a_rhs.nComp()    == ncomp);
  CH_assert(m_bCoef->nComp() == SpaceDim);

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // also need to average and sum face-centered bCoefs to cell-centers
      Box gridBox = a_rhs[dit].box();

      // approximate inverse
      a_phi[dit].copy(a_rhs[dit]);
      a_phi[dit].mult(m_lambda[dit], gridBox, 0, 0, ncomp);
    }

  relax(a_phi, a_rhs, 2);
}

void AVCAMRPoissonOp::applyOpI(LevelData<FArrayBox>&      a_lhs,
                             const LevelData<FArrayBox>& a_phi,
                             bool                        a_homogeneous )
{
  CH_TIME("AVCAMRPoissonOp::applyOpI");
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }

  applyOpNoBoundary(a_lhs, a_phi);
}

void AVCAMRPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&      a_lhs,
                                        const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("AVCAMRPoissonOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  phi.exchange(phi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit()];
      const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
      FORT_AVCCOMPUTEOP1D
#elif CH_SPACEDIM == 2
      FORT_AVCCOMPUTEOP2D
#elif CH_SPACEDIM == 3
      FORT_AVCCOMPUTEOP3D
#else
      This_will_not_compile!
#endif
	                (CHF_FRA1(a_lhs[dit],0),
	                 CHF_CONST_FRA1(phi[dit],0),
                         CHF_CONST_REAL(m_alpha),
                         CHF_CONST_FRA1((*m_aCoef)[dit],0),
                         CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                         CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                         CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                         CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                         This_will_not_compile!
#endif
                         CHF_BOX(region),
                         CHF_CONST_REAL(m_dx));
    } // end loop over boxes
}

void AVCAMRPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                      LevelData<FArrayBox>&       a_phiFine,
                                      const LevelData<FArrayBox>& a_rhsFine)
{
  CH_TIME("AVCAMRPoissonOp::restrictResidual");

  homogeneousCFInterp(a_phiFine);
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi = a_phiFine[dit];
      m_bc(phi, dblFine[dit()], m_domain, m_dx, true);
    }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       phi = a_phiFine[dit];
      const FArrayBox& rhs = a_rhsFine[dit];
      FArrayBox&       res = a_resCoarse[dit];

      const FArrayBox& thisACoef = (*m_aCoef)[dit];
      const FluxBox&   thisBCoef = (*m_bCoef)[dit];

      Box region = dblFine.get(dit());
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, 2);

      res.setVal(0.0);

#if CH_SPACEDIM == 1
      FORT_RESTRICTRESAVC1D
#elif CH_SPACEDIM == 2
      FORT_RESTRICTRESAVC2D
#elif CH_SPACEDIM == 3
      FORT_RESTRICTRESAVC3D
#else
      This_will_not_compile!
#endif
	                  (CHF_FRA1_SHIFT(res, 0, civ),
                           CHF_CONST_FRA1_SHIFT(phi, 0, iv),
                           CHF_CONST_FRA1_SHIFT(rhs, 0, iv),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_FRA1_SHIFT(thisACoef, 0, iv),
                           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                           CHF_CONST_FRA_SHIFT(thisBCoef[0], iv),
#endif
#if CH_SPACEDIM >= 2
                           CHF_CONST_FRA_SHIFT(thisBCoef[1], iv),
#endif
#if CH_SPACEDIM >= 3
                           CHF_CONST_FRA_SHIFT(thisBCoef[2], iv),
#endif
#if CH_SPACEDIM >= 4
                           This_will_not_compile!
#endif
                           CHF_BOX_SHIFT(region, iv),
                           CHF_CONST_REAL(m_dx));
      CH_assert(res.norm(0) < 1.0e+10);
    }
}



void AVCAMRPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                      const Real& a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void AVCAMRPoissonOp::setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                               const RefCountedPtr<LevelData<FluxBox  > >& a_bCoef,
                               const Real&                                 a_alpha,
                               const Real&                                 a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void AVCAMRPoissonOp::resetLambda()
{
  if (m_lambdaNeedsResetting)
  {
    Real scale = 1.0 / (m_dx*m_dx);

    // Compute it box by box, point by point
    for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       lambdaFab = m_lambda[dit];
      const FArrayBox& aCoefFab  = (*m_aCoef)[dit];
      const FluxBox&   bCoefFab  = (*m_bCoef)[dit];
      const Box& curBox = lambdaFab.box();

      // Compute the diagonal term
      lambdaFab.copy(aCoefFab);
      lambdaFab.mult(m_alpha);
      CH_assert(lambdaFab.norm(1) < 1.0e+10);
      for (int dir = 0; dir < SpaceDim; dir++)
      {
        FORT_SUMAFACES(CHF_FRA1(lambdaFab,0),
            CHF_CONST_REAL(m_beta),
            CHF_CONST_FRA(bCoefFab[dir]),
            CHF_BOX(curBox),
            CHF_CONST_INT(dir),
            CHF_CONST_REAL(scale));
      }

      // Take its reciprocal
      lambdaFab.invert(1.0);
      CH_assert(lambdaFab.norm(1) < 1.0e+10);
    }

    // Lambda is reset.
    m_lambdaNeedsResetting = false;
  }
}

// Compute the reciprocal of the diagonal entry of the operator matrix
void AVCAMRPoissonOp::computeLambda()
{
  CH_TIME("AVCAMRPoissonOp::computeLambda");

  CH_assert(!m_lambda.isDefined());

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(),m_aCoef->nComp());
  resetLambda();
}


void AVCAMRPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                            const LevelData<FArrayBox>&        a_phi,
                            LevelData<FArrayBox>&              a_residual,
                            AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("AVCAMRPoissonOp::reflux");
  
  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_TIMER("AVCAMRPoissonOp::reflux::incrementCoarse", t2);
  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& coarfab   = a_phi[dit];
    const FluxBox&   coarBCoef = (*m_bCoef)[dit];
    const Box&       gridBox   = a_phi.getBoxes()[dit];

    if (m_levfluxreg.hasCF(dit()))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        FArrayBox coarflux;
        Box faceBox = surroundingNodes(gridBox, idir);

        getFlux(coarflux, coarfab, coarBCoef, faceBox, idir);

        Real scale = 1.0;
        m_levfluxreg.incrementCoarse(coarflux, scale,dit(),
            interv, interv, idir);
      }
    }
  }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  AVCAMRPoissonOp* finerAMRPOp = (AVCAMRPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);

#if CH_SPACEDIM == 2
  //slc needs to do this properly
  phiFineRef.exchange(a_phiFine.interval());
  for (DataIterator dit(phiFineRef.dataIterator());dit.ok();++dit)
    {
      Box box = phiFineRef.disjointBoxLayout()[dit];
      //box.grow(-1);
      FORT_EXTRAPCORNER2D(CHF_FRA(phiFineRef[dit]), CHF_BOX(box));
    }
#endif 


  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // slc -- it helps here though, at least with the rate.
  phiFineRef.exchange(a_phiFine.interval());
  IntVect phiGhost = phiFineRef.ghostVect();
  int ncomps = a_phiFine.nComp();

  CH_TIMER("AVCAMRPoissonOp::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab   = a_phiFine[ditf];
      const FluxBox&   fineBCoef = (*(finerAMRPOp->m_bCoef))[ditf];
      const Box&       gridbox   = dblFine.get(ditf());

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              if (m_levfluxreg.hasCF(ditf(), sit()))
                {
                  Side::LoHiSide hiorlo = sit();
                  Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

                  FArrayBox fineflux(fluxBox,ncomps);
                  getFlux(fineflux, phifFab, fineBCoef, fluxBox, idir,
                          m_refToFiner);

                  Real scale = 1.0;
                  m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                             interv, interv, idir, hiorlo);
                }
            }
        }
    }

  CH_STOP(t3);

  Real scale = 1.0/m_dx;
  m_levfluxreg.reflux(a_residual, scale);
}


void AVCAMRPoissonOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::levelGSRB");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();

  // loop over 3 x (DIM) color scheme
#if CH_SPACEDIM > 2
  for (int kcol = 0; kcol <= 2; kcol++)
    {
#endif
#if CH_SPACEDIM > 1
      for (int jcol = 0; jcol <= 2; jcol++)
	{
#endif
	  for (int icol = 0; icol <= 2; icol++)
	    {
	      CH_TIMERS("AVCAMRPoissonOp::levelGSRB::Compute");
	      
	      // fill in intersection of ghostcells and a_phi's boxes
	      {
		CH_TIME("AVCAMRPoissonOp::levelGSRB::homogeneousCFInterp");
		homogeneousCFInterp(a_phi);
	      }
	      
	      {
		CH_TIME("AVCAMRPoissonOp::levelGSRB::exchange");
		a_phi.exchange(a_phi.interval(), m_exchangeCopier);
	      }
	      
	      {
		CH_TIME("AVCAMRPoissonOp::levelGSRB::BCs");
		// now step through grids...
		for (dit.begin(); dit.ok(); ++dit)
		  {
		    // invoke physical BC's where necessary
		    m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
		  }
	      }
	      
	      for (dit.begin(); dit.ok(); ++dit)
		{
		  const Box& region = dbl.get(dit());
		  const FluxBox& thisBCoef  = (*m_bCoef)[dit];
		  CH_assert( m_lambda[dit].norm(1) < 1.0e+10);
		  
#if CH_SPACEDIM == 1
		  FORT_GSRBHELMHOLTZAVC1D
#elif CH_SPACEDIM == 2
		    FORT_GSRBHELMHOLTZAVC2D
#elif CH_SPACEDIM == 3
		    FORT_GSRBHELMHOLTZAVC3D
#else
		    This_will_not_compile!
#endif
		    (CHF_FRA1(a_phi[dit],0),
		     CHF_CONST_FRA1(a_rhs[dit],0),
		     CHF_BOX(region),
		     CHF_CONST_REAL(m_dx),
		     CHF_CONST_REAL(m_alpha),
		     CHF_CONST_FRA1((*m_aCoef)[dit],0),
		     CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
		     CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
		     CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
		     CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
		     This_will_not_compile!
#endif
		     CHF_CONST_FRA1(m_lambda[dit],0),
		     CHF_CONST_INT(icol)
#if CH_SPACEDIM > 1
		     ,CHF_CONST_INT(jcol)
#endif
#if CH_SPACEDIM > 2
		     ,CHF_CONST_INT(kcol)
#endif
		     );
		  CH_assert( a_phi[dit].norm(0) < 1.0e+10);
		} // end loop through grids
	    } // end loop through red-black (i)
#if CH_SPACEDIM > 1
	} // end loop through red-black (j)
#endif
#if CH_SPACEDIM > 2
    } // end loop through red-black (k)
#endif
}

void AVCAMRPoissonOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                     const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::levelMultiColor");
  MayDay::Abort("AVCAMRPoissonOp::levelMultiColor - Not implemented");
}

void AVCAMRPoissonOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::looseGSRB");
#if 1
  MayDay::Abort("AVCAMRPoissonOp::looseGSRB - Not implemented");
#else
  // This implementation converges at half the rate of "levelGSRB" in
  // multigrid solves
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();

  // fill in intersection of ghostcells and a_phi's boxes
  {
    CH_TIME("AVCAMRPoissonOp::looseGSRB::homogeneousCFInterp");
    homogeneousCFInterp(a_phi);
  }

  {
    CH_TIME("AVCAMRPoissonOp::looseGSRB::exchange");
    a_phi.exchange(a_phi.interval(), m_exchangeCopier);
  }

  // now step through grids...
  for (dit.begin(); dit.ok(); ++dit)
  {
    // invoke physical BC's where necessary
    {
      CH_TIME("AVCAMRPoissonOp::looseGSRB::BCs");
      m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
    }

    const Box& region = dbl.get(dit());
    const FluxBox& thisBCoef  = (*m_bCoef)[dit];

    int whichPass = 0;

#if CH_SPACEDIM == 1
    FORT_GSRBHELMHOLTZAVC1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZAVC2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZAVC3D
#else
    This_will_not_compile!
#endif
                          (CHF_FRA(a_phi[dit]),
                           CHF_CONST_FRA(a_rhs[dit]),
                           CHF_BOX(region),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_FRA((*m_aCoef)[dit]),
                           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                           CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                           CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                           CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                           This_will_not_compile!
#endif
                           CHF_CONST_FRA(m_lambda[dit]),
                           CHF_CONST_INT(whichPass));

    whichPass = 1;

#if CH_SPACEDIM == 1
    FORT_GSRBHELMHOLTZAVC1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZAVC2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZAVC3D
#else
    This_will_not_compile!
#endif
                          (CHF_FRA(a_phi[dit]),
                           CHF_CONST_FRA(a_rhs[dit]),
                           CHF_BOX(region),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_FRA((*m_aCoef)[dit]),
                           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                           CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                           CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                           CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                           This_will_not_compile!
#endif
                           CHF_CONST_FRA(m_lambda[dit]),
                           CHF_CONST_INT(whichPass));
    CH_assert(a_phi[dit] < 1.0e+10);

    
  } // end loop through grids
#endif
}

void AVCAMRPoissonOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::overlapGSRB");
  MayDay::Abort("AVCAMRPoissonOp::overlapGSRB - Not implemented");
}

void AVCAMRPoissonOp::levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                                   const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::levelGSRBLazy");
  MayDay::Abort("AVCAMRPoissonOp::levelGSRBLazy - Not implemented");
}

void AVCAMRPoissonOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AVCAMRPoissonOp::levelJacobi");

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  LevelData<FArrayBox> resid;
  create(resid, a_rhs);

  // Get the residual
  residual(resid,a_phi,a_rhs,true);

  // Multiply by the weights
  DataIterator dit = m_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    resid[dit].mult(m_lambda[dit]);
  }

  // Do the Jacobi relaxation
  incr(a_phi, resid, 0.5);
}

void AVCAMRPoissonOp::getFlux(FArrayBox&       a_flux,
                             const FArrayBox& a_data,
                             const FluxBox&   a_bCoef,
                             const Box&       a_facebox,
                             int              a_dir,
                             int              a_ref) const
{
  CH_TIME("AVCAMRPoissonOp::getFlux");
  //MayDay::Error("AVCAMRPoissonOp::getFlux not correct");
  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox& bCoefDir = a_bCoef[a_dir];

  // reality check for bCoef
  CH_assert(bCoefDir.box().contains(a_facebox));

  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real scale = m_beta * a_ref / m_dx;

  int tdir = (a_dir + 1)%SpaceDim;
 
  for ( bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      {
	//normal
	IntVect shiftiv = BASISV(a_dir);
	IntVect ivlo = iv - shiftiv;
	IntVect ivhi = iv;
	CH_assert(a_data.box().contains(ivlo));
	CH_assert(a_data.box().contains(ivhi));
	Real phihi = a_data(ivhi);
	Real philo = a_data(ivlo);
	Real gradphi = (phihi - philo ) * scale;
	a_flux(iv) = -bCoefDir(iv, a_dir) * gradphi;
      }

#if CH_SPACEDIM == 2
      {
	//transverse
	IntVect ivmp = iv - BASISV(a_dir) + BASISV(tdir);
	IntVect ivmm = iv - BASISV(a_dir) - BASISV(tdir);
	IntVect ivpp = iv + BASISV(tdir);
	IntVect ivpm = iv - BASISV(tdir);
	Real gradphi = (a_data(ivpp) + a_data(ivmp) - a_data(ivpm) - a_data(ivmm) ) * scale * 0.25;
	a_flux(iv) += -bCoefDir(iv, tdir) * gradphi;
      }
#endif

#if CH_SPACEDIM > 2
      MayDay::Error ("AVCAMRPoissonOp::getFlux no implemented for CH_SPACEDIM > 2");
#endif

        
    }
}

//-----------------------------------------------------------------------
void
AVCAMRPoissonOp::
setTime(Real a_time)
{
  // Jot down the time.
  m_time = a_time;

  // Interpolate the b coefficient data if necessary / possible. If
  // the B coefficient depends upon the solution, the operator is nonlinear
  // and the integrator must decide how to treat it.
  if (!m_bCoefInterpolator.isNull() &&
      !m_bCoefInterpolator->dependsUponSolution())
    m_bCoefInterpolator->interpolate(*m_bCoef, a_time);

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;

  // Set the time on the boundary holder.
  m_bc.setTime(a_time);

  // Notify our observers that the time has been set.
  // FIXME: Must implement response of multigrid operators!
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

// Factory
AVCAMRPoissonOpFactory::AVCAMRPoissonOpFactory()
{
  setDefaultValues();
}

//-----------------------------------------------------------------------
//  AMR Factory define function
void AVCAMRPoissonOpFactory::define(const ProblemDomain&                           a_coarseDomain,
                                   const Vector<DisjointBoxLayout>&               a_grids,
                                   const Vector<int>&                             a_refRatios,
                                   const Real&                                    a_coarsedx,
                                   BCHolder                                       a_bc,
                                   const Real&                                    a_alpha,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                                   const Real&                                    a_beta,
                                   Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef)
{
  CH_TIME("AVCAMRPoissonOpFactory::define");

  setDefaultValues();

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;

  m_dx.resize(a_grids.size());
  m_dx[0] = a_coarsedx;

  m_domains.resize(a_grids.size());
  m_domains[0] = a_coarseDomain;

  m_exchangeCopiers.resize(a_grids.size());
  m_exchangeCopiers[0].exchangeDefine(a_grids[0], IntVect::Unit);
  //m_exchangeCopiers[0].trimEdges(a_grids[0], IntVect::Unit);

  m_cfregion.resize(a_grids.size());
  m_cfregion[0].define(a_grids[0], m_domains[0]);

  for (int i = 1; i < a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1] / m_refRatios[i-1];

      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);

      m_exchangeCopiers[i].exchangeDefine(a_grids[i], IntVect::Unit);
      //m_exchangeCopiers[i].trimEdges(a_grids[i], IntVect::Unit);

      m_cfregion[i].define(a_grids[i], m_domains[i]);
    }

  m_alpha = a_alpha;
  m_aCoef = a_aCoef;

  m_beta  = a_beta;
  m_bCoef = a_bCoef;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// AMR Factory define function, with coefficient data allocated automagically
// for operators.
void
AVCAMRPoissonOpFactory::
define(const ProblemDomain& a_coarseDomain,
       const Vector<DisjointBoxLayout>& a_grids,
       const Vector<int>& a_refRatios,
       const Real& a_coarsedx,
       BCHolder a_bc,
       const IntVect& a_ghostVect)
{
  // This just allocates coefficient data, sets alpha = beta = 1, and calls
  // the other define() method.
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(a_grids.size());
  for (int i = 0; i < a_grids.size(); ++i)
  {
    aCoef[i] = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    bCoef[i] = RefCountedPtr<LevelData<FluxBox> >(
                 new LevelData<FluxBox>(a_grids[i], 1, a_ghostVect));

    // Initialize the a and b coefficients to 1 for starters.
    for (DataIterator dit = aCoef[i]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[i])[dit()].setVal(1.0);
      for (int idir = 0; idir < SpaceDim; ++idir)
        (*bCoef[i])[dit()][idir].setVal(1.0);
    }
  }
  Real alpha = 1.0, beta = 1.0;
  define(a_coarseDomain, a_grids, a_refRatios, a_coarsedx, a_bc,
         alpha, aCoef, beta, bCoef);
}
//-----------------------------------------------------------------------

MGLevelOp<LevelData<FArrayBox> >* AVCAMRPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                  int                  a_depth,
                                                                  bool                 a_homoOnly)
{
  CH_TIME("AVCAMRPoissonOpFactory::MGnewOp");

  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref < m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
      {
        break;
      }
    }

  CH_assert(ref !=  m_domains.size()); // didn't find domain

  if (ref > 0)
  {
    dxCrse = m_dx[ref-1];
  }

  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];
  int coarsening = 1;

  for (int i = 0; i < a_depth; i++)
    {
      coarsening *= 2;
      domain.coarsen(2);
    }

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*AVCAMRPoissonOp::s_maxCoarse))
  {
    return NULL;
  }

  dx *= coarsening;

  DisjointBoxLayout layout;
  coarsen_dbl(layout, m_boxes[ref], coarsening);

  Copier ex = m_exchangeCopiers[ref];
  CFRegion cfregion = m_cfregion[ref];

  if (coarsening > 1)
  {
    ex.coarsen(coarsening);
    cfregion.coarsen(coarsening);
  }

  AVCAMRPoissonOp* newOp = new AVCAMRPoissonOp;

  newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  if (a_depth == 0)
    {
      // don't need to coarsen anything for this
      newOp->m_aCoef = m_aCoef[ref];
      newOp->m_bCoef = m_bCoef[ref];
    }
  else
    {
      // need to coarsen coefficients
      RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
      aCoef->define(layout, m_aCoef[ref]->nComp(), m_aCoef[ref]->ghostVect());
      bCoef->define(layout, m_bCoef[ref]->nComp(), m_bCoef[ref]->ghostVect());

      // average coefficients to coarser level
      // for now, do this with a CoarseAverage --
      // may want to switch to harmonic averaging at some point
      CoarseAverage averager(m_aCoef[ref]->getBoxes(),
                             layout, aCoef->nComp(), coarsening);

      CoarseAverageFace faceAverager(m_bCoef[ref]->getBoxes(),
                                     bCoef->nComp(), coarsening);

      if (m_coefficient_average_type == CoarseAverage::arithmetic)
        {
          averager.averageToCoarse(*aCoef, *(m_aCoef[ref]));
          faceAverager.averageToCoarse(*bCoef, *(m_bCoef[ref]));
        }
      else if (m_coefficient_average_type == CoarseAverage::harmonic)
        {
          averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[ref]));
          faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[ref]));
        }
      else
        {
          MayDay::Abort("AVCAMRPoissonOpFactory::MGNewOp -- bad averagetype");
        }

      newOp->m_aCoef = aCoef;
      newOp->m_bCoef = bCoef;
    }

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

AMRLevelOp<LevelData<FArrayBox> >* AVCAMRPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("AVCAMRPoissonOpFactory::AMRnewOp");

  AVCAMRPoissonOp* newOp = new AVCAMRPoissonOp;
  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++)
    {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
      {
        break;
      }
    }

  if (ref == 0)
    {
      // coarsest AMR level
      if (m_domains.size() == 1)
        {
          // no finer level
          newOp->define(m_boxes[0], m_dx[0],
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
        }
      else
        {
          // finer level exists but no coarser
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          newOp->define(m_boxes[0],  m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
        }
    }
  else if (ref ==  m_domains.size()-1)
    {
      dxCrse = m_dx[ref-1];

      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
    }
  else if ( ref == m_domains.size())
    {
      MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

    }
  else
    {
      dxCrse = m_dx[ref-1];

      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1], m_refRatios[ref],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
    }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_aCoef = m_aCoef[ref];
  newOp->m_bCoef = m_bCoef[ref];

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

int AVCAMRPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if (m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }

  if (!found)
    {
      MayDay::Abort("Domain not found in AMR hierarchy");
    }

  return retval;
}

//-----------------------------------------------------------------------
void AVCAMRPoissonOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AVCAMRPoissonOp::
finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                     int a_coarseningFactor)
{
  const AVCAMRPoissonOp& op =
    dynamic_cast<const AVCAMRPoissonOp&>(a_operator);

  // Perform multigrid coarsening on the operator data.
  LevelData<FArrayBox>& acoefCoar = *m_aCoef;
  const LevelData<FArrayBox>& acoefFine = *(op.m_aCoef);
  LevelData<FluxBox>& bcoefCoar = *m_bCoef;
  const LevelData<FluxBox>& bcoefFine = *(op.m_bCoef);
  if (a_coarseningFactor != 1)
  {
    CoarseAverage cellAverage(acoefFine.disjointBoxLayout(),
                              acoefCoar.disjointBoxLayout(),
                              1, a_coarseningFactor);
    for (DataIterator dit = acoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      acoefCoar[dit()].setVal(0.);
    cellAverage.averageToCoarse(acoefCoar, acoefFine);

    CoarseAverageFace faceAverage(bcoefFine.disjointBoxLayout(),
                                  1, a_coarseningFactor);
    for (DataIterator dit = bcoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      bcoefCoar[dit()].setVal(0.);
    faceAverage.averageToCoarse(bcoefCoar, bcoefFine);
  }

  // Handle inter-box ghost cells.
  acoefCoar.exchange();
  bcoefCoar.exchange();

  // Mark the relaxation coefficient dirty.
  m_lambdaNeedsResetting = true;

  // Notify any observers of this change.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
