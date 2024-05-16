#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IceUtility.H"
#include "IceConstants.H"
#include "amrIceF_F.H"
#include "BisiclesF_F.H"
#include "CalvingModel.H"
#include "AMRPoissonOpF_F.H"
#include "AdvectPhysicsF_F.H"
#include "QuadCFInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "CellToEdge.H"
#include "DivergenceF_F.H"
#include "PiecewiseLinearFillPatch.H"
#include "L1L2ConstitutiveRelation.H"
#include "IceThermodynamics.H"
#include "computeSum.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

/// compute RHS for velocity field solve
void IceUtility::defineRHS(Vector<LevelData<FArrayBox>* >& a_rhs,
			   const Vector<RefCountedPtr<LevelSigmaCS > >& a_CS,
			   const Vector<DisjointBoxLayout>& a_grids,
			   const Vector<RealVect>& a_dx)
			    
{

  // options that uses to get set in AmrIce, have established defaults, but only
  // apply to this function
  Real maxRhsDx = -1.0; // don't limit RHS
  bool glCorrection = true; //GL correction (one sided differences) by default
  int nLimited = 0;
  
  {
    // for backward compatibility
    ParmParse pp("amr");
    pp.query("max_rhs_dx", maxRhsDx);
    pp.query("gl_correction", glCorrection);
  }

  {
    ParmParse pp("velocity_rhs");
    pp.query("max_rhs_dx", maxRhsDx);
    pp.query("gl_correction", glCorrection);
  }


  
  CH_assert(a_rhs.size() <= a_CS.size());
  for (int lev = 0; lev < a_rhs.size(); ++lev)
    {
      CH_assert(a_rhs[lev] != NULL && a_CS[lev] != NULL);
      LevelData<FArrayBox>& levelRhs = *a_rhs[lev];
      const LevelSigmaCS& levelCS = *a_CS[lev];
      const DisjointBoxLayout& levelGrids = a_grids[lev];
      Real rhog = levelCS.iceDensity()*levelCS.gravity();
      const RealVect& levelDx = a_dx[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  FArrayBox& rhs = levelRhs[dit];
	  const FArrayBox& gradS = levelCS.getGradSurface()[dit];
	  const FArrayBox& thck = levelCS.getH()[dit];
	  
	  rhs.copy(gradS);
	  for (int dir = 0; dir <SpaceDim ; dir++)
	    {
	      rhs.mult(thck,0,dir);
	    }
	  rhs *= rhog;

	  const bool& anyFloating = levelCS.anyFloating()[dit];
	  
	  if (anyFloating && glCorrection)
	    {
	      const Box& box = levelGrids[dit];
	      const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
	      const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	      
	      for (int dir=0; dir<SpaceDim; dir++)
		{
		  FORT_GLCORRECTION(CHF_FRA1(rhs, dir),
				    CHF_CONST_FRA1(thck,0),
				    CHF_CONST_FRA1(usrf,0),
				    CHF_CONST_FIA1(mask,0),
				    CHF_INT(dir),
				    CHF_CONST_REAL(levelDx[dir]),
				    CHF_CONST_REAL(rhog),
				    CHF_BOX(box));
		}
	    }

	  if (maxRhsDx > 0.0)
	    {

	      
	      for (BoxIterator bit(rhs.box()); bit.ok(); ++bit)
		{
		  Real modRhsDx = 0.0;
		  for (int dir=0; dir<SpaceDim; dir++)
		    {
		      modRhsDx += std::pow( rhs(bit(),dir) * levelDx[dir],2) ;
		    } // end loop over dir
		  modRhsDx = std::sqrt(modRhsDx);
		  if (modRhsDx > maxRhsDx)
		    {
		      //this should be rare, so count
		      nLimited++;
		     
		      Real f = maxRhsDx/modRhsDx;
		      for (int dir=0; dir<SpaceDim; dir++)
			{
			  rhs(bit(),dir)*=f;
			} // end loop over dir
		    } // end if (modRhsDx > maxRhsDx)
		} // end loop over cells
	    } // end if (maxRhsDx > 0.0)
	} // loop over boxes
    }//loop over levels

#ifdef CH_MPI
  {
    int tmp = 1.;
    int result = MPI_Allreduce(&nLimited, &tmp, 1, MPI_INT,
			       MPI_SUM, Chombo_MPI::comm);
    
    CH_assert(result == MPI_SUCCESS);
    if (result != MPI_SUCCESS)
      {
	MayDay::Error("communication error on MPI_Allreduce");
      }
    nLimited = tmp;
  }
#endif
  if (nLimited > 0)
    pout() << "IceUtility::defineRHS: rhs was limited to " << maxRhsDx << "/dx in " << nLimited << " cells " << endl;
  
}



//apply cell-centred helmholtz operator a phi + b grad^2(phi) 
//to cell-centred phi. Assumes that boundary values have been
//set.
void IceUtility::applyHelmOp
(LevelData<FArrayBox>& a_lapPhi,
 const LevelData<FArrayBox>& a_phi, 
 const Real& a_a, const Real& a_b,
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx )
{
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      FORT_OPERATORLAP(CHF_FRA(a_lapPhi[dit]),
		       CHF_FRA(a_phi[dit]),
		       CHF_BOX(a_grids[dit]),
		       CHF_CONST_REAL(a_dx[0]),
		       CHF_CONST_REAL(a_a),
		       CHF_CONST_REAL(a_b));
    }
}



//apply cell-centred  operator div(u)
//to face-centred u. Assumes that ghost cells
//have been set
void IceUtility::applyDiv
(LevelData<FArrayBox>& a_divU,
 const LevelData<FluxBox>& a_u, 
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx
)
{
 
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      a_divU[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  FORT_DIVERGENCE(CHF_CONST_FRA(a_u[dit][dir]),
			  CHF_FRA(a_divU[dit]),
			  CHF_BOX(a_grids[dit]),
			  CHF_CONST_REAL(a_dx[dir]),
			  CHF_INT(dir));
	}
    }
}




//compute face-centered us given cell-centered
//vector u and scalar s 
//Assumes that boundary values have been set.
//USES CENTERED SCHEME
void IceUtility::computeFaceFluxCentered
(LevelData<FluxBox>& a_us,
 const LevelData<FArrayBox>& a_u, 
 const LevelData<FArrayBox>& a_s, 
 const DisjointBoxLayout& a_grids)
{
  CH_assert(a_s.nComp() == 1);
  CH_assert(a_us.nComp() == 1);
  CH_assert(a_u.nComp() == SpaceDim);
  LevelData<FArrayBox> ccus(a_grids,SpaceDim,IntVect::Unit);
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      ccus[dit].copy(a_u[dit]);
      for (int dir = 0; dir < SpaceDim; ++dir) 
	ccus[dit].mult(a_s[dit],0,dir,1);
    }
  CellToEdge(ccus,a_us);
}

//compute face-centered us given cell-centered
//vector u and scalar s 
//Assumes that boundary values have been set.
//USES FIRST ORDER UPWIND SCHEME
void IceUtility::computeFaceFluxUpwind
(LevelData<FluxBox>& a_us,
 const LevelData<FluxBox>& a_u, 
 const LevelData<FArrayBox>& a_s, 
 const DisjointBoxLayout& a_grids)
{
  CH_assert(a_s.nComp() == 1);
  CH_assert(a_us.nComp() == 1);
  CH_assert(a_u.nComp() == 1);
 
  for (DataIterator dit (a_grids); dit.ok(); ++dit)
    {
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  Box faceBox = a_grids[dit].surroundingNodes(dir);
	  FORT_UPWINDFLUXB(CHF_FRA1(a_us[dit][dir],0),
			  CHF_CONST_FRA1(a_u[dit][dir],0),
			  CHF_CONST_FRA1(a_s[dit],0),
			  CHF_CONST_INT(dir),
			  CHF_BOX(faceBox));
	}
    }
}


//apply cell-centred gradient operator grad (phi) . grad(phi)
//to cell-centred phi.Assumes that ghost cells
//have been set
void IceUtility::applyGradSq
(LevelData<FArrayBox>& a_gradPhiSq,
 const LevelData<FArrayBox>& a_phi, 
 const DisjointBoxLayout& a_grids,
 const RealVect& a_dx)
{
	
  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    {
      a_gradPhiSq[dit].setVal(0.0);
      for (int icomp = 0; icomp < a_phi.nComp(); icomp++) 
	{
	  for (int dir =0; dir < SpaceDim; dir++)
	    {
	      Real oneOnTwoDx = 1.0 / (2.0 * a_dx[dir]);
	      for (BoxIterator bit(a_grids[dit]);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  Real g = oneOnTwoDx * (a_phi[dit](iv + BASISV(dir),icomp)
					 - a_phi[dit](iv - BASISV(dir),icomp));
		    a_gradPhiSq[dit](iv) += g*g;
		    
		}
	    }
	}
    }
}

//compute A(x,y,sigma) given internal energy, geometry etc
void IceUtility::computeA
(LevelData<FArrayBox>& a_A,
 const Vector<Real>& a_sigma,
 const LevelSigmaCS& a_coordSys,
 const RateFactor* a_rateFactor,
 const LevelData<FArrayBox>& a_internalEnergy)
{
  const DisjointBoxLayout grids =  a_coordSys.grids();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      // compute A(T)
      for (int layer = 0; layer < a_sigma.size(); ++layer)
	{
	  const Box& box = a_A[dit].box();
	  
	  //1. compute pressure, using p = \rho * g * \sigma * H (should be p = - T_xx - T_yy + \rho * g * \sigma * H)
	  FArrayBox pressure(box,1);
	  pressure.copy(a_coordSys.getH()[dit]);
	  pressure *= a_coordSys.iceDensity() * a_coordSys.gravity() * a_sigma[layer];
	  
	  //2. compute temperature T (and water fraction w) from internal energy and pressure
	  FArrayBox T(box,1);
	  FArrayBox w(box,1);
	  FArrayBox layerE(box,1);
	  layerE.copy(a_internalEnergy[dit],layer,0,1); // an alias would be better?
	  IceThermodynamics::decomposeInternalEnergy(T, w, layerE , pressure, box);
	  //CH_assert(T.min() > 0.0);
	  
	  //3. compute corrected temperature 
	  FArrayBox Tpmp(box,1);
	  Tpmp.copy(pressure); Tpmp *= -IceThermodynamics::icepmeltfactor(); // Tpmp += triplepoint - TINY_NORM;
	  T -= Tpmp;

	  //4. compute A
	  FArrayBox layerA;
	  layerA.define(Interval(layer,layer),a_A[dit]);
	  a_rateFactor->computeA(layerA,T,pressure,box);
	}
    }
}


void IceUtility::computeC0(Vector<LevelData<FArrayBox>* >& a_vectC0,
			   const Vector<LevelData<FArrayBox>* >& a_vectC,
			   const Vector<DisjointBoxLayout>& a_grids,
			   const Vector<RefCountedPtr<LevelSigmaCS> >& a_coordSys,
			   const Vector<Real> a_dx, int a_finest_level)
 
{

  bool wallDrag = true; //compute additional drag due to contact with rocky walls
  Real wallDragExtra = 0.0; // assume wall drag proportional to basal drag only;
  Real thinIceDragExtra = 0.0; // mimimum *linear* drag for thin ice
  Real thinIceDragThickness = 0.0; // how thin is thin?
  bool thinIceDrag = false;
  {
    //backward compatiblity
    ParmParse ppAmr("amr");
    ppAmr.query("wallDrag",wallDrag);
    ppAmr.query("wallDragExtra",wallDragExtra);
  }
  
  {
    ParmParse pp("wall_drag");
    pp.query("basic", wallDrag);
    pp.query("extra", wallDragExtra);
  }
  
  {
    ParmParse pp("thin_ice_drag");
    pp.query("extra", thinIceDragExtra);
    pp.query("thickness", thinIceDragThickness);
    thinIceDrag = thinIceDragExtra * thinIceDragThickness > TINY_NORM;
  }
  
  for (int lev=0; lev <= a_finest_level; lev++)
    {
      const LevelSigmaCS& levelCS = *a_coordSys[lev];
      const DisjointBoxLayout& grids = a_grids[lev];
      for (DataIterator dit(grids); dit.ok(); ++dit)
        {
	  FArrayBox& thisC0 = (*a_vectC0[lev])[dit];
	  const FArrayBox& thisC = (*a_vectC[lev])[dit];
	  thisC0.setVal(0.0);
	  if (wallDrag)
	    {
	      IceUtility::addWallDrag(thisC0, levelCS.getFloatingMask()[dit], 
				  levelCS.getSurfaceHeight()[dit], levelCS.getH()[dit], 
				  levelCS.getTopography()[dit], thisC, wallDragExtra,
				  RealVect::Unit*a_dx[lev], grids[dit]);
	    }

	  if (thinIceDrag)
	    {
	      IceUtility::addThinIceDrag(thisC0, levelCS.getFloatingMask()[dit],
					 levelCS.getH()[dit], thinIceDragExtra,
					 thinIceDragThickness, grids[dit]);
	    }
	}
    }
}
			   


void IceUtility::addThinIceDrag(FArrayBox& a_drag, 
				const BaseFab<int>& a_mask,
				const FArrayBox& a_thk,
				const Real& a_extra,
				const Real& a_thin,
				const Box& a_box)
{
  
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (a_thk(iv) < a_thin)
	a_drag(iv) += a_extra;
    }
}


void IceUtility::addWallDrag(FArrayBox& a_drag, 
			      const BaseFab<int>& a_mask,
			      const FArrayBox& a_usrf,
			      const FArrayBox& a_thk,
			      const FArrayBox& a_topg,
			      const FArrayBox& a_beta,
			      const Real& a_extra, 
			      const RealVect& a_dx,
			      const Box& a_box)
{
  

  
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      for (int sign = -1; sign <= 1; sign+=2)
	{

	  for (BoxIterator bit(a_box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if (a_mask(iv) == GROUNDEDMASKVAL || a_mask(iv) == FLOATINGMASKVAL)
		{
		  const IntVect ivp = iv + sign*BASISV(dir);
		  if (a_mask(ivp) != GROUNDEDMASKVAL && a_mask(ivp) != FLOATINGMASKVAL)
		    {
		      Real contact = 
			std::min(a_thk(iv) * 0.5,  0.5*(a_topg(ivp)+a_topg(iv)) 
				 - (a_usrf(iv)- a_thk(iv)));
		      if (contact > 0.0)
			{
			  a_drag(iv) += (a_beta(iv) + a_extra) * contact / a_dx[dir];
			}
		    }
		}
	    }
	}
    }
}

/// extrapolate face centered velocity field (usually derived by cell-to-face average) 
/// to the margins
void IceUtility::extrapVelocityToMargin(LevelData<FluxBox>& a_faceVel, 
					const LevelData<FArrayBox>& a_cellVel, 
					const LevelSigmaCS& a_coordSys)
{

  //modification to fluxes at the margins, that is where mask changes to open sea or land.
  //
  // -----|-----|-----|-----|-----|
  //      |     |     |     |     |
  //   x  o  x  o  x  F     |     |
  //      |     |     |     |     |
  // -----|-----|-----|-----|-----|
  //
  //  n-2   n-2    n    n+1
  //
  // 2D case (L1L2 / SSA)
  // there are valid values of the basal velocity at points x and
  // and the z-varying velocity at points o. The points at o 
  // have been interpolated from the x, and should vary little vertically
  // We need the flux at F, but since there is no x_n+1, the interpolated
  // value makes no sense. Using the margin boundary condition to get
  // du/dx at the face is perhaps the ideal approach , but for now we just
  // take x_{n-1} and o_{n-1} and extrapolate.
  // 
  // On top of that, the face velocity is reduced by a factor 
  // f = min( (surface(n)-topography(n+1) , thickness(n) ), / thickness(n)
  // which prevents ice from flowing up vertical walls

  const DisjointBoxLayout& grids = a_cellVel.disjointBoxLayout();

  // this calculation needs at least one ghost cell in the face-centered
  // and cell centered velocity FABs
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      CH_assert(a_cellVel.ghostVect()[dir] >= 1);
      CH_assert(a_faceVel.ghostVect()[dir] >= 1);
    }

 
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      const FArrayBox& cellVel = a_cellVel[dit];
      //const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
      const FArrayBox& usrf = a_coordSys.getSurfaceHeight()[dit];
      const FArrayBox& topg = a_coordSys.getTopography()[dit];
      const FArrayBox& thk = a_coordSys.getH()[dit];

      for (int dir = 0; dir <SpaceDim; dir++)
	{
	  Box faceBox = grids[dit];
	  faceBox.surroundingNodes(dir);
	  FArrayBox& faceVel = a_faceVel[dit][dir];
	  
	  //{Real maxFaceVelocity = faceVel.norm(faceBox,0);
	  //CH_assert(maxFaceVelocity < 0.5 * HUGE_VEL);}
	  
	  Box grownFaceBox = faceBox;
	  CH_assert(faceVel.box().contains(grownFaceBox));
	  FArrayBox vface(faceBox,1);
	  FArrayBox faceVelCopy(faceVel.box(), 1); faceVelCopy.copy(faceVel);

	  FORT_EXTRAPTOMARGIN(CHF_FRA1(faceVel,0), CHF_FRA1(vface,0),
			      CHF_CONST_FRA1(faceVelCopy,0),
			      CHF_CONST_FRA1(cellVel,dir),
			      CHF_CONST_FRA1(usrf,0),
			      CHF_CONST_FRA1(topg,0),
			      CHF_CONST_FRA1(thk,0),
			      CHF_CONST_INT(dir),
			      CHF_BOX(faceBox));
	  //{Real maxFaceVelocity = faceVel.norm(faceBox,0);
	  //  CH_assert(maxFaceVelocity < HUGE_VEL);}
	}

    }
}
    

void IceUtility::computeFaceVelocity
(LevelData<FluxBox>& a_faceVelAdvection,
 LevelData<FluxBox>& a_faceVelTotal,
 LevelData<FluxBox>& a_faceDiffusivity,
 LevelData<FArrayBox>& a_cellDiffusivity,
#if BISICLES_Z == BISICLES_LAYERED
 LevelData<FluxBox>& a_layerXYFaceXYVel,
 LevelData<FArrayBox>& a_layerSFaceXYVel,
#endif
 const LevelData<FArrayBox>& a_velocity,
 const LevelSigmaCS& a_coordSys,
 const IceThicknessIBC* a_iceThicknessIBC,
 const LevelData<FArrayBox>& a_A,
#if BISICLES_Z == BISICLES_LAYERED
 const LevelData<FArrayBox>& a_sA,
 const LevelData<FArrayBox>& a_bA,
#endif			 
 const LevelData<FArrayBox>* a_crseVelocity,
 const LevelData<FArrayBox>* a_crseDiffusivity,
 int a_nRefCrse,
 const ConstitutiveRelation* a_constitutiveRelation,
 bool a_additionalVelocity,  bool a_implicitDiffusion)
{			   
  
  //We need to copy the cell-centered velocity onto
  //LevelData<FArrayBox> with 2 ghost cells, because we will 
  //need one ghost cell's worth of face-centred data.
  IntVect grownVelGhost(2*IntVect::Unit);
  const DisjointBoxLayout& grids = a_velocity.disjointBoxLayout();
  LevelData<FArrayBox> grownVel(grids,a_velocity.nComp(), grownVelGhost);
  
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      grownVel[dit].setVal(0.0);
      grownVel[dit].copy(a_velocity[dit], a_velocity[dit].box());
      a_cellDiffusivity[dit].setVal(0.0);
      for (int dir = 0; dir < SpaceDim; dir++)
	a_faceDiffusivity[dit][dir].setVal(0.0);
    }
  
  if (a_crseVelocity != NULL)
    {
      const DisjointBoxLayout& crseGrids = a_crseVelocity->disjointBoxLayout();
      PiecewiseLinearFillPatch velFiller(grids , crseGrids, a_velocity.nComp(), 
					 crseGrids.physDomain(), a_nRefCrse,
					 grownVelGhost[0]);
	
      Real time_interp_coeff = 0.0;
      velFiller.fillInterp(grownVel, *a_crseVelocity, *a_crseVelocity,
			   time_interp_coeff,0, 0, a_velocity.nComp());
	
      //slc. In L1L2 we will compute second derivatives of v, hence
      //we need quadratic coarse-fine interpolation here.
      Real dx = 1.0;
      QuadCFInterp qcfi(grids , &crseGrids, dx, a_nRefCrse, 
			a_velocity.nComp(), grids.physDomain());
      qcfi.coarseFineInterp(grownVel, *a_crseVelocity);

    }

  grownVel.exchange();
      
  //default calculation : average to faces 
  CellToEdge(grownVel, a_faceVelAdvection);

  //correct margins, where the face average does not make sense
  IceUtility::extrapVelocityToMargin(a_faceVelAdvection, grownVel, a_coordSys);

#if BISICLES_Z == BISICLES_LAYERED
  //for layered models (SSA,L1L2) assume du/dz = 0
  for (int j = 0; j < a_layerXYFaceXYVel.nComp(); ++j)
    a_faceVelAdvection.copyTo(Interval(0,0), a_layerXYFaceXYVel, Interval(j,j));
  
  for (int j = 0; j < a_layerSFaceXYVel.nComp(); j+=SpaceDim)
    grownVel.copyTo(Interval(0,SpaceDim-1), a_layerSFaceXYVel,Interval(j,j+SpaceDim-1)); 
#endif

  //allow the thickness/velocity bc to modify the face velocities 
  if (a_iceThicknessIBC != NULL)
    {
      a_iceThicknessIBC->modifyFaceVelocity(a_faceVelAdvection, a_coordSys, grids.physDomain() );
    }

  // copy faceVelAdvection into faceVelTotal - no diffusion term at this stage
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  a_faceVelTotal[dit][dir].copy(a_faceVelAdvection[dit][dir]);
	}
    }
				       
#if BISICLES_Z == BISICLES_LAYERED
  //copy the basal velocity into the vertically varying velocities.
  {
    
    for (DataIterator dit(grids); dit.ok(); ++dit)
      {

	

	FArrayBox& SFaceXYVel = a_layerSFaceXYVel[dit];
	const FArrayBox& cellVel = grownVel[dit];
	for (int ic = 0; ic < SFaceXYVel.nComp()-1; ic+=SpaceDim)
	  {
	    SFaceXYVel.copy( cellVel , 0, ic, SpaceDim);
	  }

	for (int dir = 0; dir < SpaceDim; ++dir)
	  {
	    const FArrayBox& faceVel = a_faceVelAdvection[dit][dir];
	    FArrayBox& XYFaceXYVel = a_layerXYFaceXYVel[dit][dir];
	    
	    for (int ic = 0; ic < XYFaceXYVel.nComp(); ic++)
	      {
		XYFaceXYVel.copy( faceVel , 0, ic, 1);
	      }
	  }
      }
  }
#endif 

  

  if (a_additionalVelocity)
    {
      const L1L2ConstitutiveRelation* L1L2Ptr = dynamic_cast<const L1L2ConstitutiveRelation*>(a_constitutiveRelation);
      if (L1L2Ptr != NULL)
	{
	  //L1L2Ptr->computeFaceFluxVelocity(grownVel, a_crseVelocity, a_nRefCrse, a_coordSys, 
	  //				   grids,  grids.physDomain(), a_A, a_sA, a_bA,
	  //				   a_faceVelTotal ,a_layerXYFaceXYVel, a_layerSFaceXYVel);
	 

	  
	  L1L2Ptr->modifyTransportCoefficients(grownVel, a_crseVelocity, a_crseDiffusivity, a_nRefCrse, a_coordSys, 
					       grids,  grids.physDomain(), a_A, a_sA, a_bA,
					       a_faceVelAdvection, a_faceVelTotal, a_faceDiffusivity,
					       a_cellDiffusivity, a_layerXYFaceXYVel, a_layerSFaceXYVel);
	  if (!a_implicitDiffusion)
	    {
	      // in this case, faceVelTotal contains the velocity
	      for (DataIterator dit(grids);dit.ok();++dit)
		{
		  for (int dir = 0; dir < SpaceDim; dir++)
		    { 
		      a_faceVelAdvection[dit][dir].copy(a_faceVelTotal[dit][dir] );
		    }
		}
	    }
	}
    }
  
  


  a_faceVelAdvection.exchange();
  a_faceVelTotal.exchange();
#if BISICLES_Z == BISICLES_LAYERED
  a_layerXYFaceXYVel.exchange();
  a_layerSFaceXYVel.exchange();
#endif

 

}

/// compute the cross layer velocity u^sigma (the contravariant component)
void IceUtility::computeSigmaVelocity
(LevelData<FArrayBox>& a_uSigma,
 const LevelData<FluxBox>& a_layerThicknessFlux,
 const LevelData<FArrayBox>& a_layerSFaceXYVel,
 const LevelData<FArrayBox>& a_dHdt,
 const DisjointBoxLayout a_grid,
 const LevelData<FArrayBox>& a_surfaceThicknessSource,
 const LevelData<FArrayBox>& a_basalThicknessSource,
 const Vector<Real>& a_dSigma,
 const RealVect& a_dx,
 const Real& a_dt )
{
  
  int nLayer = a_dSigma.size();
  CH_assert(nLayer == a_layerThicknessFlux.nComp());
  CH_assert((SpaceDim*(nLayer+1) == a_layerSFaceXYVel.nComp()));
  for (DataIterator dit(a_grid); dit.ok(); ++dit)
    {
      const Box& box = a_grid[dit];
      
      // this copy perhaps indicates layer should run faster than
      // dir in sFaceXYVel, but for now ...
      const FArrayBox& sFaceXYVel = a_layerSFaceXYVel[dit];
      FArrayBox uX(box, nLayer+1);
      FArrayBox uY(box, nLayer+1);
      
      for (int l = 0; l < nLayer + 1; l++)
	{
	  uX.copy(sFaceXYVel, l*SpaceDim, l);
	  uY.copy(sFaceXYVel, l*SpaceDim + 1, l);
	}
      
      FArrayBox divUHxy(box, nLayer);
      divUHxy.setVal(0.0);
      
      //const RealVect& dx = a_coordSysNew.dx(); 
      for (int dir =0; dir < SpaceDim; dir++)
	{
	  const FArrayBox& uH = a_layerThicknessFlux[dit][dir];
	  FORT_DIVERGENCE(CHF_CONST_FRA(uH),
			  CHF_FRA(divUHxy),
			  CHF_BOX(box),
			  CHF_CONST_REAL(a_dx[dir]),
			  CHF_INT(dir));
	}
      

      
      //surface and basal thickness source
      const FArrayBox& bts = a_basalThicknessSource[dit];
      const FArrayBox& sts = a_surfaceThicknessSource[dit];

      // sigma-componnet of velocity at layer faces
      FArrayBox& uSigma = a_uSigma[dit]; 
      uSigma.setVal(0.0);
      FORT_COMPUTESIGMAVEL(CHF_FRA(uSigma),
			   CHF_CONST_FRA(uX),
			   CHF_CONST_FRA(uY),
			   CHF_CONST_FRA(divUHxy),
			   CHF_CONST_VR(a_dSigma),
			   CHF_CONST_FRA1(a_dHdt[dit],0),
			   CHF_CONST_FRA1(sts,0),
			   CHF_CONST_FRA1(bts,0),
			   CHF_CONST_INT(nLayer),
			   CHF_BOX(box));
      
      
      CH_assert(uSigma.norm(0) < HUGE_NORM);
    } //end compute vertical velocity loop over boxes
      
  a_uSigma.exchange();
      
}




///Identify regions of fast ice  and eliminate them.
/**
   implies a call to IceUtility::eliminateRemoteIce iff any fast ice is removed
 */
int IceUtility::eliminateFastIce
(Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 Vector<LevelData<FArrayBox>* >& a_vel,
 Vector<LevelData<FArrayBox>* >& a_calvedIce,
 Vector<LevelData<FArrayBox>* >& a_addedIce,
 Vector<LevelData<FArrayBox>* >& a_removedIce,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<ProblemDomain>& a_domain,
 const Vector<int>& a_refRatio, Real a_crseDx,
 int a_finestLevel, int a_maxIter, Real a_thinIceTol, Real a_fastIceTol, 
 bool a_edgeOnly, int a_verbosity)
{
  CH_TIME("IceUtility::eliminateFastIce");
  if (a_verbosity > 3)
    {
      pout() << "IceUtility::eliminateFastIce" << endl;
    }
  int nEliminated = 0;
  Real fastIceTolSq = a_fastIceTol * a_fastIceTol;

  for (int iter = 0; iter < 10; iter++)
    {

  for (int lev=0; lev <= a_finestLevel ; ++lev)
    {
      for (DataIterator dit(a_grids[lev]); dit.ok(); ++dit)
	{
	  FArrayBox& H = a_coordSys[lev]->getH()[dit];
	  const BaseFab<int>& mask = a_coordSys[lev]->getFloatingMask()[dit];
	  FArrayBox& calved = (*a_calvedIce[lev])[dit];
	  FArrayBox& added = (*a_addedIce[lev])[dit];
	  FArrayBox& removed = (*a_removedIce[lev])[dit];
	  FArrayBox& u = (*a_vel[lev])[dit];
	  FArrayBox HH(H.box(),1); HH.copy(H);
	  
	  for (BoxIterator bit(a_grids[lev][dit]);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      Real usq = D_TERM(u(iv,0)*u(iv,0), + u(iv,1)*u(iv,1), u(iv,2)*u(iv,2));
	      if (usq > fastIceTolSq)
		{
		  bool elim = !a_edgeOnly;
		  int dir = 0;
		  while ( (!elim) && dir < SpaceDim) {
		    elim |= (HH(iv+BASISV(dir)) < a_thinIceTol);
		    elim |= (HH(iv-BASISV(dir)) < a_thinIceTol);
		    dir++;
		  }
		    
		  if (elim)
		    {
		      Real prevThck = H(iv);
		      H(iv) = 0.0;
		      D_DECL(u(iv,0) = 0 ,u(iv,1) = 0, u(iv,2) = 0);
		      if (a_verbosity > 5)
			{
			  pout() << " (fast) eliminated level " << lev << " iv " << iv << std::endl;
			}
		      nEliminated++;
		      // Record gain/loss of ice
		      CalvingModel::updateCalvedIce(H(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		      //removed(iv) += (prevThck-H(iv));

		    }
		}
	    }
	}
    }
  
  if (a_verbosity > 4)
	{
	  pout() << "... (this cpu) eliminated " << nEliminated << " cells " << endl;
	}
    }
#ifdef CH_MPI
  {
    int tmp = 1.;
    int result = MPI_Allreduce(&nEliminated, &tmp, 1, MPI_INT,
			       MPI_SUM, Chombo_MPI::comm);
    
    CH_assert(result == MPI_SUCCESS);
    if (result != MPI_SUCCESS)
      {
	MayDay::Error("communication error on MPI_Allreduce");
      }
    nEliminated = tmp;
  }
#endif  

  if (a_verbosity > 3)
	{
	  pout() << "... (all cpus) eliminated " << nEliminated << " cells " << endl;
	}

  if (nEliminated > 0)
    {
      
      // eliminateRemoteIce will recompute surface elevation etc
      eliminateRemoteIce(a_coordSys,a_vel,
			 a_calvedIce,a_addedIce,a_removedIce,
			 a_grids,a_domain,a_refRatio, a_crseDx,
			 a_finestLevel,a_maxIter,a_thinIceTol,  a_verbosity);
    }
  return nEliminated;
}

///Identify regions of floating ice that are remote
///from grounded ice and eliminate them.
/**
   Regions of floating ice unconnected to land 
   lead to an ill-posed problem, with a zero
   basal traction coefficient and Neumann boundary
   conditions. Here, we attempt to identify them
   by carrying out a procedure that in the worst case can be 
   O(N^2). 
*/ 
void IceUtility::eliminateRemoteIce
(Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 Vector<LevelData<FArrayBox>* >& a_vel,
 Vector<LevelData<FArrayBox>* >& a_calvedIce,
 Vector<LevelData<FArrayBox>* >& a_addedIce,
 Vector<LevelData<FArrayBox>* >& a_removedIce,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<ProblemDomain>& a_domain,
 const Vector<int>& a_refRatio, Real a_crseDx,
 int a_finestLevel, int a_maxIter, Real a_tol, int a_verbosity)
{
  CH_TIME("IceUtility::eliminateRemoteIce");
  if (a_verbosity > 3)
    {
      pout() << "IceUtility::eliminateRemoteIce" << endl;
    }
  //Define phi = 1 on grounded ice, 0 elsewhere
  Vector<LevelData<FArrayBox>* > phi(a_finestLevel + 1, NULL);
  
  for (int lev=0; lev <= a_finestLevel ; ++lev)
    {
      const DisjointBoxLayout& levelGrids = a_grids[lev];
      LevelSigmaCS& levelCS = *a_coordSys[lev];
      phi[lev] = new LevelData<FArrayBox>(levelGrids,1,IntVect::Unit);
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  //const FArrayBox& thisH = levelCS.getH()[dit];
	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  FArrayBox& thisPhi = (*phi[lev])[dit];
	  thisPhi.setVal(0.0);
	  Real a = 1.0;
	  int b = GROUNDEDMASKVAL;
	  FORT_SETONMASK(CHF_FRA1(thisPhi,0),
			 CHF_CONST_FIA1(mask,0),
			 CHF_CONST_INT(b),CHF_CONST_REAL(a),
			 CHF_BOX(levelGrids[dit]));
			 
	}
      phi[lev]->exchange();
    }
 
  int iter = 0; 
  Real sumPhi = computeSum(phi, a_refRatio ,a_crseDx, Interval(0,0), 0);
  Real oldSumPhi = 0.0;
  
  do {
    oldSumPhi = sumPhi;
    
    
    for (int lev=0; lev <= a_finestLevel ; ++lev)
      {
	LevelData<FArrayBox>& levelPhi = *phi[lev];
	LevelSigmaCS& levelCS = *a_coordSys[lev];
	const DisjointBoxLayout& levelGrids = a_grids[lev];
	
	if (lev > 0)
	  {
	    //fill ghost cells
	    PiecewiseLinearFillPatch ghostFiller(levelGrids, 
						 a_grids[lev-1],
						 1, 
						 a_domain[lev-1],
						 a_refRatio[lev-1],
						 1);
	    
	    ghostFiller.fillInterp(levelPhi, *phi[lev-1],*phi[lev-1] , 0.0, 0, 0, 1);
	    
	  }
	
	for  (DataIterator dit(levelGrids); dit.ok(); ++dit)
	  {
	    //sweep in all four directions, copying phi = 1 into cells with thickness > tol 
	    //Real tol = 1.0;
	    FORT_SWEEPCONNECTED2D(CHF_FRA1(levelPhi[dit],0),
				  CHF_CONST_FRA1(levelCS.getH()[dit],0),
				  CHF_CONST_REAL(a_tol), 
				  CHF_BOX(levelGrids[dit]));
	  }

	
	levelPhi.exchange();
      }
    
    
    sumPhi = computeSum(phi, a_refRatio, a_crseDx, Interval(0,0), 0);

    if (a_verbosity > 3)
      {
	pout() << "IceUtility::eliminateRemoteIce iteration " << iter 
	       << " connected cells = " << sumPhi << endl;
      }
    
    iter++;
  } while ( iter < a_maxIter && sumPhi > oldSumPhi );

 
  //now destroy the doomed regions, and reset a_coordSys
  for (int lev=0; lev <= a_finestLevel ; ++lev)
    {
      const LevelData<FArrayBox>& levelPhi = *phi[lev];
      LevelSigmaCS& levelCS = *a_coordSys[lev];
      const DisjointBoxLayout& levelGrids = a_grids[lev];
      
      for (DataIterator dit(levelGrids); dit.ok(); ++dit)
	{
	  const FArrayBox& thisPhi = levelPhi[dit];
	  FArrayBox& h = levelCS.getH()[dit];
	  FArrayBox& calved = (*a_calvedIce[lev])[dit];
	  FArrayBox& added = (*a_addedIce[lev])[dit];
	  FArrayBox& removed = (*a_removedIce[lev])[dit];
	  FArrayBox& u = (*a_vel[lev])[dit];

	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
	  for (BoxIterator bit(levelGrids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real prevThck = h(iv);
	      if (mask(iv) == FLOATINGMASKVAL && thisPhi(iv) < 0.5)
	      	{
		  h(iv) = 0.0;
		  D_DECL(u(iv,0) = 0 ,u(iv,1) = 0, u(iv,2) = 0);
		  if (a_verbosity > 5)
		  pout() << " (remote) eliminated level " << lev << " iv " << iv << std::endl;
	      	}
	      // Record gain/loss of ice
	      CalvingModel::updateCalvedIce(h(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	      //removed(iv) += (prevThck-h(iv));
	      /*	      if (h(iv) > prevThck)
		{
		  added(iv) += (prevThck-h(iv));
		}
	      else 
		{
		  if (mask(iv) == OPENLANDMASKVAL || mask(iv) == FLOATINGMASKVAL)
		    {
		      removed(iv) += (prevThck-h(iv));
		    }
		  else
		    {
		      removed(iv) += (prevThck-h(iv));
		    }
		} 
	      */
	    }
	}
      levelCS.getH().exchange();
      
      levelCS.getH().exchange();
      LevelSigmaCS* crseCS = (lev > 0)?&(*a_coordSys[lev-1]):NULL;
      int refRatio = (lev > 0)?a_refRatio[lev-1]:-1;
      levelCS.recomputeGeometry(crseCS, refRatio);     
    }

  for (int lev=0; lev <= a_finestLevel ; ++lev)
    {
      if (phi[lev] != NULL)
 	{
 	  delete phi[lev];
 	  phi[lev] = NULL;
 	}
    }

}

/// multiply a_u by the grounded portion of each cell
/**
   For a_subdivision > 0, the grounded portion of a cell
   is computed by (1) approximating the thickness above/below
   flotation in each quarter of each cell with a bilinear formula 
   h(x,y) (2) integrating the Heaviside function (h > 0)?1:0 over 
   each quarter  using the midpoint rule in [2^a_n]^2 subdivisions
   (making 4*[2^a_n]^2 subcells)

   a_subdivision = 2 (subcell size = dx/8) is probably plenty, 
   the cost is obviously O(a_subdivision^2)
 */
void IceUtility::multiplyByGroundedFraction
(LevelData<FArrayBox>& a_u, 
 const LevelSigmaCS& a_coords,
 const DisjointBoxLayout& a_grids,
 int a_subdivision)
{
  CH_TIME("IceUtility::multiplyByGroundedFraction");
  CH_assert(a_subdivision > 0);

  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    { 
      Box grownBox = a_grids[dit];
      grownBox.grow(1);
      FArrayBox hab(grownBox,1);
      // thickness above / under flotation (< 0 in shelf)
      Real rhoi = a_coords.iceDensity();
      Real rhoo = a_coords.waterDensity();
      Real seaLev = a_coords.seaLevel();
      Real habmin, habmax;
      FORT_HOVERUNDERFLOTATION
	( CHF_FRA1(hab,0),CHF_REAL(habmin),CHF_REAL(habmax),
	  CHF_CONST_FRA1(a_coords.getH()[dit],0),
	  CHF_CONST_FRA1(a_coords.getTopography()[dit],0),
	  CHF_CONST_REAL(rhoi),
	  CHF_CONST_REAL(rhoo),
	  CHF_CONST_REAL(seaLev),
	  CHF_BOX(grownBox));
      
      if ( (habmin < 0.0) && (habmax > 0.0))
	{
	  //int nRef = std::pow(2,a_subdivision); // std:pow(int,int) is not a real thing
	  int nRef = 1;
	  for (int i = 0; i < a_subdivision ;i++)
	    {
	      nRef *= 2;
	    }
	  //integrate (hab>0)?1:0 over each cell to approximate grounded area
	  FArrayBox ag(a_grids[dit],1);
	  FORT_INTEGRATEHEAVISIDE2D( CHF_FRA1(ag,0), CHF_CONST_FRA1(hab,0),
				     CHF_CONST_INT(nRef), CHF_BOX(a_grids[dit]));
	  
	  //weight u
	  a_u[dit] *= ag;
	}
    }
}


/// set C = 0 in floating region
/** 
    If nRef <= 0, C is set to zero according to the mask in a_levelCS

    Otherwise Hab, the thickness above/below flotation (> 0 on grounded ice, < 0 in the shelf)
    is computed, and then interpolated (using the bilinear formula : this is not
    the usual conservative Chombo formula. C is then multplied by  the integral of (Hab > 0)?1:0 
    over the each cell.  
 
 */ 
void IceUtility::setFloatingBasalFriction
(LevelData<FArrayBox>& a_C, const LevelSigmaCS& a_coords,
 const DisjointBoxLayout& a_grids)
{
  CH_TIME("IceUtility::setFloatingBasalFriction");

  // options that uses to get set in AmrIce, have established defaults, but only
  // apply to this function. 
  int subdivision = 0; // > 0 for sub-grid interpolation
  
  {
    // for backward compatibility
    ParmParse pp("amr");
    pp.query("grounding_line_subdivision", subdivision);
  }

  {
    ParmParse pp("basal_friction");
    pp.query("grounding_line_subdivision", subdivision);
  }
 
  if (subdivision > 0)
    {
      pout() << " IceUtility::setFloatingBasalFriction : interpolation... " << std::endl;
      IceUtility::multiplyByGroundedFraction(a_C,  a_coords, a_grids,  subdivision);
      pout() << " IceUtility::setFloatingBasalFriction : done interpolation " << std::endl;
    }
  
  //finally, set C = 0 in any floating cells - covers nRef == 0 and nRef > 0
  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    {
      bool anyFloating = a_coords.anyFloating()[dit]; 
      // set friction on open land and open sea to 100. Set friction on floating ice to 0
      if(anyFloating)
	{
	  FArrayBox& thisC = a_C[dit];
	  const BaseFab<int>& mask = a_coords.getFloatingMask()[dit];
	  FORT_SETFLOATINGBETA(CHF_FRA1(thisC,0),
			       CHF_CONST_FIA1(mask,0),
			       CHF_BOX(a_grids[dit]));
	  
	  // friction must be non-negative
	  CH_assert(thisC.min(a_grids[dit]) >= 0.0); 
	}  
    }
}


// ///Identify regions connected to the open ocean
// /**
//    Identify cells which are connected to open ocean via cells
//    with lsrf - topg > a_tol, where open ocean is any region where
//    a_topg < a_oceanTopg and ice is either floating or absent.
//    by carrying out a procedure that in the worst case can be 
//    O(N^2). 

//    fills the Vector<LevelData<FArrayBox>* > a_oceanDepth with
//    water depth (lower surface - bedrock) for connected cells
//    and zero elsewhere
// */ 
// void IceUtility::computeConnectedOceanDepth
// (Vector<LevelData<FArrayBox>* >& a_oceanDepth,
//  const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
//  const Vector<DisjointBoxLayout>& a_grids,
//  const Vector<ProblemDomain>& a_domain,
//  const Vector<int>& a_refRatio, Real a_crseDx,
//  int a_finestLevel, int a_maxIter, Real a_tol, Real a_topg int a_verbosity)
// {
//   CH_TIME("IceUtility::computeConnectedOceanDepth");
//   if (a_verbosity > 0)
//     {
//       pout() << "IceUtility::computeConnectedOceanDepth" << endl;
//     }
//   //fill depth=1 in areas that are open ocean by defintion topg(iv) < a_oceanTopg & not grounded
//   for (int lev=0; lev <= a_finestLevel ; ++lev)
//     {
//       const DisjointBoxLayout& levelGrids = a_grids[lev];
//       LevelSigmaCS& levelCS = *a_coordSys[lev];
//       for (DataIterator dit(levelGrids); dit.ok(); ++dit)
// 	{
// 	  const BaseFab<int>& mask = levelCS.getFloatingMask()[dit];
// 	  const FArrayBox& usrf = levelCS.getSurfaceHeight()[dit];
// 	  const FArrayBox& topg = levelCS.getTopography()[dit];
// 	  const FArrayBox& thck = levelCS.getH()[dit];
// 	  FArrayBox& depth =  (*a_oceanDepth[lev])[dit];
// 	  depth.setVal(0.0);
// 	  for (BoxIterator bit(depth.box()); bit.ok(); ++bit)
// 	    {
// 	      if (mask(iv) != GROUNDEDMASKVAL)
// 		{
// 		  if (topg(iv) < a_oceanTopg)
// 		    {
// 		      depth(iv) = 1.0;
// 		    }
// 		}
// 	    }
// 	}
//     }
  
//   //flood fill the connected regions
//   int iter = 0; 
//   Real sumOcean = computeSum(a_oceanDepth, a_refRatio ,a_crseDx, Interval(0,0), 0);
//   Real oldSumOcean = 0.0;
  
//   do {
//     oldSumOcean = sumOcean;
    
//     for (int lev=0; lev <= a_finestLevel ; ++lev)
//       {
// 	LevelData<FArrayBox>& levelOcean = *a_oceanDepth[lev];
// 	LevelSigmaCS& levelCS = *a_coordSys[lev];
// 	const DisjointBoxLayout& levelGrids = a_grids[lev];
	
// 	if (lev > 0)
// 	  {
// 	    //fill ghost cells
// 	    PiecewiseLinearFillPatch ghostFiller(levelGrids, 
// 						 a_grids[lev-1],
// 						 1, 
// 						 a_domain[lev-1],
// 						 a_refRatio[lev-1],
// 						 1);
	    
// 	    ghostFiller.fillInterp(levelOcean, *a_oceanDepth[lev-1],*a_oceanDepth[lev-1] , 0.0, 0, 0, 1);
	    
// 	  }
	
// 	for  (DataIterator dit(levelGrids); dit.ok(); ++dit)
// 	  {
// 	    //sweep in all four directions, copying phi = 1 into cells with water depth > tol 
// 	    FArrayBox wd(levelGrids[dit],1);
// 	    wd.copy(levelCS.getSurfaceHeight()[dit]);
// 	    wd.minus(levelCS.getH()[dit]);
// 	    wd.minus(levelCS.getTopography()[dit]);
	    
// 	    FORT_SWEEPCONNECTED2D(CHF_FRA1(levelOcean[dit],0),
// 				  CHF_CONST_FRA1(wd,0),
// 				  CHF_CONST_REAL(a_tol), 
// 				  CHF_BOX(levelGrids[dit]));
// 	  }

	
// 	levelOcean.exchange();
//       }
    
    
//     sumOceanDepth = computeSum(a_oceanDepth, a_refRatio, a_crseDx, Interval(0,0), 0);
    
//     if (a_verbosity > 0)
//       {
// 	pout() << "IceUtility::computeConnectedOceanDepth iteration " << iter 
// 	       << " ocean volume s = " << sumOcenDepth << endl;
//       }
    
//     iter++;
//   } while ( iter < a_maxIter && sumOceanDpeth > oldSumOcenDepth );

// }

#include "NamespaceFooter.H"
