#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "L1L2ConstitutiveRelation.H"
#include "LevelMappedDerivatives.H"
#include "ConstitutiveRelationF_F.H"
#include "L1L2ConstitutiveRelationF_F.H"
#include "BisiclesF_F.H"
#include "IceConstants.H"
#include "FineInterp.H"
#include "BoxIterator.H"
#include "ExtrapBCF_F.H"
#include "ExtrapGhostCells.H"
#include "ReflectGhostCells.H"
#include "EdgeToCell.H"
#include "CellToEdge.H"
#include "MayDay.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"
#if BISICLES_Z == BISICLES_LAYERED
#define L1L2_DEFAULT_NLAYER 8

ConstitutiveRelation* 
L1L2ConstitutiveRelation::getNewConstitutiveRelation() const
{

  L1L2ConstitutiveRelation* newPtr = new L1L2ConstitutiveRelation;
  newPtr->m_solverTol = m_solverTol;
  newPtr->m_additionalVelocitySIAGradSLimit = m_additionalVelocitySIAGradSLimit;
  newPtr->m_effectiveViscositySIAGradSLimit = m_effectiveViscositySIAGradSLimit;
  newPtr->m_additionalVelocitySIAOnly = m_additionalVelocitySIAOnly;
  newPtr->m_startFromAnalyticMu = m_startFromAnalyticMu;
  newPtr->m_layerCoarsening = m_layerCoarsening;
  GlensFlowRelation* gfr =  newPtr->getGlensFlowRelationPtr();
  gfr->setParameters(glensFlowRelation.m_n, 
		     glensFlowRelation.m_epsSqr0, 
		     glensFlowRelation.m_delta);
  return static_cast<ConstitutiveRelation*>(newPtr);

}


void
L1L2ConstitutiveRelation::parseParameters()
{
  ParmParse ppL1L2("L1L2");
  ppL1L2.query("solverTolerance", m_solverTol);
  ppL1L2.query("effectiveViscositySIAGradSLimit", m_effectiveViscositySIAGradSLimit);
  if (m_effectiveViscositySIAGradSLimit > 0.0)
    {
      m_additionalVelocitySIAGradSLimit = std::min(m_effectiveViscositySIAGradSLimit,m_additionalVelocitySIAGradSLimit);
    }
  ppL1L2.query("additionalVelocitySIAGradSLimit", m_additionalVelocitySIAGradSLimit);
  ppL1L2.query("additionalVelocitySIAOnly", m_additionalVelocitySIAOnly);
  ppL1L2.query("startFromAnalyticMu",m_startFromAnalyticMu);
  ppL1L2.query("layerCoarsening",m_layerCoarsening);
  
  //underlying glen's flow flaw parameters
  Real n = 3.0;
  ppL1L2.query("n",n);
  Real epsSqr0 = 1.0e-12;
  ppL1L2.query("epsSqr0",epsSqr0);
  Real delta = 0.0;
  ppL1L2.query("delta",delta);
  glensFlowRelation.setParameters(n, epsSqr0, delta);


  
  
}

void
L1L2ConstitutiveRelation::computeMu(LevelData<FArrayBox>& a_mu,
                                    const LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                                    const LevelData<FArrayBox>* a_crseVelPtr,
                                    int a_nRefCrse,
                                    const LevelData<FArrayBox>& a_A,
                                    const LevelSigmaCS& a_coordSys,
				    const ProblemDomain& a_domain,
                                    const IntVect& a_ghostVect) const
{

  CH_TIME("L1L2ConstitutiveRelation::computeMu");
  int nLayer = a_A.nComp();
  // either mu in each layer or average mu makes sense
  CH_assert(a_A.nComp() == a_mu.nComp() || a_mu.nComp() == 1); 
  CH_assert(a_A.nComp() == a_coordSys.getSigma().size());
  //mu at layer boundaries
  const DisjointBoxLayout& grids = a_mu.getBoxes();
  
  const Vector<Real>& sigma = a_coordSys.getSigma();
  
  LevelData<FArrayBox> gradVel(grids,SpaceDim*SpaceDim, a_ghostVect);
  LevelData<FArrayBox> epsSqr(grids,1, a_ghostVect);

  
  if (a_A.nComp() == a_mu.nComp() || a_A.nComp() == 1)
    {
      computeMuZ(a_mu, gradVel, epsSqr, sigma, a_vel, a_scale,  a_crseVelPtr, a_nRefCrse,
		 a_A, a_coordSys, a_domain, a_ghostVect);
    }
  else 
    { 
      LevelData<FArrayBox> muz(grids, nLayer, a_ghostVect);
      computeMuZ(muz, gradVel, epsSqr, sigma, a_vel, a_scale,  a_crseVelPtr, a_nRefCrse,
		 a_A, a_coordSys, a_domain, a_ghostVect); 
      
      // midpoint rule integration over sigma to find integral(mu,dz)/H
      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  FArrayBox& thisMu = a_mu[dit];
	  thisMu.setVal(0.0);
	  for (int l = 0; l < nLayer ; l++)
	    {
	      thisMu.plus(muz[dit] , a_coordSys.getDSigma()[l], l, 0);
	    }
	}// end loop over grids
    }
} 

void
L1L2ConstitutiveRelation::computeDissipation(LevelData<FArrayBox>& a_dissipation,
					     const LevelData<FArrayBox>& a_vel, 
					     const LevelData<FArrayBox>* a_crseVelPtr,
					     int a_nRefCrse,
					     const LevelData<FArrayBox>& a_A,
					     const LevelSigmaCS& a_coordSys,
					     const ProblemDomain& a_domain,
					     const IntVect& a_ghostVect) const
{
  // L1L2 dissiaption : since \sigma^{L1L2} _ij(x,y,z) = 2 * mu^{L1L2}(x,y,z) * \epsilon ^{base}_(x,y),
  // we can compute \Phi(x,y,z) = 4 \mu^{L1L2}(x,y,z} (\epsilon ^{base}_(x,y))^2 much as in the 
  // SSA case.
  CH_assert(a_dissipation.nComp() == a_A.nComp());
  for (DataIterator dit = a_dissipation.dataIterator(); dit.ok(); ++dit)
    {
      a_dissipation[dit].setVal(0.0);
    }
  computeMu(a_dissipation , a_vel, 1.0,  a_crseVelPtr, a_nRefCrse, a_A, a_coordSys, 
	    a_domain, a_ghostVect);

  LevelData<FArrayBox> epsSqr;
  epsSqr.define(a_vel.disjointBoxLayout(),1,IntVect::Zero);
  computeStrainRateInvariant(epsSqr, a_vel, a_crseVelPtr,
			     a_nRefCrse,a_coordSys, a_ghostVect);

  
  for (DataIterator dit = a_dissipation.dataIterator(); dit.ok(); ++dit)
    {
      for (int layer=0; layer < a_dissipation.nComp(); ++layer)
	{
	  a_dissipation[dit].mult(epsSqr[dit],0,layer,1);
	}
      a_dissipation[dit].mult(4.0);
    }

}

void  
L1L2ConstitutiveRelation::computeFaceMu(LevelData<FluxBox>& a_mu,
                                        LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                                        const LevelData<FArrayBox>* a_crseVelPtr,
                                        int a_nRefCrse,
                                        const LevelData<FluxBox>& a_A,
                                        const LevelSigmaCS& a_coordSys,
					const ProblemDomain& a_domain,
                                        const IntVect& a_ghostVect) const
{
  CH_TIME("L1L2ConstitutiveRelation::computeFaceMu");
  int nLayer = a_coordSys.getSigma().size();
  CH_assert(a_A.nComp() == nLayer);

  const DisjointBoxLayout& grids = a_mu.getBoxes();

  LevelData<FluxBox> mu(grids,nLayer, a_ghostVect);
  const Vector<Real>& sigma = a_coordSys.getSigma();
 
  LevelData<FluxBox> gradVel(grids, SpaceDim*SpaceDim, a_ghostVect);
  LevelData<FluxBox> epsSqr(grids, 1, a_ghostVect);

  computeFaceMuZ(mu, gradVel, epsSqr, sigma, a_vel, a_scale, a_crseVelPtr, a_nRefCrse,
                 a_A, a_coordSys, a_domain, a_ghostVect);


  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisMu = a_mu[dit];
      thisMu.setVal(0.0);
      
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
	  //midpoint rule inegration over sigma to find integral(mu,dz)/H
	  for (int l = 0; l < nLayer; l++)
	    {
	      Real s =  a_coordSys.getDSigma()[l];
	      thisMu[faceDir].plus(mu[dit][faceDir], s ,l ,0 );
	    }
          
        } // end loop over face directions
    } // end loosisp over grids

}


void 
L1L2ConstitutiveRelation::computeMuZ(LevelData<FArrayBox>& a_mu,
                                     LevelData<FArrayBox>& a_gradVel,
                                     LevelData<FArrayBox>& a_epsSqr,
                                     const Vector<Real>& a_sigma,
                                     const LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                                     const LevelData<FArrayBox>* a_crseVelPtr,
                                     int a_nRefCrse,
                                     const LevelData<FArrayBox>& a_A,
                                     const LevelSigmaCS& a_coordSys,
				     const ProblemDomain& a_domain,
                                     const IntVect& a_ghostVect) const 
{
  CH_TIME("L1L2ConstitutiveRelation::computeMuZ");

  // first, compute grad(u) and epsilon^2 (at the base)
  computeStrainRateInvariant(a_epsSqr, a_gradVel, a_vel,
                             a_crseVelPtr, a_nRefCrse, a_coordSys,
                             a_ghostVect); 
  const LevelData<FArrayBox>& H =  a_coordSys.getH();
  const LevelData<FArrayBox>& G = a_coordSys.getGradSurface();
  const DisjointBoxLayout& grids = a_mu.getBoxes();
  Box domBox = a_domain.domainBox();
  DataIterator dit = a_mu.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box muBox = grids[dit];
      muBox.grow(a_ghostVect);
      muBox &= domBox;
      
      computeEitherMuZ(a_mu[dit], a_sigma,  G[dit], a_epsSqr[dit], a_A[dit], a_scale,  
		       H[dit] ,a_coordSys.iceDensity()*a_coordSys.gravity(), muBox,m_layerCoarsening);
    }

  a_mu.exchange();

}


void
L1L2ConstitutiveRelation::computeFaceMuZ(LevelData<FluxBox>& a_mu,
                                         LevelData<FluxBox>& a_gradVel,
                                         LevelData<FluxBox>& a_epsSqr,
                                         const Vector<Real>& a_sigma,
                                         LevelData<FArrayBox>& a_vel,  const Real& a_scale,
                                         const LevelData<FArrayBox>* a_crseVelPtr,
                                         int a_nRefCrse,
                                         const LevelData<FluxBox>& a_A,
                                         const LevelSigmaCS& a_coordSys,
					 const ProblemDomain& a_domain,
                                         const IntVect& a_ghostVect) const
{

  CH_TIME("L1L2ConstitutiveRelation::computeFaceMuZ");
  // first, compute grad(u) and epsilon^2 at cell faces (at the base)
  computeStrainRateInvariantFace(a_epsSqr,  a_gradVel,  a_vel, a_crseVelPtr,
                                 a_nRefCrse, a_coordSys, a_ghostVect);
				 
  const LevelData<FluxBox>& gradSface = a_coordSys.getGradSurfaceFace();
  const LevelData<FluxBox>& faceH = a_coordSys.getFaceH();

  const DisjointBoxLayout& grids = a_mu.getBoxes();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box grownBox = grids[dit];
      grownBox.grow(a_ghostVect);
      for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
        {
          Box faceBox = grownBox;
          faceBox.surroundingNodes(faceDir);

          computeEitherMuZ(a_mu[dit][faceDir], a_sigma,  
                           gradSface[dit][faceDir],
                           a_epsSqr[dit][faceDir], a_A[dit][faceDir],  a_scale, 
                           faceH[dit][faceDir], 
			   a_coordSys.iceDensity()*a_coordSys.gravity(),
			   faceBox,m_layerCoarsening);  
        }

    }



}





void  
L1L2ConstitutiveRelation::computeEitherMuZ(FArrayBox& a_mu,
                                           const Vector<Real>& a_sigma,
                                           const FArrayBox& a_grads,
                                           const FArrayBox& a_epsSqr,
                                           const FArrayBox& a_A, const Real& a_scale,
                                           const FArrayBox& a_H,
					   const Real& a_rhog,
                                           const Box& a_box,
					   int a_layerCoarsen) const
{

  CH_TIME("L1L2ConstitutiveRelation::computeEitherMuZ");

  CH_assert(a_sigma.size() == a_mu.nComp());
  CH_assert(a_A.nComp() == 1 || a_sigma.size() == a_A.nComp());

  if (a_layerCoarsen > 0 && a_sigma.size()%2 == 0)
    {
      // Save time by computing mu in only half of the layers
      int nCrseLayer = a_sigma.size() / 2;
      FArrayBox crseMu(a_mu.box(), nCrseLayer);
      FArrayBox crseA(a_mu.box(), nCrseLayer);
      Vector<Real> crseSigma(nCrseLayer);

      for (int crseLayer = 0, layer = 0; crseLayer < nCrseLayer; crseLayer++, layer+=2)
	{
	  crseA.copy(a_A, layer, crseLayer, 1);
	  crseA.plus(a_A, layer+1, crseLayer, 1);
	  crseSigma[crseLayer] = 0.5*( a_sigma[layer] + a_sigma[layer+1]) ;
	}
      crseA *= 0.5;
      computeEitherMuZ(crseMu, crseSigma, a_grads, a_epsSqr,  crseA, a_scale, a_H, a_rhog, a_box, a_layerCoarsen-1); 

      for (int crseLayer = 0, layer = 0; crseLayer < nCrseLayer; crseLayer++, layer+=2)
	{
	  a_mu.copy(crseMu, crseLayer, layer, 1);
	  a_mu.copy(crseMu, crseLayer, layer+1, 1);
	}

      nCrseLayer = 0;
      
    }
  else
    {
      //SIA approximation to vertical contribution to second invariant of
      //deviatoric stress (rho * g * H * grad(s))^2 
      const Real& rhog = a_rhog;
      FArrayBox phiTildeSqr(a_box,1);
      int nComp = a_grads.nComp();
  
      FArrayBox grads(a_box,nComp);
      grads.copy(a_grads,0,0,nComp);
      if (m_effectiveViscositySIAGradSLimit > 0.0)
	{
	  FORT_L1L2MODLIMIT(CHF_FRA(grads),
			    CHF_CONST_REAL(m_effectiveViscositySIAGradSLimit),
			    CHF_BOX(a_box),
			    CHF_CONST_INT(nComp));
	}
      
      FORT_L1L2PHITILDESQR(CHF_FRA1(phiTildeSqr,0),
			   CHF_CONST_FRA1(a_H,0),
			   CHF_CONST_FRA(grads),
			   CHF_CONST_REAL(rhog),
			   CHF_INT(nComp),
			   CHF_BOX(a_box));
      
      
      
      FArrayBox epsSqr(a_box,1);
      epsSqr.copy(a_epsSqr);
      
      FArrayBox A(a_box,1);
      FArrayBox res(a_box,1);
      
      const int maxIter = 150;;
      
      for (unsigned int l = 0; l < a_mu.nComp(); ++l){
	
	
	A.copy(a_A,l,0);
	A *= a_scale; // scale A 
	if (l == 0)
	  {
	    // could be the top of the ice sheet (if sigma[l] == 0.0)
	    // midpoint of layer 0 or 1
	    
	    //calling FORT_COMPUTEL1L2MU with sigma = 0.0 will 
	    //get us Glen's mu 
	    Real sigma0 = 0.0;
	    
	    FORT_COMPUTEL1L2MU(CHF_FRA1(a_mu,l),
			       CHF_FRA1(res,0),
			       CHF_CONST_FRA1(A, 0),
			       CHF_CONST_FRA1(epsSqr,0),
			       CHF_CONST_FRA1(phiTildeSqr,0),
			       CHF_BOX(a_box),
			       CHF_CONST_REAL(sigma0),
			       CHF_CONST_REAL(glensFlowRelation.m_n),
			       CHF_CONST_REAL(glensFlowRelation.m_epsSqr0),
			       CHF_CONST_REAL(glensFlowRelation.m_delta),
			       CHF_CONST_REAL(m_solverTol),
			       CHF_CONST_INT(maxIter));
	  }
	
	if (a_sigma[l] > TINY_NORM)
	  {
	    // various cases below the surface of the ice, which 
	    // differ  only in the two values used to initialize the 
	    // secant solve in FORT_COMPUTEL1L2MU
	    //Real factor = 0.95;
	    if (l == 0)
	      {
		//use Glen's mu, twice - the secant solve will check and change this
		res.copy(a_mu,l,0);
	      } 
	    else if ( l == 1)
	      {
		// use mu from the previous layer, twice
		res.copy(a_mu,l-1,0); 
		a_mu.copy(a_mu,l-1,l);
		//	    a_mu.mult(0.95,l,1);
		if (m_startFromAnalyticMu)
		  FORT_ANALYTICL1L2MU(CHF_FRA1(a_mu,l),
				      CHF_CONST_FRA1(A, 0),
				      CHF_CONST_FRA1(epsSqr,0),
				      CHF_CONST_FRA1(phiTildeSqr,0),
				      CHF_BOX(a_box),
				      CHF_CONST_REAL(a_sigma[l]),
				      CHF_CONST_REAL(glensFlowRelation.m_n),
				      CHF_CONST_REAL(glensFlowRelation.m_epsSqr0));
		
	      } 
	    else 
	      {
		// use mu from the previous two layers
		res.copy(a_mu,l-2,0);
		a_mu.copy(a_mu,l-1,l);
		
		if (m_startFromAnalyticMu)
		  FORT_ANALYTICL1L2MU(CHF_FRA1(a_mu,l),
				      CHF_CONST_FRA1(A, 0),
				      CHF_CONST_FRA1(epsSqr,0),
				      CHF_CONST_FRA1(phiTildeSqr,0),
				      CHF_BOX(a_box),
				      CHF_CONST_REAL(a_sigma[l]),
				      CHF_CONST_REAL(glensFlowRelation.m_n),
				      CHF_CONST_REAL(glensFlowRelation.m_epsSqr0));
		
	      } 
	    
	    
	    FORT_COMPUTEL1L2MU(CHF_FRA1(a_mu,l),
			       CHF_FRA1(res,0),
			       CHF_CONST_FRA1(A, 0),
			       CHF_CONST_FRA1(epsSqr,0),
			       CHF_CONST_FRA1(phiTildeSqr,0),
			       CHF_BOX(a_box),
			       CHF_CONST_REAL(a_sigma[l]),
			       CHF_CONST_REAL(glensFlowRelation.m_n),
			       CHF_CONST_REAL(glensFlowRelation.m_epsSqr0),
			       CHF_CONST_REAL(glensFlowRelation.m_delta),
			       CHF_CONST_REAL(m_solverTol),
			       CHF_CONST_INT(maxIter));
	    
	    
	    
	    Real resNorm = res.norm(0);
	    CH_assert(resNorm < 2.0 * m_solverTol);
	    CH_assert(a_mu.min(a_box,l) > 0.0e0);
	  } // end general case
      } // end loop over layers
    } // end else if a_layerCoarsen = 0
}

// compute the velocities at layer faces. To do so, we need A
// at the layer interfaces (for now, rethink this) 
void 
L1L2ConstitutiveRelation::computeFaceFluxVelocity(LevelData<FArrayBox>& a_vel,
                                                  const LevelData<FArrayBox>* a_crseVelPtr,
                                                  int a_nRefCrse,
						  const LevelData<FArrayBox>& a_A,
						  const LevelData<FArrayBox>& a_thickness,
						  const RealVect& a_dx,
						  LevelData<FluxBox>& a_fluxVel,
						  LevelData<FluxBox>& a_layerXYFaceXYVel,
						  LevelData<FArrayBox>& a_layerSFaceXYVel,
                                                  const LevelSigmaCS& a_coordSys,
						  const ProblemDomain& a_domain,
						  const IntVect& a_cellGhost) const
{
   CH_TIME("L1L2ConstitutiveRelation::computeFaceFluxVelocity");
   const DisjointBoxLayout& grids = a_vel.getBoxes();  
   DataIterator dit = grids.dataIterator();   
   const Vector<Real>& sigma = a_coordSys.getFaceSigma();
   int nLayer = sigma.size() - 1;

   CH_assert(a_A.nComp() == nLayer + 1); 
   // need the flow law coefficient (Glenn's A) at (i + 1/2 , j, l + 1/2) and (i , j + 1/2, l + 1/2)
   IntVect extraGhost = a_cellGhost + IntVect::Unit;
  
   // effective viscosity mu(sigma,x[,y]) at cell centres,
   // plus grad(u)(1,x,y) and epsilon^2(1,x,y)
   LevelData<FArrayBox> mu(grids, sigma.size(), extraGhost);
   LevelData<FArrayBox> gradVel(grids, SpaceDim*SpaceDim, extraGhost);
   LevelData<FArrayBox> epsSqr(grids, 1, extraGhost);
   for (dit.begin(); dit.ok(); ++dit)
     {
       mu[dit].setVal(0.0);       
       gradVel[dit].setVal(0.0);
     }


   computeMuZ(mu, gradVel, epsSqr, sigma, a_vel, 1.0,  a_crseVelPtr, a_nRefCrse, a_A, 
              a_coordSys,  a_domain, extraGhost);
   
   // We need to do something at the domain boundaries
   // and reflection is more common than other choices.
   for (int dir = 0; dir < SpaceDim; ++dir)
     {
       ReflectGhostCells(mu, a_domain, dir, Side::Lo );
       ReflectGhostCells(mu, a_domain, dir, Side::Hi );
       ReflectGhostCells(gradVel, a_domain, dir, Side::Lo );
       ReflectGhostCells(gradVel, a_domain, dir, Side::Hi );
       ReflectGhostCells(epsSqr, a_domain, dir, Side::Lo );
       ReflectGhostCells(epsSqr, a_domain, dir, Side::Hi );
     }
  
   const LevelData<FluxBox>& faceG = a_coordSys.getGradSurfaceFace();
   const LevelData<FluxBox>& faceH = a_coordSys.getFaceH();

   // grad(u) at cell faces
   LevelData<FluxBox> faceGradU(grids, SpaceDim*SpaceDim, a_cellGhost);
   //epsilion^2 at cell faces
   LevelData<FluxBox> faceEpsSqr(grids, 1,a_cellGhost);
   // mu at cell faces
   LevelData<FluxBox> faceMu(grids,sigma.size(), a_cellGhost);
   // A at cell faces
   LevelData<FluxBox> faceA(grids, a_A.nComp(), extraGhost);
   CellToEdge(a_A, faceA );
   
   
   computeFaceMuZ(faceMu,faceGradU,faceEpsSqr,sigma,
                  a_vel, 1.0, a_crseVelPtr,a_nRefCrse,
                  faceA,a_coordSys,a_domain,a_cellGhost);
       
   for (dit.begin(); dit.ok(); ++dit)
     {
       Box cellBox = grids[dit];
       cellBox.grow(a_cellGhost);
       
       //const FArrayBox& thisH = a_thickness[dit]; 
       const FluxBox& thisFaceH = faceH[dit];
       const FluxBox& thisFaceG = faceG[dit];

       Box outerCellBox = cellBox;
       outerCellBox.grow(1);
       
       CH_assert(SpaceDim == 2);//ought to work for 3(2) and 2(1) at some point
      
       for (int faceDir=0; faceDir<SpaceDim; faceDir++)
         {
           Box faceBox = cellBox;
           faceBox.surroundingNodes(faceDir);
	   // vertical part of the stress tensor T, phi_zi 
	   // at the i = a_faceDir cell face
	   // dir runs faster than layer

	   FArrayBox phi(faceBox, SpaceDim * sigma.size());
           phi.setVal(0.0);
	   
	   if (!m_additionalVelocitySIAOnly)
	     {
	       // z-derivative of (non-gravitational) phi_zi at the i cell face
	       FArrayBox dphi(faceBox, SpaceDim*sigma.size());
	       dphi.setVal(0.0);
	       const RealVect& dx = a_dx;
	       
	       const int nn = derivComponent(faceDir,faceDir);
	       int transDir = (faceDir + 1)%SpaceDim;
	       const int tt = derivComponent(transDir,transDir);
	       const int nt = derivComponent(faceDir,transDir);
	       const int tn = derivComponent(transDir,faceDir);
	       
	       const FArrayBox& thisGradVel = gradVel[dit];
	       const FArrayBox& thisMu = mu[dit];
	       
	       for (int l = 0; l < sigma.size(); ++l)
		 {
		   FArrayBox ldphi; // dphi at this sigma
		   ldphi.define(Interval(SpaceDim*l,SpaceDim*l+SpaceDim-1),dphi);
		   
		   //this should probably go into a FORTRAN kernel.
		   for (BoxIterator bit(faceBox); bit.ok(); ++bit)
		     {
		       const IntVect& i = bit(); 
		       IntVect iw = i - BASISV(faceDir);
		       if (true){
			 
			 ldphi(i,faceDir) = 4.0/dx[faceDir] * 
			   ( thisMu(i,l) * thisGradVel(i,nn)  - 
			     thisMu(iw,l) * thisGradVel(iw,nn));
		       }
                   
		       if (true && SpaceDim == 2)
			 {
			   
			   ldphi(i,faceDir) += 2.0/dx[faceDir] * 
			     (thisMu(i,l) * thisGradVel(i,tt)  - 
			      thisMu(iw,l) * thisGradVel(iw,tt));
                       
			   
			   IntVect inw = iw + BASISV(transDir);
			   IntVect isw = iw - BASISV(transDir);
			   IntVect in = i + BASISV(transDir);
			   IntVect is = i - BASISV(transDir);
			   
			   ldphi(i,faceDir) += 0.25/dx[faceDir] * 
			     (thisMu(in,l) * (thisGradVel(in,nt)+thisGradVel(in,tn))
			      + thisMu(inw,l) * (thisGradVel(inw,nt)+thisGradVel(inw,tn))
			      -thisMu(is,l) * (thisGradVel(is,nt)+thisGradVel(is,tn))
                          - thisMu(isw,l) * (thisGradVel(isw,nt)+thisGradVel(isw,tn)));
			   
                       
			   ldphi(i,transDir) =  1.0/dx[faceDir] * 
			     ( thisMu(i,l) * (thisGradVel(i,nt) + thisGradVel(i,tn))   - 
			       thisMu(iw,l) * (thisGradVel(iw,nt) + thisGradVel(iw,tn)));
			   
			   ldphi(i,transDir) += 1.0/dx[faceDir] * 
			     (thisMu(in,l) * thisGradVel(in,nn)
			      + thisMu(inw,l) * thisGradVel(inw,nn)
			      -thisMu(is,l) * thisGradVel(is,nn) 
			      -  thisMu(isw,l) * thisGradVel(isw,nn));
                       
			   ldphi(i,transDir) += 0.5/dx[faceDir] * 
			     (thisMu(in,l) * thisGradVel(in,tt)
			      + thisMu(inw,l) * thisGradVel(inw,tt)
			      -thisMu(is,l) * thisGradVel(is,tt) 
			      - thisMu(isw,l) * thisGradVel(isw,tt));
                       
			 }
		       
		     }
		   const FArrayBox& thisFaceHDir = thisFaceH[faceDir];
		   ldphi.mult(thisFaceHDir, 0, faceDir);
		   ldphi.mult(thisFaceHDir, 0, transDir);
		   
		 }
	    
	       // use the trapezium rule to integrate dphi
	       for (int l = 1; l < sigma.size(); ++l){
		 Real s = 0.5 * (sigma[l] - sigma[l-1]);
		 phi.plus(dphi, s , SpaceDim*(l - 1), SpaceDim * l , SpaceDim);
		 phi.plus(dphi, s , SpaceDim*(l), SpaceDim * l , SpaceDim);
		 phi.plus(phi, 1.0, SpaceDim*(l - 1), SpaceDim * l , SpaceDim);
	       }
           
	     } // end if !m_additionalVelocitySIAOnly

           // now add rho*g*ds/dx[i] to phi_zi
           FArrayBox faceGdir(faceBox, SpaceDim);
	   faceGdir.copy(thisFaceG[faceDir],0,0,SpaceDim);
	   if (m_additionalVelocitySIAGradSLimit > 0.0)
	     {
	       FORT_L1L2MODLIMIT(CHF_FRA(faceGdir),
				 CHF_CONST_REAL(m_additionalVelocitySIAGradSLimit),
				 CHF_BOX(faceBox),
				 CHF_CONST_INT(SpaceDim));
	     }
           const FArrayBox& faceHdir = thisFaceH[faceDir];

           for (int l = 0; l < sigma.size(); ++l){
             Real rgSigma =  a_coordSys.iceDensity() * a_coordSys.gravity() * sigma[l];
             FArrayBox lphi;
             lphi.define(Interval(SpaceDim*l,SpaceDim*l+(SpaceDim-1)),phi);
             
             FORT_L1L2ADDGSTRESS(CHF_FRA(lphi),
                                 CHF_CONST_FRA(faceGdir),
                                 CHF_CONST_FRA1(faceHdir,0),
                                 CHF_CONST_REAL(rgSigma),
                                 CHF_BOX(faceBox),
                                 CHF_CONST_INT(SpaceDim));

           }
           
           //can now fill dphi with du[dir]/dz
           FArrayBox du(faceBox, SpaceDim*sigma.size());
           du.setVal(0.0);
           for (int l = nLayer; l >= 0; --l){
             FArrayBox ldu; 
             ldu.define(Interval(SpaceDim*l,SpaceDim*l+(SpaceDim-1)),du);
             FArrayBox lphi; 

             lphi.define(Interval(SpaceDim*l,SpaceDim*l+(SpaceDim-1)),phi);
             
             FORT_L1L2UIGRAND(CHF_FRA(ldu),
                              CHF_CONST_FRA1(faceHdir,0),
                              CHF_CONST_FRA1(faceA[dit][faceDir],l),
                              CHF_CONST_FRA1(faceMu[dit][faceDir],l),
                              CHF_CONST_FRA1(faceEpsSqr[dit][faceDir],0),
                              CHF_CONST_FRA(lphi),
                              CHF_CONST_REAL(glensFlowRelation.m_n),
                              CHF_CONST_INT(SpaceDim),
                              CHF_BOX(faceMu[dit][faceDir].box()));
             
           } // end loop over layers
           

         
           // use the trapezium rule to integrate du and get u - u_base
           FArrayBox vel(du.box(),du.nComp());
	   //FArrayBox& vel = a_layerFaceVel[dit][faceDir];
           vel.setVal(0.0);
           for (int l = nLayer - 1; l >= 0; --l){
             Real s = 0.5 * (sigma[l+1] - sigma[l]);
             vel.plus(du, s , SpaceDim*(l + 1), SpaceDim * l , SpaceDim);
             vel.plus(du, s , SpaceDim*(l), SpaceDim * l , SpaceDim);  
             vel.plus(vel, 1.0 , SpaceDim*(l + 1), SpaceDim * l , SpaceDim);	     
           }

	   Real vel_max = vel.norm(vel.box(), 0, 0,  nLayer);
	   if (vel_max > 2.0e+2)
	     {
	       pout() << "L1L2ConstititiveRelation::computeFaceFlux  max(|u - u_base|) = "
		 << vel_max << std::endl; 
	     }
	   
	   //a_layerSFaceXYVel contains the basal velocity at cell centers
	   //on entry, so average u - u_base from cell faces and add
	   for (int l = 0; l < nLayer + 1 ; l++)
	     {
	       int comp = SpaceDim*l + faceDir;
	       FArrayBox ccVel;
	       ccVel.define(Interval(comp,comp),a_layerSFaceXYVel[dit]);
	       FArrayBox fcExtra;
	       fcExtra.define(Interval(comp,comp),vel);
	       Box ccBox = fcExtra.box();
	       ccBox.growHi(faceDir,-1);
	       for (BoxIterator bit(ccBox); bit.ok(); ++bit)
		 {
		   const IntVect& iv = bit();
		   ccVel(iv) += 0.5* (fcExtra(iv)+fcExtra(iv+BASISV(faceDir)));
		 }

	     }

           //FArrayBox& a_fluxVel contains the basal velocity
           //at cell faces on entry.
           //CH_assert(a_fluxVel.norm(0,0,a_fluxVel.nComp()) < HUGE_NORM)
	   for (int l = 0; l < nLayer; l++)
	     {
	       a_layerXYFaceXYVel[dit][faceDir].copy(a_fluxVel[dit][faceDir],0,l);
	       a_layerXYFaceXYVel[dit][faceDir].plus(vel,0.5,l*SpaceDim + faceDir,l,1);
	       a_layerXYFaceXYVel[dit][faceDir].plus(vel,0.5,(l+1)*SpaceDim + faceDir,l,1);
	       CH_assert(a_layerXYFaceXYVel[dit][faceDir].norm(faceBox,0,l,1) < HUGE_NORM);
	     }
	   
           for (int l = 0; l < nLayer; l++)
             {
               Real s = 0.5 * (sigma[l+1] - sigma[l]);
               a_fluxVel[dit][faceDir].plus(vel,faceBox,faceBox,s,
                                        l*SpaceDim + faceDir, 0, 1);
               a_fluxVel[dit][faceDir].plus(vel,faceBox,faceBox,s,
                                        (l+1)*SpaceDim + faceDir, 0, 1);
             }
	   


           CH_assert(a_fluxVel[dit][faceDir].norm(faceBox,0,0,a_fluxVel.nComp()) < HUGE_NORM);

         } // end loop over face directions

     } // end loop over grids
}



void
L1L2ConstitutiveRelation::modifyTransportCoefficients
(const LevelData<FArrayBox>& a_cellVel,
 const LevelData<FArrayBox>* a_crseVelPtr,
 const LevelData<FArrayBox>* a_crseDiffusivityPtr,
 int a_nRefCrse,
 const LevelSigmaCS& a_coordSys,
 const DisjointBoxLayout& a_grids,
 const ProblemDomain& a_domain,
 const LevelData<FArrayBox>& a_A,
 const LevelData<FArrayBox>& a_sA,
 const LevelData<FArrayBox>& a_bA,
 LevelData<FluxBox>& a_faceVelAdvection,
 LevelData<FluxBox>& a_faceVelTotal,
 LevelData<FluxBox>& a_faceDiffusivity,
 LevelData<FArrayBox>& a_cellDiffusivity,
 LevelData<FluxBox>& a_layerXYFaceXYVel,
 LevelData<FArrayBox>& a_layerSFaceXYVel) const
{
  CH_TIME("L1L2ConstitutiveRelation::modifyTransportCoefficients");

  //assumes that the usual interpolation of the base velocity to cell faces has been done
  //so that a_faceVelAdvection and a_faceVelTotal contain a_cellVel averaged to faces
  
  IntVect grownGhost = IntVect::Zero;
  grownGhost += 2*IntVect::Unit;
  LevelData<FArrayBox> sigmaFaceA(a_grids,a_coordSys.getFaceSigma().size(),grownGhost);
  LevelData<FArrayBox> grownThickness(a_grids, 1, grownGhost);
  const LevelData<FArrayBox>& H = a_coordSys.getH();
  
  // (DFM) if cellVel's ghosting is greater than grownGhost, we can 
  // skip this copy, but putting it here for now for simplicity..
  LevelData<FArrayBox> grownCellVel(a_grids, SpaceDim, grownGhost);
  // do fab-by-fab copies to preserve ghosting
  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      grownThickness[dit].copy(H[dit]);
      grownCellVel[dit].copy(a_cellVel[dit]);
      //a_cellVel.copyTo(grownCellVel);
    }
  
  //followed by L1L2 corrections
  for (DataIterator dit(a_grids) ; dit.ok(); ++dit)
    {
      int nLayer = a_coordSys.getSigma().size();
      //ice surface
      sigmaFaceA[dit].copy(a_sA[dit],0,0);
      for (int l = 1; l < nLayer; ++l)
	{
	  //interior : average mid-layer flow law coefficient (Glenn's A)s to layer interfaces
	  sigmaFaceA[dit].copy(a_A[dit],l-1,l);
	  sigmaFaceA[dit].plus(a_A[dit],l,l);
	  sigmaFaceA[dit].mult(0.5,l);
	}
      //ice base
      sigmaFaceA[dit].copy(a_bA[dit],0,nLayer);
    }
  
 
  RealVect dx = a_coordSys.dx();
  IntVect ghostVect = IntVect::Zero;

  LevelData<FluxBox> vertAverageVelFace(a_grids, 1 , IntVect::Zero);
  for (DataIterator dit(a_grids);dit.ok();++dit)
    {
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  vertAverageVelFace[dit][dir].copy(a_faceVelAdvection[dit][dir]);
	}
    }

  //find u(z) and \bar{u}
  computeFaceFluxVelocity(grownCellVel, a_crseVelPtr, a_nRefCrse,
			  sigmaFaceA, grownThickness, dx,  
			  vertAverageVelFace, a_layerXYFaceXYVel,  a_layerSFaceXYVel,
			  a_coordSys, a_domain, ghostVect);

  //find average u(z) - u(b) =  \bar{u'} - u(z=b)
  for (DataIterator dit(a_grids);dit.ok();++dit)
    {
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  vertAverageVelFace[dit][dir].minus(a_faceVelAdvection[dit][dir]);
	}
    }

  //leave u_advection = u_base, set D = \bar{u'}H / grad(H),  u_total = u_base + \bar{u'}
  //LevelData<FArrayBox> cellDiffusivity(a_grids, SpaceDim , IntVect::Unit);
  if (a_crseDiffusivityPtr != NULL)
    {
      FineInterp  interpolator(a_grids, 1, a_nRefCrse, a_domain);
      interpolator.interpToFine(a_cellDiffusivity, *a_crseDiffusivityPtr);
    }

  for (DataIterator dit(a_grids);dit.ok();++dit)
    {
      FArrayBox& D = a_cellDiffusivity[dit];
      D.setVal(0.0);
      FluxBox& uvavg = vertAverageVelFace[dit];
      const FArrayBox& s = a_coordSys.getSurfaceHeight()[dit];

      //set D = \bar{u'}H / grad(H)

      FORT_L1L2COMPUTEDIFFUSIVITY(CHF_FRA1(D,0),
				  CHF_CONST_FRA1(uvavg[0],0),
				  CHF_CONST_FRA1(uvavg[1],0),
				  CHF_CONST_FRA1(H[dit],0),
				  CHF_CONST_FRA1(H[dit],0),
				  CHF_CONST_INT(SpaceDim),
				  CHF_CONST_REAL(dx[0]),
				  CHF_CONST_REAL(dx[1]),
				  CHF_BOX(a_grids[dit]));
      //set D = \bar{u'}H / grad(s)
      FORT_L1L2COMPUTEDIFFUSIVITY(CHF_FRA1(D,0),
				  CHF_CONST_FRA1(uvavg[0],0),
				  CHF_CONST_FRA1(uvavg[1],0),
				  CHF_CONST_FRA1(H[dit],0),
				  CHF_CONST_FRA1(s,0),
				  CHF_CONST_INT(SpaceDim),
				  CHF_CONST_REAL(dx[0]),
				  CHF_CONST_REAL(dx[1]),
				  CHF_BOX(a_grids[dit]));



      CH_assert(D.norm(0) < HUGE_NORM*HUGE_NORM);
    }
  a_cellDiffusivity.exchange();
  //CellToEdge(cellDiffusivity, a_faceDiffusivity);
  for (DataIterator dit(a_grids);dit.ok();++dit)
    {
      FArrayBox& Dcell = a_cellDiffusivity[dit];
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  FArrayBox& Dface = a_faceDiffusivity[dit][dir];
	  Box fbox = a_grids[dit].surroundingNodes(dir);
	  Real tiny = 1.0/HUGE_NORM;
	  FORT_L1L2CELLTOFACEHARMONIC(CHF_FRA1(Dface,0),
				      CHF_CONST_FRA1(Dcell,0),
				      CHF_CONST_INT(dir),
				      CHF_CONST_REAL(tiny),
				      CHF_BOX(fbox));
	  CH_assert(Dface.norm(0) < HUGE_NORM*HUGE_NORM);

	 
	  //modify grounded ice advection and floating ice diffusion 
          //to account for the relationships between surface, topography and thickness
	  { 
	    const FArrayBox& b = a_coordSys.getTopography()[dit];
	    const BaseFab<int>& mask = a_coordSys.getFloatingMask()[dit];
	    Real ffactor = (1.0 - a_coordSys.iceDensity()/a_coordSys.waterDensity());


	    for (BoxIterator bit(fbox);bit.ok();++bit)
	      {
		const IntVect& iv = bit();
		const IntVect& ivm = iv - BASISV(dir);
	
		if (mask(ivm)==FLOATINGMASKVAL || mask(iv)==FLOATINGMASKVAL)
		  {
		    Dface(iv)*=ffactor;
		  }
		else if (mask(ivm)==GROUNDEDMASKVAL || mask(iv)==GROUNDEDMASKVAL)
		  {
		    Real db = (b(iv)-b(ivm))/dx[dir];
		    Real ua = - 2.0 * db * Dface(iv) / (H[dit](iv) + H[dit](ivm));
		    a_faceVelAdvection[dit][dir](iv) += ua;
		  }
	      }
	  }


	  a_faceVelTotal[dit][dir] += vertAverageVelFace[dit][dir];			      
	}
    }


 
}



#endif
#include "NamespaceFooter.H"
