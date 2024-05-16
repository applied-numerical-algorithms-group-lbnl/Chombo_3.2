#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "NonlinearViscousTensor.H"
#include "QuadCFInterp.H"
#include "CornerCopier.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "BisiclesF_F.H"
#include "LinearSolver.H"
#include "RelaxSolver.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "CGSolver.H"
#include "IceUtility.H"
#include <sstream>
#include "NamespaceHeader.H"


IceNonlinearViscousTensor::~IceNonlinearViscousTensor()
{
  delete m_bcPtr;
}

NonlinearViscousTensor* 
IceNonlinearViscousTensor::newNonlinearViscousTensor()
{
  IceNonlinearViscousTensor* ptr = new IceNonlinearViscousTensor(*this);
  ptr->setFaceViscCoef(m_muCoef);
  return static_cast<NonlinearViscousTensor*>(ptr);
}

void IceNonlinearViscousTensor::setupCoeffs()
{
  m_mu.resize(m_finestLevel + 1);
  m_lambda.resize(m_finestLevel + 1);
  m_alpha.resize(m_finestLevel + 1);
  m_muCoef.resize(m_finestLevel + 1);

  for (int lev =0; lev <= m_finestLevel; ++lev)
    {
      DisjointBoxLayout levelGrids = m_grids[lev];

      m_mu[lev] = RefCountedPtr<LevelData<FluxBox> >
	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Unit));
      
      m_muCoef[lev] = NULL;
      
      m_lambda[lev] = RefCountedPtr<LevelData<FluxBox> >
	(new LevelData<FluxBox>(levelGrids, 1, IntVect::Unit));	
      
      m_alpha[lev] = RefCountedPtr<LevelData<FArrayBox> >
	(new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit));

      // (DFM - 10/29/13) sending these values into the constructor 
      // uninitialized is giving valgrind fits. Initialize everything 
      // to bogus values here
      bool initializeValues = true;
      if (initializeValues)
        {
          // negative viscosities and friction should be suitably unphysical
          Real bogusValue = -1;
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*m_mu[lev])[dit].setVal(bogusValue);
              (*m_lambda[lev])[dit].setVal(bogusValue);
              (*m_alpha[lev])[dit].setVal(bogusValue);
            }
        }
    }

   Real alpha = -1.0;
   Real beta = 1.0;
   
   // for the moment, at least, this only works for dx = dy:
   if (SpaceDim > 1) CH_assert(m_dxs[0][0] == m_dxs[0][1]);

   m_velSolveBC = m_bcPtr->velocitySolveBC();
   m_opFactoryPtr = 
     RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
     (new ViscousTensorOpFactory(m_grids, m_mu, m_lambda, m_alpha, alpha, 
				 beta, m_refRatio , m_domains[0], m_dxs[0][0], 
				 m_velSolveBC, m_vtopSafety, m_vtopRelaxTol, 
				 m_vtopRelaxMinIter));

}


IceNonlinearViscousTensor::IceNonlinearViscousTensor
(const IceNonlinearViscousTensor& a)
  :m_u(a.m_u), m_C(a.m_C), m_C0(a.m_C0), m_grids(a.m_grids), m_refRatio(a.m_refRatio),
   m_domains(a.m_domains), m_dxs(a.m_dxs),m_finestLevel(a.m_finestLevel), 
   m_coordSys(a.m_coordSys), m_constRelPtr(a.m_constRelPtr), 
   m_basalFrictionRelPtr(a.m_basalFrictionRelPtr), m_bcPtr(a.m_bcPtr->new_thicknessIBC()),
   m_A(a.m_A), m_faceA(a.m_faceA) , m_time(a.m_time), m_vtopSafety(a.m_vtopSafety),
   m_vtopRelaxMinIter(a.m_vtopRelaxMinIter),m_vtopRelaxTol(a.m_vtopRelaxTol),
   m_muMin(a.m_muMin),m_muMax(a.m_muMax),m_scale(a.m_scale),
   m_artificialDragCoef(a.m_artificialDragCoef),
   m_artificialDragPower(a.m_artificialDragPower)
{
  setupCoeffs();
}


IceNonlinearViscousTensor::IceNonlinearViscousTensor
(const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio,
 const Vector<ProblemDomain>& a_domains,
 const Vector<RealVect>& a_dxs,
 const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 const Vector<LevelData<FArrayBox>*>& a_u,
 const Vector<LevelData<FArrayBox>*>& a_C,
 const Vector<LevelData<FArrayBox>*>& a_C0,
 const int a_finestLevel,
 const ConstitutiveRelation& a_constRel,
 const BasalFrictionRelation& a_basalFrictionRel,
 IceThicknessIBC& a_bc,
 const Vector<LevelData<FArrayBox>*>& a_A,
 const Vector<LevelData<FluxBox>*>& a_faceA,
 Real a_time,
 Real a_vtopSafety,
 int a_vtopRelaxMinIter,
 Real a_vtopRelaxTol,
 Real a_muMin ,
 Real a_muMax,
 Real a_scale, Real a_artificialDragCoef, Real a_artificialDragPower)
:m_u(a_u), m_C(a_C), m_C0(a_C0), m_grids(a_grids), m_refRatio(a_refRatio),
 m_domains(a_domains), m_dxs(a_dxs),m_finestLevel(a_finestLevel), 
 m_coordSys(a_coordSys), m_constRelPtr(&a_constRel), 
 m_basalFrictionRelPtr(&a_basalFrictionRel), m_bcPtr(a_bc.new_thicknessIBC()),
 m_A(a_A), m_faceA(a_faceA) , m_time(a_time), m_vtopSafety(a_vtopSafety),
 m_vtopRelaxMinIter(a_vtopRelaxMinIter),m_vtopRelaxTol(a_vtopRelaxTol),
 m_muMin(a_muMin),m_muMax(a_muMax),m_scale(a_scale),
 m_artificialDragCoef(a_artificialDragCoef),m_artificialDragPower(a_artificialDragPower)
{
  setupCoeffs();
}


void IceNonlinearViscousTensor::computeViscousTensorFace(const Vector<LevelData<FluxBox>*>& a_viscousTensor)
{
  if (m_opFactoryPtr == NULL)
    MayDay::Error("IceNonlinearViscousTensor::computeViscousTensorFace missing m_opFactoryPtr");
  
  ViscousTensorOpFactory* vtopf = static_cast<ViscousTensorOpFactory*>(&(*m_opFactoryPtr));

  //compute cell centered grad(vel)
  Vector<LevelData<FArrayBox>* > gradVel(m_finestLevel + 1, NULL);
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      const DisjointBoxLayout& grids = m_grids[lev];
      gradVel[lev] = new LevelData<FArrayBox>(grids , SpaceDim*SpaceDim, IntVect::Unit);
      LevelData<FArrayBox>& vel = *m_u[lev];
      LevelData<FArrayBox>& grad = *gradVel[lev];
      
      if (lev > 0)
	{
	  const DisjointBoxLayout& crseGrids = m_grids[lev-1];
	  const ProblemDomain& pd = grids.physDomain();
	  TensorCFInterp cf(grids,&crseGrids,m_dxs[lev][0],m_refRatio[lev-1],SpaceDim,pd);
	  cf.coarseFineInterp(vel,grad,*m_u[lev-1]);
	}
      ViscousTensorOp* vtop = vtopf->AMRnewOp(m_domains[lev]);
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	{
	  vtop->cellGrad(grad[dit],vel[dit],grids[dit]);
	}
      grad.exchange();
      delete vtop;
    }
  
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      ViscousTensorOp* vtop = vtopf->AMRnewOp(m_domains[lev]);
      const DisjointBoxLayout& grids = m_grids[lev];
      LevelData<FluxBox>& flux = *a_viscousTensor[lev];
      LevelData<FArrayBox>& vel = *m_u[lev];
      LevelData<FArrayBox>& grad = *gradVel[lev];
      
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      const FArrayBox& muFace  =  (*m_mu[lev])[dit][dir];
	      const FArrayBox& lambdaFace  = (*m_lambda[lev])[dit][dir];
	      Box faceBox =  surroundingNodes(grids[dit],dir);
	      Box domFaceBox = surroundingNodes(m_domains[lev].domainBox(), dir);
	      faceBox &= domFaceBox;
	      FArrayBox tmpFlux(faceBox,SpaceDim); //  vtop->getFlux tinkers with the FAB box
	      vtop->getFlux(tmpFlux, vel[dit], grad[dit] , muFace, lambdaFace, faceBox, dir, 1);
	      flux[dit][dir].setVal(0.0);
	      flux[dit][dir].copy(tmpFlux,faceBox);
	    }
	}
      delete vtop;
    }
  
  for (int lev=0; lev <= m_finestLevel; lev++)
    {
      if (gradVel[lev] != NULL)
	{
	  delete gradVel[lev];
	  gradVel[lev] = NULL; 
	}
    }
}


//store a_u and set the coeffients mu(u), lambda(u) = 2*mu(u) and alpha(u)
//in L[u] = div( mu * (grad(u) + grad(u)^T) + lambda * div(u)*I) - alpha*u
void IceNonlinearViscousTensor::setState(const Vector<LevelData<FArrayBox>*>& a_u)
{
  m_u = a_u;
  
  CH_TIME("IceNonlinearViscousTensor::setState");

  pout() << "IceNonlinearViscousTensor::setState (recomputes viscous tensor coeffs)" << endl;
  
  for (int lev=0; lev <= m_finestLevel ; lev++)
    {
      ProblemDomain levelDomain = m_domains[lev];
      const DisjointBoxLayout& levelGrids = m_grids[lev];
      LevelData<FArrayBox>& levelVel = *a_u[lev];
      LevelData<FArrayBox>* crseVelPtr = NULL;
      int nRefCrse = -1;
      if (lev > 0) 
        {
          crseVelPtr = a_u[lev-1];
          nRefCrse = m_refRatio[lev-1];
        }
      LevelData<FluxBox>& levelMu = *m_mu[lev];
     
      LevelSigmaCS& levelCoords = *m_coordSys[lev];
      const LevelData<FArrayBox>& levelC = *m_C[lev];
      LevelData<FArrayBox>& levelAlpha = *m_alpha[lev];

      // // first thing, if there is a finer level, average-down
      // // the current velocity field
      if (lev < m_finestLevel )
        {
          LevelData<FArrayBox>& finerVel = *a_u[lev+1];

          CoarseAverage averager(finerVel.getBoxes(),
                                 levelGrids,
                                 finerVel.nComp(),
                                 m_refRatio[lev]);

          averager.averageToCoarse(levelVel, finerVel);
        }

      // just in case, add an exchange here
      levelVel.exchange();
      
      // first set BC's on vel
      m_bcPtr->velocityGhostBC(levelVel,
                               levelCoords,
                               levelDomain, m_time);
      

      DataIterator dit = levelGrids.dataIterator();
      dit.reset();

      if (lev > 0) 
        {
          QuadCFInterp qcfi(levelGrids, &m_grids[lev-1],
                            m_dxs[lev][0], m_refRatio[lev-1], 
                            SpaceDim, levelDomain);
          qcfi.coarseFineInterp(levelVel, *a_u[lev-1]);
        }

      //slc : qcfi.coarseFineInterp fills the edges of lev > 0 cells
      //but not the corners. We need them filled to compute the
      //rate-of-strain invariant, so here is a bodge for now
      if (SpaceDim == 2)
      	{
      	  for (dit.begin(); dit.ok(); ++dit)
      	    {
      	      Box sbox = levelVel[dit].box();
      	      sbox.grow(-1);
      	      FORT_EXTRAPCORNER2D(CHF_FRA(levelVel[dit]),
      				  CHF_BOX(sbox));
      	    }
	  
      	}


      // actually need to use a cornerCopier, too...
      CornerCopier cornerCopier(levelGrids, levelGrids, 
                                levelDomain,levelVel.ghostVect(),
                                true);
      levelVel.exchange(cornerCopier);
				
      (*m_constRelPtr).computeFaceMu(levelMu,
                                     levelVel, m_scale, 
                                     crseVelPtr,
                                     nRefCrse,
                                     *m_faceA[lev],
                                     levelCoords,
				     levelDomain,
                                     IntVect::Zero);

      // now limit and  multiply by ice thickness H
      const LevelData<FluxBox>& faceH = levelCoords.getFaceH();
      LevelData<FluxBox>* muCoefPtr = m_muCoef[lev];
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
	      FArrayBox& thisMu = levelMu[dit][dir];
	      const Box& box = thisMu.box();
	       
	      FORT_MAXFAB1(CHF_FRA(thisMu),
	       		  CHF_CONST_REAL(m_muMin),
	       		  CHF_BOX(box));

	      if (muCoefPtr!=NULL)
		{
		  levelMu[dit][dir].mult( (*muCoefPtr)[dit][dir] );
		}

	      Box facebox = levelGrids[dit].surroundingNodes(dir);
	      CH_assert(levelMu[dit][dir].min(facebox) >= 0.0);

              levelMu[dit][dir].mult(faceH[dit][dir],
                                    levelMu[dit][dir].box(),0,0,1);
	    
	      FORT_MINFAB1(CHF_FRA(thisMu),
	       		  CHF_CONST_REAL(m_muMax),
	       		  CHF_BOX(box));

	    }
        
	  // also update alpha
          const Box& gridBox = levelGrids[dit];

	  m_basalFrictionRelPtr->computeAlpha
		      (levelAlpha[dit], levelVel[dit], levelC[dit] , m_scale, 	
	     levelCoords, dit, lev, gridBox);
	  levelAlpha[dit] += (*m_C0[lev])[dit];

	  //artificial drag
	  if (m_artificialDragCoef > 0.0)
	    {
	      const Real& a = m_artificialDragCoef;
	      Real asq = a*a; 
	      Real p = 0.5*(m_artificialDragPower-1.0);
	      const FArrayBox& u = levelVel[dit];
	      for ( BoxIterator bit(gridBox); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  Real usq = D_TERM(u(iv,0)*u(iv,0), + u(iv,1)*u(iv,1), + u(iv,2)*u(iv,2));
		  levelAlpha[dit](iv) += std::pow(asq * usq,p)*a; // alpha = a^(n-1) |u|^(n-1) a  => |alpha*u| = a^n |u|^n
 		}
	    }
	  CH_assert(levelAlpha[dit].min(gridBox) >= 0.0);

	  
	} // end loop over grids		
    }

  //coarse average mu and alpha
  for (int lev = m_finestLevel; lev > 0 ; lev--)
  {
      CoarseAverage averageOp(m_grids[lev],1,m_refRatio[lev-1]);
      averageOp.averageToCoarse(*m_alpha[lev-1], *m_alpha[lev]);

      CoarseAverageFace averageOpFace(m_grids[lev],1,m_refRatio[lev-1]);
      averageOpFace.averageToCoarse(*m_mu[lev-1], *m_mu[lev] );
  }
  // lambda = 2*mu
  for (int lev=0; lev <= m_finestLevel ; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_grids[lev];
      LevelData<FluxBox>& levelLambda = *m_lambda[lev];
      LevelData<FluxBox>& levelMu = *m_mu[lev];
      LevelData<FArrayBox>& levelAlpha = *m_alpha[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{

#if CH_SPACEDIM==2
	  {
	    Real mu0 = 1.0;
	    Real C0 = 1.0;
	    
	    FORT_ENFORCEWELLPOSEDCELL
	      (CHF_FRA1(levelAlpha[dit],0),
	       CHF_FRA1(levelMu[dit][0],0),
	       CHF_FRA1(levelMu[dit][1],0),
	       CHF_CONST_REAL(mu0),
	       CHF_CONST_REAL(C0),
	       CHF_BOX(levelGrids[dit]));

	  }
#endif


	  FluxBox& lambda = levelLambda[dit];
	  for (int dir=0; dir<SpaceDim; dir++)
	    {
	      lambda[dir].copy(levelMu[dit][dir]);
	      lambda[dir] *= 2.0;
	    }
	} // end loop over grids
    } // end loop over levels

}





