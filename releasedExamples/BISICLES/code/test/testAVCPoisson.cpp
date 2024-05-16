#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;


#include "AVCAMRPoissonOp.H"
#include "BiCGStabSolver.H"
#include "LoadBalance.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "BoxIterator.H"
int testAVCPoisson();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  
  int status = testAVCPoisson();

  pout() << "AVCPoisson test";

  if( status == 0 )
    pout() << " passed." << endl ;
  else
    pout() << " failed with return code "
         << status << endl ;

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int testAVCPoisson()
{
  int rc = -1;
  int nCrse = 64;
  int max_box = 32;
  int block_factor = 1;
  IntVect lo = IntVect::Zero;
  IntVect hi = (nCrse-1) * IntVect::Unit;
  Box crseDomBox(lo,hi);
  Real crseDx = 1.0/Real(nCrse);

  int finest_level = 1;
  Vector<ProblemDomain> domain(finest_level + 1);
  Vector<DisjointBoxLayout> grids(finest_level + 1);
  Vector<int> ratio(finest_level + 1,2);
  Vector<Real> dx(finest_level + 1,-crseDx);
  domain[0] = ProblemDomain(crseDomBox);
  {
    //0-level set up
    Vector<Box> boxes;
    domainSplit(domain[0],boxes,max_box, block_factor);
    Vector<int> proc(boxes.size());
    LoadBalance(proc,boxes);
    grids[0] = DisjointBoxLayout(boxes, proc, domain[0]);
    dx[0] = crseDx;
  }
  for (int lev = 1; lev < finest_level + 1; lev++)
    {
      Vector<Box> boxes;
      IntVect reflo = IntVect::Zero;
      IntVect refhi = (nCrse-1) * IntVect::Unit;
      Box refDomBox(lo,hi);
      domainSplit(ProblemDomain(refDomBox),boxes,max_box, block_factor);
      Vector<int> proc(boxes.size());
      LoadBalance(proc,boxes);
      domain[lev] = refine(domain[lev-1],ratio[lev-1]);
      grids[lev] = DisjointBoxLayout(boxes, proc, domain[lev]);
      dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
    }
      
  Real alpha = 1.0;
  Real beta  = 0.1;
  Real gamma = beta * 10.0;
  Real nx = std::sqrt(1.0/2.0);
  Real ny = std::sqrt(1.0 - nx*nx);
  Real amplitude = 1.0 / std::pow(0.5, 4);


  Vector<LevelData<FArrayBox>* > phi(finest_level + 1, NULL);
  Vector<LevelData<FArrayBox>* > out(finest_level + 1, NULL);
  Vector<LevelData<FArrayBox>* > rhs(finest_level + 1, NULL);
  Vector< RefCountedPtr<LevelData<FArrayBox> > > acoef(finest_level + 1);
  Vector< RefCountedPtr<LevelData<FluxBox> > > bcoef(finest_level + 1);

  //set the rhs, initial guess, and coeffs

  for (int lev = 0; lev < finest_level + 1; lev++)
    {
      phi[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      rhs[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Zero);
      out[lev] = new LevelData<FArrayBox>(grids[lev],3,IntVect::Zero);
      acoef[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(grids[lev],1,IntVect::Zero));
      bcoef[lev] = RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(grids[lev],SpaceDim,IntVect::Zero));

      for (DataIterator dit(grids[lev]);dit.ok();++dit)
	{
	  (*phi[lev])[dit].setVal(0.0);
	  (*acoef[lev])[dit].setVal(alpha);
	  (*bcoef[lev])[dit][0].setVal(beta + gamma * nx * nx ,0);
	  (*bcoef[lev])[dit][0].setVal(gamma * nx * ny,1);
	  (*bcoef[lev])[dit][1].setVal(beta + gamma * nx * nx,1);
	  (*bcoef[lev])[dit][1].setVal(gamma * nx * ny,0);

	  for (BoxIterator bit(grids[lev][dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real x = (Real(iv[0]) + 0.5)*dx[lev];
	      Real y = (Real(iv[1]) + 0.5)*dx[lev];

	      (*rhs[lev])[dit](iv) = amplitude * (
		+ 2.0 * (1.0-x)*x*(ny*ny*gamma+beta)
		+ 2.0 * (1.0-y)*y*(nx*nx*gamma+beta)
		- 2.0 * nx*ny*(x*y-(1.0-x)*y-x*(1.0-y)+(1.0-x)*(1.0-y))*gamma
		+ 1.0 * alpha*(1.0-x)*x*(1.0-y)*y );

	      (*out[lev])[dit](iv,0) =  amplitude * x * (1.0-x) * y * (1.0-y);
	      (*phi[lev])[dit](iv,0) = (*out[lev])[dit](iv,0);
	    }
	}
    }

  

  //Dirichlett boundary conditions, u = 0 on all edges
  BCHolder bc(ConstDiriNeumBC(IntVect::Unit, RealVect::Zero,
			      IntVect::Unit, RealVect::Zero));

  AVCAMRPoissonOpFactory opf;
  opf.define(domain[0],  grids , ratio, crseDx, bc, 1.0 , acoef,  1.0 , bcoef );


  AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  mgSolver.define(domain[0], opf, &bottomSolver, finest_level+1);
  mgSolver.m_eps = 1.0e-30;
  mgSolver.m_normThresh = 1.0e-8;
  mgSolver.m_iterMax = 20;
  mgSolver.m_hang = 0.0;
  int numMGSmooth = 2;
  mgSolver.m_pre = numMGSmooth;
  mgSolver.m_post = numMGSmooth;
  mgSolver.m_bottom = numMGSmooth;
  mgSolver.m_verbosity = 10;
  mgSolver.solve(phi, rhs, finest_level, 0,  false);

  for (int lev = 0; lev < finest_level + 1; lev++)
    {
      phi[lev]->copyTo(Interval(0,0),*out[lev],Interval(1,1));
      rhs[lev]->copyTo(Interval(0,0),*out[lev],Interval(2,2));
    }

  Vector<std::string> names;
  names.push_back("syn");
  names.push_back("sol");
  names.push_back("rhs");

  std::string file("testAVCPoisson.2d.hdf5");
  WriteAMRHierarchyHDF5(file ,grids, out ,names, domain[0].domainBox(),
  			dx[0], 0.0, 0.0, ratio, finest_level + 1);

  rc = 0;
  return rc;
}
