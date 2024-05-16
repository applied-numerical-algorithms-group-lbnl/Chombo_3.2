#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// glfaces.cpp
// read in a bisicles plotfile and dump out a list of x,y coordinates 
// which sit on the grounding line.
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "JFNKSolver.H"
#include "L1L2ConstitutiveRelation.H"
#include "BasalFrictionRelation.H"
#include "IceThicknessIBC.H"
#include "BasicThicknessIBC.H"
#include "FortranInterfaceIBC.H"
#include "VieliPayneIBC.H"
#include "MarineIBC.H"
#include "HumpIBC.H"
#include "LevelDataIBC.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"

#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif

int main(int argc, char* argv[]) {

#ifdef CH_MPI
#warning "Parallel glfaces might well not work";
  MPI_Init(&argc, &argv);
#endif

{ // Begin nested scope

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int rank, number_procs;
#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank=0;
  number_procs=1;
#endif
  
  if(argc != 3) 
    { 
      std::cerr << " usage: " << argv[0] << " <plot_file> <config_file> " << std::endl; 
      exit(0); 
    }
  char* plot_file = argv[1];
  char* in_file = argv[2];
  
  ParmParse pp(0,0,NULL,in_file);
  ParmParse pp2("main");
  
  Box domainBox;
  
  
  Vector<std::string> name;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout > grids;
  Vector<int > ratio;
  int numLevels;
  
  
  Real dt ,crseDx, time;
  ReadAMRHierarchyHDF5(std::string(plot_file), grids, data, name , 
                       domainBox, crseDx, dt, time, ratio, numLevels);
  
  
  Vector<ProblemDomain> domain(numLevels,domainBox);
  Vector<RealVect> vdx(numLevels,RealVect::Unit*crseDx);
  Vector<Real> dx(numLevels,crseDx);
  for (int lev=1;lev<numLevels;++lev)
    {
      dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
      vdx[lev] = vdx[lev-1] / Real(ratio[lev-1]);
      domain[lev] = domain[lev-1];
      domain[lev].refine(ratio[lev-1]);
    }
  
  //build valid region masks
  Vector<LevelData<BaseFab<int> >* > mask(numLevels,NULL);
  for (int lev = numLevels - 1; lev >= 0; lev--)
    {
      mask[lev] = new LevelData<BaseFab<int> >(grids[lev],1,IntVect::Zero);
      for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
        {
          (*mask[lev])[dit].setVal(1);
	  
          if (lev < numLevels - 1)
            {
              for (DataIterator fit(grids[lev+1]); fit.ok(); ++fit)
                {
                  Box covered = grids[lev+1][fit];
                  covered.coarsen(ratio[lev]);
                  covered &= grids[lev][dit];
                  if (!covered.isEmpty())
                    {
                      (*mask[lev])[dit].setVal(0,covered,0,1);
                    }
                }
            }
        }	
    }
  
  
  
  //constitutive relation and rate factor
  ConstitutiveRelation* constRelPtr = ConstitutiveRelation::parse("main");
  if (constRelPtr == NULL)
    {
      MayDay::Error("undefined constitutiveRelation in inputs");
    }
  
  Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
  {
    ParmParse ppc("constants");
    ppc.query("seconds_per_unit_time",seconds_per_unit_time);
  }
  std::string rateFactorType = "constRate";
  pp2.query("rateFactor", rateFactorType);
  RateFactor* rateFactorPtr; 
  if (rateFactorType == "constRate")
    {
      ParmParse crPP("constRate");
      Real A = 9.2e-18 * seconds_per_unit_time/SECONDS_PER_TROPICAL_YEAR;
      crPP.query("A", A);
      ConstantRateFactor* crf = new ConstantRateFactor(A);
      rateFactorPtr = static_cast<RateFactor*>(crf);
    }
  else if (rateFactorType == "arrheniusRate")
    {
      ArrheniusRateFactor* arf = new ArrheniusRateFactor(seconds_per_unit_time);
      ParmParse arPP("ArrheniusRate");
      rateFactorPtr = static_cast<RateFactor*>(arf);
    }
  else if (rateFactorType == "patersonRate")
    {
      PatersonRateFactor* prf =  new PatersonRateFactor(seconds_per_unit_time);
      ParmParse arPP("PatersonRate");
      rateFactorPtr = static_cast<RateFactor*>(prf);
    }
  else if (rateFactorType == "zwingerRate")
    {
      ZwingerRateFactor* zrf = new ZwingerRateFactor(seconds_per_unit_time);
      ParmParse arPP("ZwingerRate");
      rateFactorPtr = static_cast<RateFactor*>(zrf);
    }
  
  
  //basal friction relation
  BasalFrictionRelation* basalFrictionRelationPtr 
    = BasalFrictionRelation::parse("main",0);
  
  //we only need this for the velocity boundary condition
  IceThicknessIBC* thicknessIBC;
  std::string problem_type;
  ParmParse geomPP("geometry");
  geomPP.get("problem_type", problem_type);
  if (problem_type == "basic")
    {
      thicknessIBC = new BasicThicknessIBC;
    }
  else if (problem_type == "VieliPayne")
    {
      VieliPayneIBC* ibcPtr = new VieliPayneIBC;     
      thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
    }
  else if (problem_type == "marineIceSheet")
    {
      MarineIBC* ibcPtr = new MarineIBC;
      thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
    }
  else if (problem_type == "hump")
    {
      HumpIBC* ibcPtr = new HumpIBC;
      ParmParse humpPP("hump");
      thicknessIBC = static_cast<IceThicknessIBC*>(ibcPtr);
    }
  else if (problem_type == "LevelData")
    {
      LevelDataIBC* ptr = new LevelDataIBC();
      thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
    }
  else if (problem_type =="fortran")
    {
      FortranInterfaceIBC* ptr = new FortranInterfaceIBC;
      thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
    }
#ifdef HAVE_PYTHON
  else if (problem_type == "Python")
    {	 
      ParmParse pyPP("PythonIBC");
      std::string module;
      pyPP.get("module",module);
      std::string thckFuncName = "thickness";
      pyPP.query("thicknessFunction",thckFuncName);
      std::string topgFuncName = "topography";
      pyPP.query("topographyFunction",topgFuncName);
      std::string rhsFuncName = "";
      pyPP.query("RHSFunction",rhsFuncName);
      std::string faceVelFuncName = "";
      pyPP.query("faceVelFunction",faceVelFuncName);
      PythonInterface::PythonIBC* ptr = new PythonInterface::PythonIBC
        (module, thckFuncName, topgFuncName, rhsFuncName,faceVelFuncName);
      thicknessIBC = static_cast<IceThicknessIBC*>( ptr);
    }
#endif
  else 
    {
      MayDay::Error("bad problem type");
    }
  
  ParmParse ppCon("constants");
  // set reasonable defaults for constants
  // (as of 8/6/20 using the defaults in SigmaCS)
  Real iceDensity = 910.0;
  ppCon.query("ice_density",iceDensity);
  Real seaWaterDensity = 1028.0;
  ppCon.query("sea_water_density",seaWaterDensity);
  Real gravity = 9.81;
  ppCon.query("gravity",gravity);
  
  ParmParse ppAmr("amr");
  Vector<int> ancells(3);
  ppAmr.getarr("num_cells", ancells, 0, ancells.size());
  
  //periodicity information is not stored in plot files
  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; dir++)
    is_periodic[dir] = false;
  Vector<int> is_periodic_int(SpaceDim, 0);
  ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++) 
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }
  
  
  int nLayers = ancells[2];
  Vector<Real> faceSigma(nLayers+1);
  Real dsigma = 1.0 / Real(nLayers);
  for (unsigned int l = 0; l < faceSigma.size(); ++l)
    faceSigma[l] = dsigma * (Real(l));
  ppAmr.queryarr("sigma",faceSigma,0,faceSigma.size());
  
  //extract the topography, basal friction, thickness and velocity fields
  Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
  Vector<LevelData<FArrayBox>* > C(numLevels,NULL); // basal friction coefficient
  Vector<LevelData<FArrayBox>* > C0(numLevels,NULL); // basal friction coefficien
  Vector<LevelData<FArrayBox>* > vel(numLevels,NULL); // velocity field
  Vector<LevelData<FArrayBox>* > cellA(numLevels,NULL); // cell-centered rate factor
  Vector<LevelData<FluxBox>* > faceA(numLevels,NULL); // face-centered rate factor
  Vector<LevelData<FluxBox>* > faceVT(numLevels,NULL); // face-centered viscous tensor
  
  IntVect sigmaCSGhost(2*IntVect::Unit);
  
  for (int lev = 0; lev < numLevels; lev++)
    {
      coords[lev] = RefCountedPtr<LevelSigmaCS> 
        (new LevelSigmaCS(grids[lev], vdx[lev], sigmaCSGhost));
      coords[lev]->setIceDensity(iceDensity);
      coords[lev]->setWaterDensity(seaWaterDensity);
      coords[lev]->setGravity(gravity);
      coords[lev]->setFaceSigma(faceSigma);
      
      
      C[lev] = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      C0[lev] = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      vel[lev] = new LevelData<FArrayBox>(grids[lev], SpaceDim, IntVect::Unit);
      cellA[lev] = new LevelData<FArrayBox>(grids[lev], faceSigma.size()-1, IntVect::Unit);
      faceA[lev] = new LevelData<FluxBox>(grids[lev], faceSigma.size()-1, IntVect::Zero);
      faceVT[lev] = new LevelData<FluxBox>(grids[lev], SpaceDim, IntVect::Zero);
      if (lev > 0)
        {
          coords[lev]->interpFromCoarse(*coords[lev-1],ratio[lev-1]);
        }
      data[lev]->exchange();
      for (int j = 0; j < name.size(); j++)
        {
          if (name[j] == "Z_base")
            {
              data[lev]->copyTo(Interval(j,j),coords[lev]->getTopography(),Interval(0,0));
            }
          else if (name[j] == "thickness")
            {
              data[lev]->copyTo(Interval(j,j),coords[lev]->getH(),Interval(0,0));
            }
          else if (name[j] == "basal_friction")
            {
              data[lev]->copyTo(Interval(j,j),*C[lev],Interval(0,0));
            }
          else if (name[j] == "xVel")
            {
              data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(0,0));
            }
          else if (name[j] == "yVel")
            {
              data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(1,1));
            }
          
          for (DataIterator dit(grids[lev]);dit.ok();++dit)
            {
              (*cellA[lev])[dit].setVal(3.1536e-18); //yikes....
              (*C0[lev])[dit].setVal(0.0);
            }
        }
      
      thicknessIBC->setGeometryBCs(*coords[lev],domain[lev],vdx[lev],time,dt);
      
      if (lev > 0)
        {
          coords[lev]->recomputeGeometry(coords[lev-1],ratio[lev-1]);
        }
      else
        {
          coords[lev]->recomputeGeometry(NULL,0);
        }
      CellToEdge(*cellA[lev],*faceA[lev]);
    }
  
  //these parameters don't matter because we don't solve anything. 
  Real vtopSafety = 1.0;
  int vtopRelaxMinIter = 4;
  Real vtopRelaxTol = 1.0;
  Real muMin = 0.0;
  Real muMax = 1.23456789e+300;
  IceNonlinearViscousTensor state(grids, ratio, domain, vdx, coords, vel, C, C0,
                                  numLevels-1, *constRelPtr, *basalFrictionRelationPtr, 
                                  *thicknessIBC, cellA, faceA, time, 
                                  vtopSafety, vtopRelaxMinIter, vtopRelaxTol, 
                                  muMin, muMax);
  state.setState(vel);
  state.computeViscousTensorFace(faceVT);
  
  for (int lev = 0; lev < numLevels; lev++)
    {
      for (DataIterator dit(grids[lev]);dit.ok();++dit)
        {
          //const FArrayBox& thisC = (*C[lev])[dit];
          
          const BaseFab<int>& validMask = (*mask[lev])[dit];
          const BaseFab<int>& floatingMask = coords[lev]->getFloatingMask()[dit];
          const FArrayBox& thk = coords[lev]->getH()[dit];
          const FArrayBox& topg = coords[lev]->getTopography()[dit];
          const FArrayBox& u = (*vel[lev])[dit];
          const Real& rhoi = coords[lev]->iceDensity();
          const Real& rhow = coords[lev]->waterDensity();
          const Real& gravity = coords[lev]->gravity();
          const Real flotationThicknessFactor = -rhow/rhoi;
          
          for (int fdir = 0; fdir < SpaceDim; ++fdir)
            {
              const FArrayBox& fvt =  (*faceVT[lev])[dit][fdir];
              const FArrayBox& fthk =  coords[lev]->getFaceH()[dit][fdir];
              
              Box b = grids[lev][dit];
              //b.grow(fdir,-1);
              for (BoxIterator bit(b);bit.ok();++bit)
                {
                  IntVect iv = bit();
                  if (validMask(iv) == 1 && floatingMask(iv) == GROUNDEDMASKVAL)
                    {
                      
                      for (int sign = -1; sign <=1; sign+=2)
                        {
                          const IntVect ivp = iv+sign*BASISV(fdir);
                          if ((floatingMask(ivp) == FLOATINGMASKVAL) || (floatingMask(ivp) == OPENSEAMASKVAL))
                            {
                              //face centred point between iv amd ivp
                              int off = (sign<0)?0:1;
                              IntVect ivf = iv + off*BASISV(fdir); 
                              
                              pout() << "gl " << sign << " " << ((fdir==0)?"x ":"y ") << lev << " ";
                              //grid index on face centred grids
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << ivf[dir] << " ";
                              //grid index of the grounded cell
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << iv[dir] << " ";
                              //grid index of the floating cell
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << ivp[dir] << " ";
                              //(x,y) coords at the face
                              for (int dir = 0; dir < SpaceDim; dir++)
                                {
                                  if (dir == fdir)
                                    pout() << Real(ivf[dir])*dx[lev] << " ";
                                  else
                                    pout() << (Real(ivf[dir]) + 0.5)*dx[lev] << " ";
                                }
                              //(x,y) coords of the grounded cell
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << (Real(iv[dir]) + 0.5)*dx[lev] << " ";
                              //(x,y) coords of the floating cell
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << (Real(iv[dir]) + 0.5)*dx[lev] << " ";
                              //time
                              pout() << time << " ";
                              // thickness by averaging
                              pout() << fthk(ivf) << " ";
                              // thickness by flotation
                              Real ftopg = 0.5*(topg(iv)+topg(ivp));
                              pout() << flotationThicknessFactor * ftopg << " ";
                              // velocity at the grounding line
                              for (int dir = 0; dir < SpaceDim; dir++)
                                pout() << 0.5*(u(iv,dir)+u(ivp,dir)) << " ";
                              //viscous tensor components;
                              //for (int dir = 0; dir < SpaceDim; dir++)
                              //  pout() << fvt(ivf,dir) << " ";
                              pout() << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
  
  for (int lev = 0; lev < numLevels; lev++)
    {
      if (C[lev] != NULL)
        {
          delete C[lev];
          C[lev] = NULL;
        }
      
      if (C0[lev] != NULL)
        {
          delete C0[lev];
          C0[lev] = NULL;
        }
      
      if (vel[lev] != NULL)
        {
          delete C0[lev];
          vel[lev] = NULL;
        }
      
      if (cellA[lev] != NULL)
        {
          delete cellA[lev];
          cellA[lev] = NULL;
        }
      
      if (faceA[lev] != NULL)
        {
          delete faceA[lev];
          faceA[lev] = NULL;
        }
      
      if (faceVT[lev] != NULL)
        {
          delete faceVT[lev];
          faceVT[lev] = NULL;
        }
      
      if (mask[lev] != NULL)
        {
          delete mask[lev];
          mask[lev] = NULL;
        }
      
      if (data[lev] != NULL)
        {
          delete data[lev];
        }
    }
  if (rateFactorPtr != NULL)
    {
      delete rateFactorPtr;
      rateFactorPtr = NULL;
    }
  
  if (basalFrictionRelationPtr != NULL)
    {
      delete basalFrictionRelationPtr;
    }
  
  if (thicknessIBC != NULL)
    {
      delete thicknessIBC;
    }
  
 }  // end nested scope
 CH_TIMER_REPORT();
 
#ifdef CH_MPI
 MPI_Finalize();
#endif
 
 return 0;
}
