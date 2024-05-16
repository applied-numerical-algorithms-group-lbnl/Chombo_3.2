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
// ppismip6.cpp
// Post-processing for ISMIP6. Reads in a bisicles plotfile and creates a
// new plot file that includes ISMIP6 fields that are not already present
//
// ligroundf  : flux across grounding line.
//              Only for grid cells immediately upstream of the GL
// licalvf    : flux across the calving front (assuming fixed front) 
// sftgrf     : Fraction of grid cell covered by grounded ice
// sftlflf    : Fraction of grid cell covered by grounded ice
// strbasemag : Basal traction
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
#include "PiecewiseLinearFillPatch.H"
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
  
  if(argc != 4) 
    { 
      std::cerr << " usage: " << argv[0] << " <plot_file> <config_file> <post_process_file>" << std::endl; 
      exit(0); 
    }
  char* plot_file = argv[1];
  char* in_file = argv[2];
  char* out_file = argv[3];

  ParmParse pp(0,0,NULL,in_file);
  
  // READ constants
  ParmParse ppc("constants");
  Real seconds_per_unit_time = SECONDS_PER_TROPICAL_YEAR;
  ppc.query("seconds_per_unit_time",seconds_per_unit_time);
  Real iceDensity = 910.0;
  ppc.query("ice_density",iceDensity);
  Real seaWaterDensity = 1028.0;
  ppc.query("sea_water_density",seaWaterDensity);
  Real gravity = 9.81;
  ppc.query("gravity",gravity);
   
  Box domainBox;
  Vector<std::string> name;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout > grids;
  Vector<int > ratio;
  int numLevels;
  Real dt, crseDx, time;

  //Create a vector called newData that will store the output data
  Vector<LevelData<FArrayBox>*> newData;
#define SFTGRF 0
#define SFTFLF 1
#define LIGROUNDF 2
#define LICALVF 3
#define STRBASEMAG 4
#define NEW_DATA_NUM_COMP 5
  Vector<std::string> newNames(NEW_DATA_NUM_COMP);
  newNames[SFTGRF] = "sftgrf"; //grounded ice sheet area fraction
  newNames[SFTFLF] = "sftflf"; //floating ice sheet area fraction
  newNames[LIGROUNDF] = "ligroundf"; //grounding line flux
  newNames[LICALVF] = "licalvf"; // calving front flux
  newNames[STRBASEMAG] = "strbasemag"; // basal drag

  //Read in the HDF5 plotfile
  ReadAMRHierarchyHDF5(std::string(plot_file), grids, data, name, domainBox, crseDx, dt, time, ratio, numLevels);
  
  //resize the vector depending on the number of levels the plotfile has
  //e.g if the plotfile has 4 levels, newData has a size of 4, each NULL valued
  newData.resize(numLevels, NULL);
  
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
  //extract the topography, basal friction, thickness and velocity fields
  Vector<RefCountedPtr<LevelSigmaCS > > coords(numLevels);
 
  Vector<LevelData<FArrayBox>* > dragCoef(numLevels,NULL); // basal drag coefficien  
  Vector<LevelData<FArrayBox>* > vel(numLevels,NULL); // velocity field

  IntVect sigmaCSGhost(2*IntVect::Unit);

  int newDataNumComps = 4;
  
  for (int lev = 0; lev < numLevels; lev++)
    {
      coords[lev] = RefCountedPtr<LevelSigmaCS> (new LevelSigmaCS(grids[lev], vdx[lev], sigmaCSGhost));
      coords[lev]->setIceDensity(iceDensity);
      coords[lev]->setWaterDensity(seaWaterDensity);
      coords[lev]->setGravity(gravity);
      newData[lev] = new LevelData<FArrayBox>(grids[lev], NEW_DATA_NUM_COMP, IntVect::Zero);
      dragCoef[lev] = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
      vel[lev] = new LevelData<FArrayBox>(grids[lev], SpaceDim, IntVect::Unit);

      if (lev > 0)
        {
          coords[lev]->interpFromCoarse(*coords[lev-1],ratio[lev-1]);
	  PiecewiseLinearFillPatch pwlfill (grids[lev],grids[lev-1],SpaceDim,
					   grids[lev-1].physDomain(),ratio[lev-1], 1);
	  pwlfill.fillInterp(*vel[lev],*vel[lev-1],*vel[lev-1],0.0,0,0,SpaceDim);
	  
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
          else if (name[j] == "dragCoef")
            {
              data[lev]->copyTo(Interval(j,j),*dragCoef[lev],Interval(0,0));
            }          
          else if (name[j] == "xVel")
            {
              data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(0,0));
            }
          else if (name[j] == "yVel")
            {
              data[lev]->copyTo(Interval(j,j),*vel[lev],Interval(1,1));
            }
          
        } 
      if (lev > 0)
        {
          coords[lev]->recomputeGeometry(coords[lev-1],ratio[lev-1]);
        }
      else
        {
          coords[lev]->recomputeGeometry(NULL,0);
        }
    }
  

  //on each level
  for (int lev = 0; lev < numLevels; lev++)
    {
      //for each box
      for (DataIterator dit(grids[lev]);dit.ok();++dit)
        {
          //const FArrayBox& thisC = (*C[lev])[dit];
          
          //create thisNewData which is a specific box at a certain level
          FArrayBox& thisNewData = (*newData[lev])[dit];
          
          const BaseFab<int>& validMask = (*mask[lev])[dit];
          const BaseFab<int>& floatingMask = coords[lev]->getFloatingMask()[dit];
          const FArrayBox& thk = coords[lev]->getH()[dit];
          const FArrayBox& topg = coords[lev]->getTopography()[dit];
          const FArrayBox& u = (*vel[lev])[dit];
          const FArrayBox& thisDragCoef = (*dragCoef[lev])[dit];
          const Real& rhoi = coords[lev]->iceDensity();
          const Real& rhow = coords[lev]->waterDensity();
          const Real& gravity = coords[lev]->gravity();
          const Real flotationThicknessFactor = -rhow/rhoi;
          
          //initializing the data in each farraybox (data on each box) to be zero everywhere
          thisNewData.setVal(0.0);
          
          for (int fdir = 0; fdir < SpaceDim; ++fdir)
            {
              //const FArrayBox& fthk =  coords[lev]->getFaceH()[dit][fdir];
              
              Box b = grids[lev][dit];
              for (BoxIterator bit(b);bit.ok();++bit)
                {
                  IntVect iv = bit();

                  // set basalDrag = dragCoef * mag(vel)
                  thisNewData(iv,STRBASEMAG) = thisDragCoef(iv,0)*sqrt((u(iv,0)*u(iv,0)) + u(iv,1)*u(iv,1));
                  
                  //Set the floating fraction value at ivp to be 1.0 if the ice is floating
                  if (validMask(iv) == 1 && floatingMask(iv) == FLOATINGMASKVAL)
                    {
                      thisNewData(iv, SFTFLF) = 1.0;
                    }
		  //grounded ice fraction and grounding line flux
		  else if (validMask(iv) == 1 && floatingMask(iv) == GROUNDEDMASKVAL)
                    {
		      thisNewData(iv, SFTGRF) = 1.0;
                    }  
                      
                  // GL and CF fluxes
		  for (int sign = -1; sign <=1; sign+=2) // loop for (left or right)(sign == -1) or (up or down)(sign == 1) 
		    {
		      const IntVect ivp = iv+sign*BASISV(fdir);
		      // grounding line face at sheet/shelf interface
		      if ((floatingMask(iv) == GROUNDEDMASKVAL) && (floatingMask(ivp) == FLOATINGMASKVAL))
			{
			  Real fvel = 0.5*(u(iv,fdir)+u(ivp,fdir));
			  Real ftopg = min(0.0, 0.5*(topg(iv)+topg(ivp)));
			  thisNewData(iv, LIGROUNDF) += sign* fvel * flotationThicknessFactor * ftopg;
			}
		      Real tol = TINY_THICKNESS;
		      // grounding line face at sheet / no ice interface
		      if ((floatingMask(iv) == GROUNDEDMASKVAL) &&  (floatingMask(ivp) == OPENSEAMASKVAL))
			{
			  Real ftopg = min(0.0, 0.5*(topg(iv)+topg(ivp)));
			  thisNewData(iv, LIGROUNDF) += max(sign*u(iv,fdir),0.0) * flotationThicknessFactor * ftopg;
			}
		      
		      //calving front face
		      if ((thk(iv) >= tol) && (floatingMask(ivp) == OPENSEAMASKVAL))
			{
			  thisNewData(iv, LICALVF) += thk(iv)*max(sign*u(iv,fdir),0.0);
			}
		    } //end loop for (left or right)(sign == -1) or (up or down)(sign == 1) 
                } //end loop over cells
            } //end loop over direction (fdir= 0 (x) or fdir=1 (y))
        } //end loop over grids
    } //end loop over levels

  WriteAMRHierarchyHDF5(out_file, grids, newData, newNames, domainBox, crseDx, dt, time, ratio, numLevels);

  for (int lev = 0; lev < numLevels; lev++)
    {
      if (newData[lev] != NULL)
        {
          delete newData[lev];
          newData[lev] = NULL;
        }	  
     
      if (vel[lev] != NULL)
        {
          delete vel[lev];
          vel[lev] = NULL;
        }
     
      if (mask[lev] != NULL)
        {
          delete mask[lev];
          mask[lev] = NULL;
        }
      
      if (data[lev] != NULL)
        {
          delete data[lev];
	  data[lev] = NULL;
        }
    }
  

  
 }  // end nested scope
 CH_TIMER_REPORT();
 
#ifdef CH_MPI
 MPI_Finalize();
#endif
 
 return 0;
}
