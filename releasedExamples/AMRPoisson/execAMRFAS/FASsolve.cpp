#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRNonLinearPoissonOp.H"
#include "AMRPoissonOp.H"
#include "VCAMRPoissonOp2.H"
#include "AMRFASMultiGrid.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "computeNorm.H"
#include "FABView.H"

#include "UsingNamespace.H"

int s_verbosity = 1;

enum probTypes {exactPrb,
                inexact,
                gaussians,
                numProbTypes};


//int s_probtype = exactPrb;
int s_probtype = gaussians;

//  -----------------------------------------
// boundary condition stuff
//  -----------------------------------------
///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void
setACoef(LevelData<FArrayBox>& a_aCoef)
{
  DataIterator dit = a_aCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& aCoef = a_aCoef[dit];
      ForAllXBNN(Real, aCoef, aCoef.box(), 0, aCoef.nComp());
      {
        // constant-coefficient
        aCoefR = 1.0;
      }EndFor;
    } // end loop over grids
}

void
setBCoef(LevelData<FluxBox>& a_bCoef)
{
  DataIterator dit = a_bCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisBCoef = a_bCoef[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& dirFlux = thisBCoef[dir];
          const Box& dirBox = dirFlux.box();
          ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
            {
              // constant-coefficient
              dirFluxR = -1.0;
            }EndFor
        } // end loop over directions
    }
}

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values,
                Real* b_values,
                Real* c_values)
{
  a_values[0]=0.;
  b_values[0]=0.;
  c_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             RealVect a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                }
            } // end if is not periodic in ith direction
        }
    }
}

void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<RealVect>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setRHS");

  Real gamma = 0;
  ParmParse pp("solver");
  pp.query("gamma", gamma);

  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
      const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

      // rhs is cell-centered...
      RealVect ccOffset = 0.5*a_amrDx[lev];

      DataIterator levelDit = levelGrids.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          FArrayBox& thisRhs = levelRhs[levelDit];

          if (s_probtype == exactPrb)
            {
              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  D_TERM(Real x = loc[0];, Real y = loc[1];, Real z = loc[2];)
                  Real mult_arg = D_TERM((x-x*x),*(y-y*y),*(z-z*z));
                  Real plus_arg = D_TERM((x-x*x),+(y-y*y),+(z-z*z));

                  thisRhs(iv, 0) = ( 2*plus_arg + gamma*mult_arg*exp(mult_arg));
                }
            }
          else if (s_probtype == inexact)
            {
              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  Real x=1;
                  Real y=0;
                  Real z=1;
                  D_TERM(x = loc[0];, y = loc[1];,  z = loc[2];)
                  thisRhs(iv, 0) = ((9*M_PI*M_PI + gamma*exp((x*x-x*x*x)*sin(3*M_PI*(y))))*(x*x-x*x*x) +6*x -2)*sin(3*M_PI*y);

                }
            }  
          else if (s_probtype == gaussians)
            {
              int numGaussians = 3;
              Vector<RealVect> center(numGaussians,RealVect::Zero);
              Vector<Real> scale(numGaussians, 1.0);
              Vector<Real> strength(numGaussians, 1.0);

              for (int n=0; n<numGaussians; n++)
                {
                  if (n==0)
                    {
                      strength[0] = 1.0;
                      scale[0] = 1.0e-2;
                      center[0] = 0.25*RealVect::Unit;
                    }
                  else if (n == 1)
                    {
                      strength[1] = 3.0;
                      scale[1] = 1.0e-2;
                      center[1] = RealVect(D_DECL(0.5,0.75, 0.75));
                    }
                  else if (n == 2)
                    {
                      strength[2] = 2.0;
                      scale[2] = 1.0e-2;
                      center[2] = RealVect(D_DECL(0.75,0.5, 0.5));
                    }
                  else
                    {
                      MayDay::Error("too many Gaussian sources attempted");
                    }
                }

              thisRhs.setVal(0.0);

              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  for (int n=0; n<numGaussians; n++)
                    {
                      RealVect dist = loc - center[n];
                      Real radSqr = D_TERM(dist[0]*dist[0],
                                           +dist[1]*dist[1],
                                            +dist[2]*dist[2]);

                       Real val = strength[n]*exp(-radSqr/scale[n]);
                       thisRhs(iv,0) += val;
                     }
                 }
             }
           else
             {
               MayDay::Error("undefined problem type");
             }
         } // end loop over grids on this level
     } // end loop over levels
 }


void setExact(Vector<LevelData<FArrayBox>* > a_rhs,
              Vector<ProblemDomain>& a_amrDomains,
              Vector<int>& a_refRatios,
              Vector<RealVect>& a_amrDx,
              int a_finestLevel)
{
  CH_TIME("setExact");

  Real gamma = 0;
  ParmParse pp("solver");
  pp.query("gamma", gamma);

  for (int lev=0; lev<=a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

    // rhs is cell-centered...
    RealVect ccOffset = 0.5*a_amrDx[lev];

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisRhs = levelRhs[levelDit];

      if (s_probtype == exactPrb)
      {
        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;


          D_TERM(Real x = loc[0];, Real y = loc[1];, Real z = loc[2];)
          Real arg = D_TERM((x-x*x),*(y-y*y),*(z-z*z));
          thisRhs(iv, 0) = arg;

        }
      }
      else if (s_probtype == inexact)
      {
        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;

          Real x=1;
          Real y=0;
          Real z=1;
          D_TERM(x = loc[0];, y = loc[1];, z = loc[2];)
          thisRhs(iv, 0) = (x*x - x*x*x)*sin(3*M_PI*y);

        }
      }
      else if (s_probtype == gaussians)
      {
            thisRhs.setVal(0.0);
      }
      else
      {
        MayDay::Error("undefined problem type");
      }
    } // end loop over grids on this level
  } // end loop over levels
}



void
setupGrids(Vector<DisjointBoxLayout>& a_amrGrids,
           Vector<ProblemDomain>& a_amrDomains,
           Vector<int>& a_refRatios,
           Vector<RealVect>& a_amrDx,
           int& a_finestLevel)
{
  CH_TIME("setupGrids");

  a_finestLevel = 0;
  ParmParse ppGrids("grids");

  // get grid generation parameters
  int maxLevel, maxBoxSize, blockFactor;
  Real fillRatio;

  ppGrids.get("max_level", maxLevel);

  ppGrids.get("max_box_size",maxBoxSize);

  ppGrids.get("block_factor", blockFactor);

  ppGrids.get("fillRatio", fillRatio);

  // note that there only need to be numLevels-1 refinement ratios
  a_refRatios.resize(maxLevel);
  ppGrids.getarr("ref_ratio", a_refRatios, 0, maxLevel);

  Vector<int>  is_periodic_int;
  bool is_periodic[SpaceDim];
  ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

  IntVect numCells;
  Vector<int> incells(SpaceDim);
  ppGrids.getarr("num_cells", incells, 0, SpaceDim);
  numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                             incells[3],incells[4],incells[5]) );

  RealVect domainSize = RealVect::Unit;
  if (ppGrids.contains("domain_size"))
    {
      Vector<Real> insize(SpaceDim);
      ppGrids.getarr("domain_size", insize, 0, SpaceDim);
      domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                              insize[3],insize[4],insize[5]) );
    }

  // resize dataholders
  int maxNumLevels = maxLevel +1;
  a_amrGrids.resize(maxNumLevels);
  a_amrDomains.resize(maxNumLevels);
  a_amrDx.resize(maxNumLevels);
  a_finestLevel = 0;

  // assumes anisotropic mesh possible
  a_amrDx[0][0] = domainSize[0]/numCells[0];
  a_amrDx[0][1] = domainSize[1]/numCells[1];
  if (CH_SPACEDIM == 3) {
      a_amrDx[0][2] = domainSize[2]/numCells[2];
  }

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = numCells - IntVect::Unit;

  ProblemDomain baseDomain(domLo, domHi, is_periodic);
  a_amrDomains[0] = baseDomain;

  // set up refined domains, etc
  for (int lev=1; lev<= maxLevel; lev++)
    {
      a_amrDomains[lev] = a_amrDomains[lev-1];
      a_amrDomains[lev].refine(a_refRatios[lev-1]);
      a_amrDx[lev][0] = a_amrDx[lev-1][0]/a_refRatios[lev-1];
      a_amrDx[lev][1] = a_amrDx[lev-1][1]/a_refRatios[lev-1];
      if (CH_SPACEDIM == 3) {
          a_amrDx[lev][2] = a_amrDx[lev-1][2]/a_refRatios[lev-1];
      }
    }

  Vector<Vector<Box> > vectBoxes(maxLevel+1);

  // local scope. for base-level grid generation
  {
    CH_TIME("BaseGridCreation");
    // generate base level grids

    domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

    Vector<int> procAssign(vectBoxes[0].size(), 0);

    LoadBalance(procAssign, vectBoxes[0]);

    DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

    a_amrGrids[0] = baseGrids;
  }


  if (maxLevel > 0)
    {
      bool read_grids = false;
      ppGrids.query("read_in_grids", read_grids);
      if (read_grids)
        {
          for (int ilev = 1; ilev <= maxLevel; ilev++)
            {
              const ProblemDomain& levDomain = a_amrDomains[ilev];

              Vector<Box>   boxes;
              char boxCountVar[100];
              int boxCount;
              sprintf(boxCountVar, "level_%d_box_count", ilev);
              ppGrids.get(boxCountVar, boxCount);
              boxes.resize(boxCount);
              for (int ibox = 0; ibox < boxCount; ibox++)
                {
                  char boxLoVar[100];
                  char boxHiVar[100];
                  sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
                  sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
                  Vector<int> boxLo, boxHi;
                  ppGrids.getarr(boxLoVar, boxLo, 0, SpaceDim);
                  ppGrids.getarr(boxHiVar, boxHi, 0, SpaceDim);
                  IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
                  IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
                  boxes[ibox] = Box(ivLo, ivHi);
                  if (!levDomain.contains(boxes[ibox]))
                    {
                      MayDay::Error("box outside of domain");
                    }
                }
              //check to see if level 0 domain is covered
              if (ilev == 0)
                {
                  IntVectSet ivDom(levDomain.domainBox());
                  for (int ibox = 0; ibox < boxes.size(); ibox++)
                    {
                      ivDom -= boxes[ibox];
                    }
                  if (!ivDom.isEmpty())
                    {
                      MayDay::Error("level 0 boxes must cover the domain");
                    }
                }
              Vector<int>  proc(boxes.size());
              LoadBalance(proc,boxes);
              a_amrGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
              a_finestLevel++;
            }

        }
      else
        {
          // tag on grad(rhs)
          int bufferSize = 1;
          //ppGrids.query("buffer_size", bufferSize);
          BRMeshRefine meshGen(a_amrDomains[0],
                               a_refRatios,
                               fillRatio,
                               blockFactor,
                               bufferSize,
                               maxBoxSize);

          // to be used by MeshRefine...
          Vector<Vector<Box> > oldMeshes(maxLevel+1);
          oldMeshes[0] = vectBoxes[0];
          for (int lev=1; lev<oldMeshes.size(); lev++)
            {
              oldMeshes[lev].push_back(a_amrDomains[lev].domainBox());
            }

          Real refineThresh;
          ppGrids.get("refine_threshold", refineThresh);

          Real threshSqr = refineThresh*refineThresh;

          bool moreLevels = true;
          while (moreLevels)
            {
              // tag based on grad(rhs)
              // first need to allocate RHS
              Vector<LevelData<FArrayBox>* > tempRHS(a_finestLevel+1, NULL);
              for (int lev=0; lev<= a_finestLevel; lev++)
                {
                  // note that we add a ghost cell to simplify gradients
                  tempRHS[lev] = new LevelData<FArrayBox>(a_amrGrids[lev],
                                                          1, IntVect::Unit);
                }

              setRHS(tempRHS, a_amrDomains, a_refRatios, a_amrDx, 
                     a_finestLevel);

              Vector<IntVectSet> tags(a_finestLevel+1);

              for (int lev=0; lev<a_finestLevel+1; lev++)
                {
                  const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
                  const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
                  IntVectSet& levelTags = tags[lev];

                  // compute mag(gradient)
                  DataIterator dit = levelGrids.dataIterator();
                  for (dit.begin(); dit.ok(); ++dit)
                    {
                      const FArrayBox& rhsFab = levelRHS[dit];
                      // local storage for gradient
                      FArrayBox gradFab(levelGrids[dit],1);
                      gradFab.setVal(0.0);
                      Real thisGrad;

                      BoxIterator bit(levelGrids[dit]);
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv=bit();
                          for (int dir=0; dir<SpaceDim; dir++)
                            {
                              // use mag(undivided gradient)
                              IntVect hi = iv + BASISV(dir);
                              IntVect lo = iv - BASISV(dir);
                              thisGrad = rhsFab(hi,0) - rhsFab(lo,0);
                              gradFab(iv,0) += (thisGrad*thisGrad);
                            } // end loop over directions
                        } // end loop over cells

                      //gradFab now has mag(grad*dx)^2

                      // tag where mag(gradient) > tolerance^2
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv = bit();
                          if (gradFab(iv,0) > threshSqr)
                            {
                              levelTags |= iv;
                            }
                        } // end loop over cells
                    } // end loop over grids on this level

                } // end loop over levels


              // call meshRefine.
              for (int lev=1; lev<=a_finestLevel; lev++)
                {
                  oldMeshes[lev] = vectBoxes[lev];
                }

              int topLevel = a_finestLevel;
              int newFinestLevel =  meshGen.regrid(vectBoxes,
                                                   tags,
                                                   0,
                                                   topLevel,
                                                   oldMeshes);


              // define new grids if necessary and test to see if we're done
              if (newFinestLevel > a_finestLevel)
                {
                  a_finestLevel = newFinestLevel;

                  // setup new grid hierarchy
                  for (int lev=1; lev<=a_finestLevel; lev++)
                    {
                      Vector<int> procAssign(vectBoxes[lev].size(),0);
                      LoadBalance(procAssign, vectBoxes[lev]);
                      DisjointBoxLayout levelGrids(vectBoxes[lev],
                                                   procAssign,
                                                   a_amrDomains[lev]);
                      a_amrGrids[lev] = levelGrids;
                    }
                  if (s_verbosity>2) pout() << "setupGrids: "<< a_finestLevel <<") size " << a_amrGrids[a_finestLevel].size() << endl;
                }
              else
                {
                  moreLevels = false;
                }

              if (a_finestLevel == maxLevel)
                {
                  moreLevels = false;
                }

              // clean up before starting again
              for (int lev=0; lev<tempRHS.size(); lev++)
                {
                  delete tempRHS[lev];
                }

            } // end while (moreLevels)

        }

      // fill in remaining levels with empty DisjointBoxLayouts
      for (int lev= a_finestLevel+1; lev<=maxLevel; lev++)
        {
          a_amrGrids[lev] = DisjointBoxLayout();
        }

    }


}


void
setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<RealVect>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setupSolver");

  ParmParse ppSolver("solver");

  int poissonOp = 1;
  ppSolver.query("poissonOp", poissonOp);

  bool FASmultigrid = true;
  ppSolver.query("FASmultigrid", FASmultigrid);

  int numLevels = a_finestLevel+1;

  if (poissonOp == 1) {
      AMRNonLinearPoissonOpFactory opFactory;

      // solving nonlinear poisson problem here
      Real alpha = 0.0;
      Real beta =  1.0;
      Real gamma = 0.0;

      ppSolver.query("gamma", gamma);

      opFactory.define(a_amrDomains[0],
                       a_amrGrids,
                       a_refRatios,
                       a_amrDx[0],
                       &ParseBC, 
                       alpha, beta, gamma, FASmultigrid);


      AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

      a_amrSolver->define(a_amrDomains[0], castFact,
                          &a_bottomSolver, numLevels);

  } else if (poissonOp == 2) {
      AMRPoissonOpFactory opFactory;

      // solving poisson problem here
      Real alpha = 0.0;
      Real beta = -1.0;
   
      opFactory.define(a_amrDomains[0],
                       a_amrGrids,
                       a_refRatios,
                       a_amrDx[0],
                       &ParseBC, alpha, beta, FASmultigrid);

      AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

      a_amrSolver->define(a_amrDomains[0], castFact,
                          &a_bottomSolver, numLevels);
  } else if (poissonOp == 3) {
      VCAMRPoissonOp2Factory opFactory;

      // solving nonlinear poisson problem here
      Real alpha = 0.0;
      Real beta  = -1.0;

      int maxNumLevels = a_amrGrids.size();

      Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef;
      Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef;
      aCoef.resize(maxNumLevels);
      bCoef.resize(maxNumLevels);
      for (int lev=0; lev<numLevels; lev++) {
         aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_amrGrids[lev], 1, IntVect::Zero));
         bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(a_amrGrids[lev], 1, IntVect::Zero));

         setACoef(*aCoef[lev]);
         setBCoef(*bCoef[lev]);
      }
      for (int lev=numLevels; lev<maxNumLevels; lev++) {
         aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>);
         bCoef[lev] = RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>);
      }
      opFactory.define(a_amrDomains[0],
                       a_amrGrids,
                       a_refRatios,
                       a_amrDx[0],
                       &ParseBC, alpha, aCoef,
                       beta, bCoef, FASmultigrid);

      AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

      a_amrSolver->define(a_amrDomains[0], castFact,
                          &a_bottomSolver, numLevels);
  } else {
      MayDay::Error("wrong choice for poisson Operator");
  }

  // multigrid solver parameters
  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ppSolver.get("num_smooth", numSmooth);
  ppSolver.get("num_mg",     numMG);
  ppSolver.get("max_iterations", maxIter);
  ppSolver.get("tolerance", eps);
  ppSolver.get("hang",      hang);

  Real normThresh = 1.0e-30;
 
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh, !FASmultigrid);// last param is homogeneous BC
  a_amrSolver->m_verbosity = s_verbosity-1;

  // optional parameters
  ppSolver.query("num_pre", a_amrSolver->m_pre);
  ppSolver.query("num_post", a_amrSolver->m_post);
  ppSolver.query("num_bottom", a_amrSolver->m_bottom);
}


int runSolver()
 {
   CH_TIME("runSolver");

   int status = 0;
   ParmParse ppMain("main");
   ParmParse ppSolver("solver");

   ppMain.query("verbosity", s_verbosity);
   ppMain.query("problem", s_probtype);

   // set up grids&
   Vector<DisjointBoxLayout> amrGrids;
   Vector<ProblemDomain> amrDomains;
   Vector<int> refRatios;
   Vector<IntVect> refRatios_anys;
   Vector<RealVect> amrDx;
   int finestLevel;

   setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);

   pout() << "max level vs finest level " << refRatios.size() << " " << finestLevel << endl;

   refRatios_anys.resize(refRatios.size());
   for (int lev=0; lev<refRatios.size(); lev++) {
       refRatios_anys[lev] = refRatios[lev] * IntVect::Unit;
   }

  // initialize solver
  bool FASmultigrid = true;
  ppSolver.query("FASmultigrid", FASmultigrid);

  // allocate solution and RHS, initialize RHS
  int numLevels = amrGrids.size();
  Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
  // this is for convenience
  Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > exact(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > error(numLevels, NULL);

   for (int lev=0; lev<=finestLevel; lev++)
     {
       const DisjointBoxLayout& levelGrids = amrGrids[lev];
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       exact[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       error[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
     }

   setExact(exact, amrDomains, refRatios, amrDx, finestLevel );
   setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );

   // do solve
   int iterations = 3;
   ppMain.get("iterations", iterations);

   for (int iiter = 0; iiter < iterations; iiter++)
     {
         // initialize solver
         AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
         if (FASmultigrid) {
             amrSolver = new AMRFASMultiGrid<LevelData<FArrayBox> >();
         } else {
             amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
         }
         BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
         bottomSolver.m_verbosity = s_verbosity-2;
         setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
                     refRatios, amrDx, finestLevel);

         bool zeroInitialGuess = true;
         // if prob type is exact
         //ppSolver.query("zeroInitialGuess", zeroInitialGuess);
         //setExact(phi, amrDomains, refRatios, amrDx, finestLevel );
         amrSolver->solve(phi, rhs, finestLevel, 0, zeroInitialGuess);

         if (iiter == iterations-1) {
             pout() << "End iterations. norm=" << amrSolver->computeAMRResidual(resid,phi,rhs,finestLevel,0) << endl;
         }

         delete amrSolver;
      }

   // Compute error -- if prob type is exact
   if ( (s_probtype == exactPrb) || (s_probtype == inexact) ) {
       for (int lev=0; lev<=finestLevel; lev++)
       {
         phi[lev]->copyTo(*error[lev]);
         for (DataIterator dit = error[lev]->dataIterator(); dit.ok(); ++dit)
         {
           (*error[lev])[dit].minus((*exact[lev])[dit]);
         }
       }

       // Compute error metrics
       Real max, L1, L2;
       max = computeNorm(error,   refRatios_anys, amrDx[0], Interval(0,0), 0);
       L1 = computeNorm(error,    refRatios_anys, amrDx[0], Interval(0,0), 1);
       L2 = computeNorm(error,    refRatios_anys, amrDx[0], Interval(0,0), 2);

       char errStr[1024];
       sprintf(errStr, "Error = %1.2e (max), %1.2e (L1), %1.2e (L2)", max, L1, L2 );
       pout() << errStr << endl;
   }

   // write results to file

   bool writePlots = true;
   ppMain.query("writePlotFiles", writePlots);

   int numPlotComps = 5;
#ifdef CH_USE_HDF5

   if (writePlots)
     {
       int numLevels = finestLevel +1;
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);
       
       
       for (int lev=0; lev<numLevels; lev++)
         {
           plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                                    numPlotComps, IntVect::Zero);

           Interval phiInterval(0,0);
           phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
           Interval rhsInterval(1,1);
           rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
           Interval resInterval(2,2);
           resid[lev]->copyTo(phiInterval, *plotData[lev], resInterval);
           Interval exactInterval(3,3);
           exact[lev]->copyTo(phiInterval, *plotData[lev], exactInterval);
           Interval errInterval(4,4);
           error[lev]->copyTo(phiInterval, *plotData[lev], errInterval);
         }

       string fname = "FASOut.";

       char suffix[30];
       sprintf(suffix, "%dd.hdf5",SpaceDim);
       fname += suffix;

       Vector<string> varNames(numPlotComps);
       varNames[0] = "phi";
       varNames[1] = "rhs";
       varNames[2] = "res";
       varNames[3] = "exact";
       varNames[4] = "error";

       Real bogusVal = 1.0;

       //WriteAnisotropicAMRHierarchyHDF5(fname,
       WriteAMRHierarchyHDF5(fname,
                             amrGrids,
                             plotData,
                             varNames,
                             amrDomains[0].domainBox(),
                             //amrDx[0],
                             amrDx[0][0],
                             bogusVal,
                             bogusVal,
                             //refRatios_anys,
                             refRatios,
                             numLevels);

       // clean up
       for (int lev=0; lev<plotData.size(); lev++)
         {
           delete plotData[lev];
         }
     } // end if writing plots
#endif // end if HDF5

   // clean up
   for (int lev=0; lev<phi.size(); lev++)
     {
       delete phi[lev];
       delete rhs[lev];
       delete resid[lev];
       delete exact[lev];
       delete error[lev];
     }


   return status;
}


/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runSolver();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
