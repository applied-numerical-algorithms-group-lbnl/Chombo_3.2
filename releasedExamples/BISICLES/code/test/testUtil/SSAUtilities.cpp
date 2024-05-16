#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "SSAUtilities.H"
#include "functionsF_F.H"
#include "PoissProbF_F.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "CoarseAverage.H"

#include "NamespaceHeader.H"

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}


void neumannBC(FArrayBox& a_state,
               const Box& a_valid,
               const ProblemDomain& a_domain,
               Real a_dx,
               bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  //Real bcVal = 0.0;
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}


void dirichletBC(FArrayBox& a_state,
                 const Box& a_valid,
                 const ProblemDomain& a_domain,
                 Real a_dx,
                 bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}


extern void ParseBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      if(!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }
      
      Box grownValid = a_valid;
      // do this so that we can fill corner clels
      grownValid.grow(1);
      for(int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              // do this so that we can get corner cells filled
              // effect of grownValid + this is a box which is the 
              // valid box grown by one in all transverse directions
              Box valid = grownValid;
              valid.grow(i,-1);
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if(GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const Neumann bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);
                    }
                  else if(GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const Dirichlet bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);

                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if(GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const Neumann bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if(GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if(!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const Dirichlet bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if is not periodic in ith direction
        }
    }
  

  // hardwire to Neumann BC's for the moment
  //neumannBC(a_state, a_valid, a_domain, a_dx, a_homogeneous);
}


void
outputData(const Vector<RefCountedPtr<LevelData<FArrayBox> > >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain&     domain,
           const Vector<int>& vectRatio,
           RealVect dxCoarsest,
           int numlevels,
           string filename,
           string varname)
{
#ifdef CH_USE_HDF5
  Vector<string> vectName(vectPhi[0]->nComp(), varname);
  for (int n=0; n<vectName.size(); n++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "%s_%d", varname.c_str(), n);
      string label(labelChSt);
      vectName[n] = label;
    }
  Real time = 1;  //placeholder

  // need to remove RefCountedPtr wrapper...
  Vector<LevelData<FArrayBox>* > tempPhi(vectPhi.size(), NULL);
  for (int lev=0; lev<vectPhi.size(); lev++)
    {
      tempPhi[lev] = const_cast<LevelData<FArrayBox>*>(vectPhi[lev].operator->());
    }
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        tempPhi,
                        vectName,
                        domain.domainBox(),
                        dxCoarsest[0], time, time,
                        vectRatio,
                        numlevels);
#endif
}


void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain&     domain,
           const Vector<int>& vectRatio,
           RealVect dxCoarsest,
           int numlevels,
           string filename,
           string varname)
{
#ifdef CH_USE_HDF5
  Vector<string> vectName(vectPhi[0]->nComp(), varname);
  for (int n=0; n<vectName.size(); n++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "%s_%d", varname.c_str(), n);
      string label(labelChSt);
      vectName[n] = label;
    }

  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        vectName,
                        domain.domainBox(),
                        dxCoarsest[0], time, time,
                        vectRatio,
                        numlevels);
#endif
}

int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           PoissonParameters&              a_params,
           int a_phiType,
           int a_numLevels)

{
  int numlevels = a_numLevels;
  if(numlevels < 0)
    {
      numlevels = a_params.numLevels;
    }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real rhono, rno;
  int iprob;
  ParmParse pp;
  pp.get("iprob", iprob);

  rhono = 0.75;
  rno = 0.5;

  RealVect dxLev = a_params.coarsestDx;
  ProblemDomain domlev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      if(iprob == -1)
        {
          // prob for convergence tests
          setLOfPhi(*vectRhs[ilev], dxLev, a_params, a_phiType);
        }
      else
        {
          LevelData<FArrayBox>& rhsLD = *vectRhs[ilev];
          DataIterator dit =  rhsLD.dataIterator();
          for(dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& rhsFab = rhsLD[dit()];
              rhsFab.setVal(0.);
              bool useEBGrids;
              pp.get("use_eb_grids", useEBGrids);

              if(useEBGrids)
                {
                  rhsFab.setVal(1.0);
                }
              else
                {
                  FORT_GETRHSPOIS(CHF_FRA(rhsFab),
                                  CHF_BOX(rhsFab.box()),
                                  CHF_BOX(domlev.domainBox()),
                                  CHF_CONST_REALVECT(dxLev),
                                  CHF_CONST_REAL(rhono),
                                  CHF_CONST_REAL(rno),
                                  CHF_CONST_INT(iprob));
                }
            }
        }
      dxLev /= a_params.refRatio[ilev];
      domlev.refine(a_params.refRatio[ilev]);

    }
  return 0;
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<RealVect>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels)
{
  for (int lev=baseLevel; lev!= numLevels; lev++)
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;

      maxRHS = norm(levelRhs, levelRhs.interval(), 0);

      Real tagVal = maxRHS * refine_thresh;

      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level

      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;

      tagVect[lev] = local_tags;

    } // end loop over levels
}

/*
  Set grid hierarchy from input file
 */
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<RealVect>&              vectDx,
                         PoissonParameters&         a_params)
{

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(    a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for(int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      vectDx[ilev] = vectDx[ilev-1]/a_params.refRatio[ilev-1];
    }


  vectDomain[0] = a_params.coarsestDomain;
  for(int ilev = 1;ilev < a_params.numLevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],a_params.refRatio[ilev-1]);
    }
}

int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             PoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<RealVect>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp;
  bool useEBGrids;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  pp.query("use_eb_grids", useEBGrids);
  if(pp.contains("read_in_grids"))
    {
      pp.get("read_in_grids", readInGrids);
    }

  if(readInGrids)
    {

      ProblemDomain levDomain = a_params.coarsestDomain;
      for(int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          Vector<Box>   boxes;
          char boxCountVar[100];
          int boxCount;
          sprintf(boxCountVar, "level_%d_box_count", ilev);
          pp.get(boxCountVar, boxCount);
          boxes.resize(boxCount);
          for(int ibox = 0; ibox < boxCount; ibox++)
            {
              char boxLoVar[100];
              char boxHiVar[100];
              sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
              sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
              Vector<int> boxLo, boxHi;
              pp.getarr(boxLoVar, boxLo, 0, SpaceDim);
              pp.getarr(boxHiVar, boxHi, 0, SpaceDim);
              IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
              IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
              boxes[ibox] = Box(ivLo, ivHi);
              if(!levDomain.contains(boxes[ibox]))
                {
                  MayDay::Error("box outside of domain");
                }
            }
          //check to see if level 0 domain is covered
          if(ilev == 0)
            {
              IntVectSet ivDom(levDomain.domainBox());
              for(int ibox = 0; ibox < boxes.size(); ibox++)
                {
                  ivDom -= boxes[ibox];
                }
              if(!ivDom.isEmpty())
                {
                  MayDay::Error("level 0 boxes must cover the domain");
                }
            }
          Vector<int>  proc(a_params.numLevels);
          LoadBalance(proc,boxes);
          vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
          levDomain.refine(a_params.refRatio[ilev]);
        }

    }
  else if(useEBGrids)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple  grids" << endl;
      Vector<int> proc(1, 0);
      Box coarBox = vectDomain[0].domainBox();
      Vector<Box> coarBoxes(1, coarBox);
      vectGrids[0] = DisjointBoxLayout(coarBoxes, proc, vectDomain[0]);

      for(int ilev = 1; ilev < numlevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if(iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(a_params.refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          vectGrids[ilev] = DisjointBoxLayout(refBoxes, proc,
                                              vectDomain[ilev]);
        }
    }
  else
    {
      pout() << "tagging on gradient of RHS" << endl;
      int maxLevel = numlevels-1;
      Vector<Vector<Box> > newBoxes(numlevels);
      Vector<Vector<Box> > oldBoxes(numlevels);

      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);
      int ncomps = 1;

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize, a_params.blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);

      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);

      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps,
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);

      int nesting_radius = 2;
      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio,
                              a_params.fillRatio,
                              a_params.blockFactor, nesting_radius,
                              a_params.maxGridSize);

      while (moreLevels)
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy
          setRHS(vectRHS, a_params, topLevel+1);

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain,
                   a_params.refineThresh,
                   tags_grow, baseLevel, topLevel+1);

          int new_finest = meshrefine.regrid(newBoxes, tagVect,
                                             baseLevel,
                                             topLevel, oldBoxes);

          if (new_finest > topLevel)
            {
              topLevel++;
            }

          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++)
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              delete vectRHS[lev];
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomps,
                                                      IntVect::Zero);
            } // end loop over levels for initialization

          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop
      // clean up temp storage
      for (int ilev=0; ilev <vectRHS.size(); ilev++)
        {
          if (vectRHS[ilev] != NULL)
            {
              delete vectRHS[ilev];
              vectRHS[ilev] = NULL;
            }
        }
    }
  return 0;
}


void
defineMGSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
               const Vector<DisjointBoxLayout>&             a_grids,
               Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
               Vector<RefCountedPtr<LevelData<SigmaCS> > >&  a_coordSys,
               LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
               const PoissonParameters&                     a_params)
{
  ParmParse pp2;

  // first try homogeneous Dirichlet BC's
  IntVect loBCType = IntVect::Unit;
  IntVect hiBCType = IntVect::Unit;
  RealVect loBCVal = RealVect::Zero;
  RealVect hiBCVal = RealVect::Zero;
  
  RefCountedPtr<BCFunction> bcPtr = ConstDiriNeumBC(loBCType, loBCVal,
                                                    hiBCType, hiBCVal);
  
  BCHolder bc(bcPtr);
  
  SSAVelocityOpFactory opFactory(a_grids,
                                 a_mu,
                                 a_coordSys,
                                 a_params.refRatio,
                                 a_params.coarsestDomain,
                                 a_params.coarsestDx,
                                 bc);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory.m_coefficient_average_type = a_params.coefficient_average_type;
    }

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = a_params.verbosity;

}

// returns pointer to multilevelLinearOp so that we don't leak memory.
// this is obviously not a great thing, but it's only a test code.
MultilevelLinearOp<FArrayBox>* 
defineBiCGStabSolver(BiCGStabSolver<Vector<LevelData<FArrayBox>* > >&  a_solver,
                     const Vector<DisjointBoxLayout>&             a_grids,
                     Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
                     Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
                     AMRMultiGrid<LevelData<FArrayBox> >& a_precondSolver,
                     const PoissonParameters&                     a_params)
{
  ParmParse pp2;
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory
    = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
    (defineOperatorFactory(a_grids, a_mu, a_coordSys, a_params));
  
  int lBase = 0;
  MultilevelLinearOp<FArrayBox>* mlOpPtr = new MultilevelLinearOp<FArrayBox>;
  MultilevelLinearOp<FArrayBox>& mlOp = *mlOpPtr;
  int numMGIter = 1;
  pp2.query("num_mg", numMGIter);

  mlOp.m_num_mg_iterations = numMGIter;
  int numMGSmooth = 4;
  pp2.query("num_smooth", numMGSmooth);
  mlOp.m_num_mg_smooth = numMGSmooth;
  int preCondSolverDepth = -1;
  pp2.query("preCondSolverDepth", preCondSolverDepth);
  mlOp.m_preCondSolverDepth = preCondSolverDepth;
  
  Real tolerance = 1.0e-7;
  pp2.query("tolerance", tolerance);
  
  int max_iter = 10;
  pp2.query("max_iterations", max_iter);
  
  Vector<RealVect> vectDx(a_grids.size());
  Vector<ProblemDomain> vectDomain(a_grids.size());
  vectDx[0] = a_params.coarsestDx;
  vectDomain[0] = a_params.coarsestDomain;

  for (int lev=1; lev<vectDx.size(); lev++)
    {
      vectDx[lev] = vectDx[lev-1]/a_params.refRatio[lev-1];
      vectDomain[lev] = vectDomain[lev-1];
      vectDomain[lev].refine(a_params.refRatio[lev-1]);
    }


  
  mlOp.define(a_grids, a_params.refRatio, vectDomain,
              vectDx, opFactory, lBase);
  bool homogeneousBC = false;
  a_solver.define(&mlOp, homogeneousBC);
  a_solver.m_verbosity = a_params.verbosity;
  a_solver.m_normType = 0;
  a_solver.m_eps = tolerance;
  a_solver.m_imax = max_iter;
  
  return mlOpPtr;
}


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_mu,
                      Vector<RefCountedPtr<LevelData<SigmaCS> > >& a_coordSys,
                      const PoissonParameters&                     a_params)
{

  ParmParse pp2;
  // first try homogeneous Dirichlet BC's
  IntVect loBCType = IntVect::Unit;
  IntVect hiBCType = IntVect::Unit;
  RealVect loBCVal = RealVect::Zero;
  RealVect hiBCVal = RealVect::Zero;
  
  RefCountedPtr<BCFunction> bcPtr = ConstDiriNeumBC(loBCType, loBCVal,
                                                    hiBCType, hiBCVal);
  
  BCHolder bc(bcPtr);
  SSAVelocityOpFactory* opFactory = new SSAVelocityOpFactory(a_grids,
                                                             a_mu,
                                                             a_coordSys,
                                                             a_params.refRatio,
                                                             a_params.coarsestDomain,
                                                             a_params.coarsestDx,
                                                             bc);;

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;

}


#if 0
extern void
defineSimpleOp(NewPoissonOp& a_simpleOp,
               const PoissonParameters& a_params)
{
  RealVect dx = a_params.coarsestDx;
  ProblemDomain dom(a_params.coarsestDomain);

  a_simpleOp.define(dx, dom, &ParseBC);
}



#endif




/********/
void setPhi(LevelData<FArrayBox>&    a_phi,
            const RealVect&          a_dx,
            const PoissonParameters& a_params,
            int a_phitype)
{
  CH_assert(a_phi.nComp() == 2);

  //int comp = 0;

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curPhiFAB = a_phi[dit()];
      Box curPhiBox = curPhiFAB.box();
      
      if (a_phitype == constZero)
        {
          curPhiFAB.setVal(0.0);
        }
      else if (a_phitype == constOne)
        {
          curPhiFAB.setVal(1.0);
        }
      else
        {
          MayDay::Error("phiType not implemented yet in getPhi");
#if 0
          const RealVect&     trig = getRV();

          FORT_GETPHI(CHF_FRA1(curPhiFAB,comp),
                      CHF_CONST_REALVECT(trig),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_params.probLo),
                      CHF_CONST_REALVECT(a_params.probHi),
                      CHF_BOX(curPhiBox));
#endif
        } // end switching for phiType
    } // end loop over grids
      
}

/********/
void setLOfPhi(LevelData<FArrayBox>&    a_LOfPhi,
               const RealVect&          a_dx,
               const PoissonParameters& a_params,
               int a_phiType)
{
  CH_assert(a_LOfPhi.nComp() == 2);

  //int comp = 0;
  //const RealVect&  trig = getRV();

  for (DataIterator dit = a_LOfPhi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curLOfPhiFAB = a_LOfPhi[dit()];

      if (a_phiType == constZero)
        {
          curLOfPhiFAB.setVal(0.0);
        }
      else if (a_phiType == constOne)
        {
          curLOfPhiFAB.setVal(0.0);
        }
      else
        {
          MayDay::Error("phiType not implemented in setLOfPhi");
#if 0
          Box curLOfPhiBox = curLOfPhiFAB.box();
          FORT_GETLOFPHI(CHF_FRA1(curLOfPhiFAB,comp),
                         CHF_CONST_REALVECT(trig),
                         CHF_CONST_REALVECT(a_dx),
                         CHF_CONST_REALVECT(a_params.probLo),
                         CHF_CONST_REALVECT(a_params.probHi),
                         CHF_CONST_REAL(a_params.alpha),
                         CHF_CONST_REAL(a_params.beta),
                         CHF_BOX(curLOfPhiBox));
#endif
        } // end switch on phiType

    }


}
/******/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/******/
void compareError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                  const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const PoissonParameters&                 a_paramsFine,
                  const PoissonParameters&                 a_paramsCoar,
                  const string& a_testName)
{
  const Vector<int> refRat = a_paramsCoar.refRatio;
  int lbase = a_paramsFine.baseLevel;
  const int ncomp = a_errorFine[lbase]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  Real* orders    = new Real[ncomp*nnorm];
  for(int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for(int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[iscript(icomp, inorm, ncomp)] = 0;
          normsFine[iscript(icomp, inorm, ncomp)] = 0;
        }
    }
  ParmParse pp;
  pout() << "==============================================" << endl;
  for(int comp = 0; comp < ncomp; comp++)
    {
      pout() << "Comparing error in variable  " << comp << endl;
      pout() << "==============================================" << endl;
      for(int inorm = 0; inorm <= 2; inorm++)
        {

          if(inorm == 0)
            {
                pout() << endl << "Using max norm." << endl;
            }
          else
            {
              pout() << endl << "Using L-" << inorm << " norm." << endl;
            }
          Real dxCoar =refRat[0];
          Real dxFine = 1.0;
          Interval comps(comp,comp);
          Real coarnorm = computeNorm(a_errorCoar, refRat, dxCoar, comps, inorm, lbase);

          Real finenorm = computeNorm(a_errorFine, refRat, dxFine, comps, inorm, lbase);

          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;

          normsCoar[iscript(comp,inorm,ncomp)] = coarnorm;
          normsFine[iscript(comp,inorm,ncomp)] = finenorm;

          if((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
            {
              Real order = log(Abs(coarnorm/finenorm))/log(2.0);
              //pout() << "Order of scheme = " << order << endl;
              orders[iscript(comp,inorm,ncomp)] = order;
            }
        }
      pout() << "==============================================" << endl ;;
    }


  //output in latex format to be safe
  int nfine = a_paramsFine.coarsestDomain.size(0);
  pout() << setw(12)
         << setprecision(6)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific) ;

  for (int inorm = 0; inorm <= 2; inorm++)
    {
      pout() << "\\begin{table}[p]" << endl;
      pout() << "\\begin{center}" << endl;
      pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
      pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
      pout() << "\\hline \\hline " << endl;
      for(int icomp = 0; icomp < ncomp; icomp++)
        {
          int iindex = iscript(icomp,inorm,ncomp);
          pout() << "var" << icomp << " &    \t "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsCoar[iindex]  << " & "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsFine[iindex] << " & "
                 << setw(12)
                 << setprecision(2)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << orders[iindex];
          pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
        }
      pout() << "\\end{tabular}" << endl;
      pout() << "\\end{center}" << endl;
      pout() << "\\caption{";
      pout() << a_testName ;
      pout() << " convergence rates using L-" << inorm << " norm. " << endl;
      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
      pout() << "\\end{table}" << endl;
      pout() << endl << endl;
    }
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;
}

/******/
void compareError(const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_errorFine,
                  const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const PoissonParameters&                 a_paramsFine,
                  const PoissonParameters&                 a_paramsCoar,
                  const string& a_testName)
{
  // strip out pointers from RefCountedPtrs for errors
  Vector<LevelData<FArrayBox>* > fineTemp(a_errorFine.size(), NULL);
  Vector<LevelData<FArrayBox>* > crseTemp(a_errorCoar.size(), NULL);

  for (int lev=a_paramsFine.baseLevel; lev<fineTemp.size(); lev++)
    {
      fineTemp[lev] = const_cast<LevelData<FArrayBox>*>(a_errorFine[lev].operator->());
      crseTemp[lev] = const_cast<LevelData<FArrayBox>*>(a_errorCoar[lev].operator->());
    }

  compareError(fineTemp, crseTemp, a_gridsFine, a_gridsCoar,
               a_paramsFine, a_paramsCoar, a_testName);
}


/********/
void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                              const Vector<DisjointBoxLayout>& a_gridsFine,
                              const PoissonParameters&         a_paramsCoar)
{
  int nlevels = a_paramsCoar.numLevels;
  a_gridsCoar.resize(nlevels);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      CH_assert(a_gridsFine[ilev].coarsenable(2));
      coarsen(a_gridsCoar[ilev], a_gridsFine[ilev], 2);
    }

}
/********/
void PoissonParameters::coarsen(int a_factor)
{
  coarsestDx *= a_factor;
  coarsestDomain.coarsen(a_factor);
}
/********/
void PoissonParameters::refine(int a_factor)
{
  coarsestDx /= a_factor;
  coarsestDomain.refine(a_factor);
}
/********/
void getPoissonParameters(PoissonParameters&  a_params)
{
  ParmParse pp;

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }

  Vector<int> is_periodic(SpaceDim, false);
  pp.queryarr("periodic", is_periodic, 0, SpaceDim);

  pp.get("refine_threshold",a_params.refineThresh);
  pp.get("block_factor",a_params.blockFactor);
  pp.get("fill_ratio",a_params.fillRatio);
  pp.get("buffer_size",a_params.bufferSize);

  // set to a bogus default value, so we only break from solver
  // default if it's set to something real
  a_params.coefficient_average_type = -1;
  if (pp.contains("coefficient_average_type"))
    {
      std::string tempString;
      pp.get("coefficient_average_type", tempString);
      if (tempString == "arithmetic")
        {
          a_params.coefficient_average_type = CoarseAverage::arithmetic;
        }
      else if (tempString == "harmonic")
        {
          a_params.coefficient_average_type = CoarseAverage::harmonic;
        }
      else
        {
          MayDay::Error("bad coefficient_average_type in input");
        }
    } // end if an average_type is present in inputs

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
  a_params.verbosity = 3;
  pp.query("verbosity", a_params.verbosity);
  a_params.baseLevel = 0;
  pp.query("baseLevel", a_params.baseLevel);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  Box crseDomBox(lo,hi);
  ProblemDomain crseDom(crseDomBox);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      crseDom.setPeriodic(dir, is_periodic[dir]);
    }
  a_params.coarsestDomain = crseDom;

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainSize[idir] = dLArray[idir];
    }

  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  for (int dir=0; dir<SpaceDim; dir++)
    {
      a_params.coarsestDx[dir] = a_params.domainSize[dir]/a_params.nCells[dir];
    }

  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainSize;


}

#include "NamespaceFooter.H"
