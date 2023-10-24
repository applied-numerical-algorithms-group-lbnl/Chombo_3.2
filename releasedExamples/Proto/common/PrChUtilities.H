#ifndef __PR_CH_UTILIITIES_
#define __PR_CH_UTILIITIES_
#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "Proto.H"
#include "memusage.H"

namespace Chombo
{
  ///
  /**
     This may look like yet another class just to keep the syntax clean and it may be that.
     This class is also a pretty wrapper around proto/include/transition/Proto_Transition_Utilities (PTU).   
     PTU on this side of the fence is a plain-old-data interface.   
     This class prettifies that somewhat and provides some code reuse for the test codes.
   **/
  temmplate <int ncomp>
  class PrChUtilities
  {
  public:
    ///no data. all static functions
    PrChUtilities(){;}
    ~PrChUtilities(){;}
    typedef        LevelData<FArrayBox>          ch_ldf;
    typedef        ProblemDomain                 ch_dom;
    typedef Proto::ProblemDomain                 pr_dom;
    typedef        DisjointBoxLayout             ch_dbl;
    typedef Proto::DisjointBoxLayout             pr_dbl;
    
    ///just always tag IntVect::Zero;
    static void
    setupLLCornerGrids(Vector<ch_dbl>& a_amrGrids,
                       Vector<ch_dom>& a_amrDomains,
                       Vector<int>& a_refRatios,
                       Vector<Real>& a_amrDx,
                       int& a_finestLevel)
    {
      CH_TIME("PrChUtilities::setupLLCornerGrids");

      a_finestLevel = 0;
      ParmParse ppGrids("setupLLCornerGrids");

      // get grid generation parameters
      int maxLevel, maxBoxSize, blockFactor;
      Real fillRatio;

      ppGrids.get("max_level", maxLevel);

      ppGrids.get("max_box_size",maxBoxSize);

      ppGrids.get("block_factor", blockFactor);

      ppGrids.get("fillRatio", fillRatio);
      int maxNumLevels = maxLevel +1;

      a_refRatios.resize(maxNumLevels);
      ppGrids.getarr("ref_ratio", a_refRatios, 0, maxNumLevels);

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
      a_amrGrids.resize(maxNumLevels);
      a_amrDomains.resize(maxNumLevels);
      a_amrDx.resize(maxNumLevels,-1);
      a_finestLevel = 0;

      // assumes dx=dy=dz
      a_amrDx[0] = domainSize[0]/numCells[0];

      IntVect domLo = IntVect::Zero;
      IntVect domHi  = numCells - IntVect::Unit;

      ch_dom baseDomain(domLo, domHi, is_periodic);
      a_amrDomains[0] = baseDomain;

      // set up refined domains, etc
      for (int lev=1; lev<= maxLevel; lev++)
      {
        a_amrDomains[lev] = a_amrDomains[lev-1];
        a_amrDomains[lev].refine(a_refRatios[lev-1]);
        a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
      }

  
      Vector<Vector<Box> > vectBoxes(maxLevel+1);
      {
        CH_TIME("setupGrids: BaseGridCreation");
        // generate base level grids

        domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

        Vector<int> procAssign(vectBoxes[0].size(), 0);

        LoadBalance(procAssign, vectBoxes[0]);

        ch_dbl baseGrids(vectBoxes[0], procAssign, baseDomain);

        a_amrGrids[0] = baseGrids;
      }


      if (maxLevel > 0)
      {
        CH_TIME("setupGrids: meshRefine and all that");
        int bufferSize = 0;
        BRMeshRefine meshGen(a_amrDomains[0],
                             a_refRatios,
                             fillRatio,
                             blockFactor,
                             bufferSize,
                             maxBoxSize);
      
        // to be used by MeshRefine...
        Vector<Vector<Box> > old_meshes(maxLevel+1);
        old_meshes[0] = vectBoxes[0];
        for (int lev=1; lev< old_meshes.size(); lev++)
        {
          Box whole_dom = a_amrDomains[lev].domainBox();
          old_meshes[lev] = Vector<Box>(1, whole_dom);
        }

        //just always tag IntVect::Zero;
        IntVectSet taglev(IntVect::Zero);
        Vector<IntVectSet> tags(a_finestLevel+1, taglev);
        Vector<Vector<Box> > new_meshes;

        meshGen.regrid(new_meshes, tags,  0, a_finestLevel,  old_meshes);

        for (int ilev=0; ilev < a_finestLevel; ilev++)
        {
          const auto& boxes = new_meshes[ilev];
          Vector<int> procs;
          LoadBalance(procs, boxes);
          a_amrGrids[ilev] = ch_dbl(boxes, procs);
        }
      }
    } //end setupLLCornerGrids


    ///always tag domain highest point (nx-1, ny-1, nz-1)
    static void
    setupURCornerGrids(Vector<ch_dbl>& a_amrGrids,
                       Vector<ch_dom>& a_amrDomains,
                       Vector<int>& a_refRatios,
                       Vector<Real>& a_amrDx,
                       int& a_finestLevel)
    {
      CH_TIME("PrChUtilities::setupURCornerGrids");

      a_finestLevel = 0;
      ParmParse ppGrids("setupURCornerGrids");

      // get grid generation parameters
      int maxLevel, maxBoxSize, blockFactor;
      Real fillRatio;

      ppGrids.get("max_level", maxLevel);

      ppGrids.get("max_box_size",maxBoxSize);

      ppGrids.get("block_factor", blockFactor);

      ppGrids.get("fillRatio", fillRatio);
      int maxNumLevels = maxLevel +1;

      a_refRatios.resize(maxNumLevels);
      ppGrids.getarr("ref_ratio", a_refRatios, 0, maxNumLevels);

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
      a_amrGrids.resize(maxNumLevels);
      a_amrDomains.resize(maxNumLevels);
      a_amrDx.resize(maxNumLevels,-1);
      a_finestLevel = 0;

      // assumes dx=dy=dz
      a_amrDx[0] = domainSize[0]/numCells[0];

      IntVect domLo = IntVect::Zero;
      IntVect domHi  = numCells - IntVect::Unit;

      ch_dom baseDomain(domLo, domHi, is_periodic);
      a_amrDomains[0] = baseDomain;

      // set up refined domains, etc
      for (int lev=1; lev<= maxLevel; lev++)
      {
        a_amrDomains[lev] = a_amrDomains[lev-1];
        a_amrDomains[lev].refine(a_refRatios[lev-1]);
        a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
      }

  
      Vector<Vector<Box> > vectBoxes(maxLevel+1);
      {
        CH_TIME("setupGrids: BaseGridCreation");
        // generate base level grids

        domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

        Vector<int> procAssign(vectBoxes[0].size(), 0);

        LoadBalance(procAssign, vectBoxes[0]);

        ch_dbl baseGrids(vectBoxes[0], procAssign, baseDomain);

        a_amrGrids[0] = baseGrids;
      }


      if (maxLevel > 0)
      {
        CH_TIME("setupGrids: meshRefine and all that");
        int bufferSize = 0;
        BRMeshRefine meshGen(a_amrDomains[0],
                             a_refRatios,
                             fillRatio,
                             blockFactor,
                             bufferSize,
                             maxBoxSize);
      
        // to be used by MeshRefine...
        Vector<Vector<Box> > old_meshes(maxLevel+1);
        old_meshes[0] = vectBoxes[0];
        for (int lev=1; lev< old_meshes.size(); lev++)
        {
          Box whole_dom = a_amrDomains[lev].domainBox();
          old_meshes[lev] = Vector<Box>(1, whole_dom);
        }

        //just always tag (nx-1, ny-1, nz-1)
        Vector<IntVectSet>   tags(maxNumLevels);
        for(int ilev = 0; ilev <  maxNumLevels; ilev++)
        {
          const auto&  wholeDom = a_amrDomains[ilev];
          IntVect highpt = wholeDom.domainBox().bigEnd();
          tags[ilev] = IntVectSet(highpt);
        }
        Vector<Vector<Box> > new_meshes;

        meshGen.regrid(new_meshes, tags,  0, a_finestLevel,  old_meshes);

        for (int ilev=0; ilev < a_finestLevel; ilev++)
        {
          const auto& boxes = new_meshes[ilev];
          Vector<int> procs;
          LoadBalance(procs, boxes);
          a_amrGrids[ilev] = ch_dbl(boxes, procs);
        }
      }
    }   //end setupURCornerGrids
  };    //end class PrChUtilties
}       //end namespace chombo
#endif //include guard
