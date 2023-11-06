#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VCLOCALFUNCS_H_
#define _VCLOCALFUNCS_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "VCAMRPoissonOp2.H"
#include "ParmParse.H"
#include "FArrayBox.H"

namepsace Chombo
{

  /// all data, no functions
  class VCPoissonParameters
  {
  public:
    IntVect       nCells;
    int           maxGridSize;
    int           blockFactor;
    int           bufferSize;
    Real          fillRatio;
    Real          refineThresh;
    string        domainBC;
    int           maxLevel;
    int           numLevels;
    Vector<int>   refRatio;
    ProblemDomain coarsestDomain;
    Real          coarsestDx;
    RealVect      domainLength;
    RealVect      probLo;
    RealVect      probHi;
    Real          alpha;
    Real          beta;
    int           verbosity;
  };
  ///all functions, no data.
  class VCLocalFunctions
  {
  public:
    /// all static funcs, no data
    VCLocalFunctions(){;}
    ~VCLocalFunctions(){;}

    ///not sure this can be in a class
    static void
    ParseBC(FArrayBox& a_state,
            const Box& a_valid,
            const ProblemDomain& a_domain,
            Real a_dx,
            bool a_homogeneous)
    {
      string domBC;
      ParmParse pp;
      pp.get("domain_bc", domBC);
      Real bc_value = 0;
  
      if (!a_domain.domainBox().contains(a_state.box()))
      {

        Box valid = a_valid;
        for(SideIterator sit(); sit.ok(); ++sit)
        {
          for (int face_dir=0; face_dir < SpaceDim; face_dir++)
          {
            Box ghostBox = adjCellBox(valid, face_dir, sit(), 1);
            if (!a_domain.domainBox().contains(ghostBox))
            {
              if (domBC == string("Neumann"))
              {

                NeumBC(a_state,
                       valid,
                       a_dx,
                       a_homogeneous,
                       bc_value,
                       face_dir,
                       sit());
              }
              else if(domBC == string("Dirichlet"))
              {
                DiriBC(a_state,
                       valid,
                       a_dx,
                       a_homogeneous,
                       bc_value,
                       face_dir,
                       sit());

              }
              else
              {
                MayDay::Error("ParseBC:bogus DomBC");
              }
            }//end if there are cells outside domain on this side
          }    //end loop over face directions
        }      //end loop over sides
      }        //end if the data has ghost cells outside domain
    }            //end function ParseBC
    ///
    static  RealVect
    cellLocation(const IntVect & a_iv,
                 const RealVect& a_dx )
    {
      RealVect retval;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        retval[idir] = (Real(a_iv[idir]) + 0.5)*a_dx[idir];
      }
      return retval;
    }
    
    ///
    static   RealVect
    faceLocation(const IntVect&  a_iv,
                 const RealVect& a_dx ,
                 int             a_faceDir)
    {
      RealVect retval;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_faceDir == idir)
        {
          retval[idir] = (Real(a_iv[idir]))*a_dx[idir];
        }
        else
        {
          retval[idir] = (Real(a_iv[idir]) + 0.5)*a_dx[idir];
        }
      }
      return retval;
    }

    ///
    static  void
    setACoef(LevelData<FArrayBox>     & a_aCoef,
             const VCPoissonParameters& a_params,
             const RealVect           & a_dx)
    {
      int num;
      DataIterator dit = a_aCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& aCoef = a_aCoef[dit];
        for(BoxIterator boxit(aCoef.box()); boxit.ok(); ++boxit)
        {
          auto iv = boxit();
          RealVect pos = cellLocation(iv, a_dx);
          aCoef(iv, 0) = pos[0];
        } //end loop over cells
      }   // end loop over grids
    }     //end function setACoef

    ///
    static void
    setBCoef(LevelData<FluxBox>       & a_bCoef,
             const VCPoissonParameters& a_params,
             const RealVect           & a_dx)
    {
      DataIterator dit = a_bCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = a_bCoef[dit];
        for (int face_dir=0; face_dir < SpaceDim; face_dir++)
        {
          FArrayBox& bcofab = thisBCoef[dir];
          const Box& dirBox = dirFlux.box();
          RealVect pos;
          for(BoxIterator boxit(dirBox); boxit.ok(); ++boxit)
          {
            IntVect   iv  = boxit();
            RealVect  pos = faceLocation(iv, a_dx, face_dir);
            //bcoef = x + y + z in the amrpoissonop example
            Real sum = pos.sum();
            bcofab(iv,  0) = sum;
          } //end loop over cells
        }   //end loop over face direction
      }     //end loop over boxes
    }       //end function setBCoef

    /********/
    static void
    getPoissonParameters(VCPoissonParameters&  a_params)
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
      pp.get("alpha",a_params.alpha);
      pp.get("beta", a_params.beta);

      pp.get("max_level", a_params.maxLevel);
      a_params.numLevels = a_params.maxLevel + 1;
      pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
      a_params.verbosity = 3;
      pp.query("verbosity", a_params.verbosity);

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
        a_params.domainLength[idir] = dLArray[idir];
      }

      pp.get("max_grid_size",a_params.maxGridSize);

      //derived stuff
      a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];

      a_params.probLo = RealVect::Zero;
      a_params.probHi = RealVect::Zero;
      a_params.probHi += a_params.domainLength;
    }

    /********/
    static void setRHS(LevelData<FArrayBox>&    a_rhs,
                       const RealVect&          a_dx,
                       const VCPoissonParameters& a_params)
    {
      CH_assert(a_rhs.nComp() == 1);

      int comp = 0;
      const RealVect&  trig = getTrigRV();

      for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisRHS = a_rhs[dit()];
        thisRHS.setVal(1.);

      }
    }


    static 
    RefCountedPtr< VCAMRPoissonOp2Factory >
    defineOperatorFactory(const Vector<DisjointBoxLayout>                    &   a_grids,
                          const Vector<ProblemDomain>                        &   a_vectDomain,
                          const Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_aCoef,
                          const Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_bCoef,
                          const VCPoissonParameters                          &   a_params)
    {
      RefCountedPtr< VCAMRPoissonOp2Factory > retval(new VCAMRPoissonOp2Factory());

      retval->define(a_params.coarsestDomain,
                     a_grids,
                     a_params.refRatio,
                     a_params.coarsestDx,
                     &ParseBC,
                     a_params.alpha,
                     a_aCoef,
                     a_params.beta,
                     a_bCoef);

      return RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >( opFactory);
    }
  };
}
#endif
