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
#include "ViscousTensorOp.H"
#include "ResistivityOp.H"
#include "ParmParse.H"
#include "FArrayBox.H"

namespace Chombo
{

  /// all data, no functions
  class VCPoissonParameters
  {
  public:
    //string        m_domainBC; always dirichlet
    int           m_maxLevel;
    int           m_numLevels;
    Vector<int>   m_refRatio;
    ProblemDomain m_coarsestDomain;
    Real          m_coarsestDx;
    Real          m_alpha;
    Real          m_beta;
  };

  ///sending in a bunch of stuff so I can keep functions consistent
  inline VCPoissonParameters
  getVCParameters(int                   a_finestLevel,
                  Real                  a_coarsestDx,
                  ProblemDomain         a_coarsestDomain,
                  Vector<int>           a_refRatio)
  {
    VCPoissonParameters  retval;

    ///have to use the same one that i use in the bc function
    retval.m_maxLevel       = a_finestLevel;
    retval.m_numLevels      = a_finestLevel + 1;
    retval.m_refRatio       = a_refRatio;
    retval.m_coarsestDomain = a_coarsestDomain;
    retval.m_coarsestDx     = a_coarsestDx;
    
    //solver to match other things
    ParmParse ppGetVCP("solver");
    ppGetVCP.get( "alpha",   retval.m_alpha);
    ppGetVCP.get( "beta",    retval.m_beta);
    return retval;
  }
  ///all functions, no data.
  class VCLocalFunctions
  {
  public:
    /// all static funcs, no data
    VCLocalFunctions(){;}
    ~VCLocalFunctions(){;}

    ///
    static void
    ParseBC(FArrayBox& a_state,
            const Box& a_valid,
            const ProblemDomain& a_domain,
            Real a_dx,
            bool a_homogeneous)
    {
      if (!a_domain.domainBox().contains(a_state.box()))
      {
        Box valid = a_valid;
        for(SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide  side = sit();
          for (int face_dir=0; face_dir < SpaceDim; face_dir++)
          {
            Box ghostBox = adjCellBox(valid, face_dir, sit(), 1);
            if (!a_domain.domainBox().contains(ghostBox))
            {
              HomogeneousDirichletBC(a_state,
                                     valid,
                                     face_dir,
                                     side);
            }  //end if there are cells outside domain on this side
          }    //end loop over face directions
        }      //end loop over sides
      }        //end if the data has ghost cells outside domain
    }          //end function ParseBC
    ///
    static  RealVect
    cellLocation(const IntVect & a_iv,
                 const Real    & a_dx )
    {
      RealVect retval;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        retval[idir] = (Real(a_iv[idir]) + 0.5)*a_dx;
      }
      return retval;
    }
    
    ///
    static   RealVect
    faceLocation(const IntVect&  a_iv,
                 const Real   &  a_dx ,
                 int             a_faceDir)
    {
      RealVect retval;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_faceDir == idir)
        {
          retval[idir] = (Real(a_iv[idir]))*a_dx;
        }
        else
        {
          retval[idir] = (Real(a_iv[idir]) + 0.5)*a_dx;
        }
      }
      return retval;
    }

    ///
    static 
    RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > >
    getConductivityFactory( Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_aCoef,
                            Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_bCoef,
                            Vector<DisjointBoxLayout>                    &   a_grids,
                            VCPoissonParameters                          &   a_params)
    {
      VCAMRPoissonOp2Factory* derived_ptr = (new VCAMRPoissonOp2Factory());

      derived_ptr->define(a_params.m_coarsestDomain,
                          a_grids,
                          a_params.m_refRatio,
                          a_params.m_coarsestDx,
                          &ParseBC,
                          a_params.m_alpha,
                          a_aCoef,
                          a_params.m_beta,
                          a_bCoef);

      AMRLevelOpFactory<LevelData<FArrayBox> >* base_ptr
        = static_cast<AMRLevelOpFactory<LevelData<FArrayBox> >*>(derived_ptr);
      return RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(base_ptr);
    } //end function getConductivityFactory

    ///
    static 
    RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > >
    getViscousTensorFactory( Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_acoef,
                             Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_eta,
                             Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_lam,
                             Vector<DisjointBoxLayout>                    &   a_grids,
                             VCPoissonParameters                          &   a_params)
    {
      ViscousTensorOpFactory* derived_ptr
        = (new ViscousTensorOpFactory(a_grids, a_eta, a_lam, a_acoef, 
                                      a_params.m_alpha,
                                      a_params.m_beta,
                                      a_params.m_refRatio,
                                      a_params.m_coarsestDomain,
                                      a_params.m_coarsestDx,
                                      &ParseBC));

      AMRLevelOpFactory<LevelData<FArrayBox> >* base_ptr
        = static_cast<AMRLevelOpFactory<LevelData<FArrayBox> >*>(  derived_ptr);
      return RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(base_ptr);
    } //end function getViscousTensorFactory


    ///
    static 
    RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > >
    getResistivityOpFactory( Vector<RefCountedPtr<LevelData<FArrayBox> > >&   a_acoef,
                             Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_eta,
                             Vector<RefCountedPtr<LevelData<FluxBox> > >  &   a_lam,
                             Vector<DisjointBoxLayout>                    &   a_grids,
                             VCPoissonParameters                          &   a_params)
    {
      ResistivityOpFactory* derived_ptr
        = (new ResistivityOpFactory(a_grids, a_eta, a_lam, a_acoef, 
                                      a_params.m_alpha,
                                      a_params.m_beta,
                                      a_params.m_refRatio,
                                      a_params.m_coarsestDomain,
                                      a_params.m_coarsestDx,
                                      &ParseBC));

      AMRLevelOpFactory<LevelData<FArrayBox> >* base_ptr
        = static_cast<AMRLevelOpFactory<LevelData<FArrayBox> >*>(  derived_ptr);
      return RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(base_ptr);
    } //end function getResistivityOpFactory

  };  //end class VCLocalFunctions
}     //end namespace chombo

#endif
