#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "convert3Dto2Dplot.H"

#define SAVE_SPACEDIM CH_SPACEDIM
#include "LevelData.H.multidim"
#include "FArrayBox.H.multidim"
#include "Vector.H"
#include "fabncio.H"
#include "AMRIO.H.multidim"

// At this point, CH_SPACEDIM=0; the .H.multidim files do that.

#include "Slicing.H.transdim"
#include "Injection.H.transdim"

//#define CH_SPACEDIM SAVE_SPACEDIM
//#undef SAVE_SPACEDIM
//#include "UsingNamespace.H"

using namespace Chombo;
using namespace CH_MultiDim;
using std::endl;

int convert3DNCTo2DPlot(const string& fname3d,
                        const string& fname2d)
{
  
  int status =0;

  // first, define 3D dataholders
  D3::Vector<D3::DisjointBoxLayout> vectGrids3d;
  D3::Vector<D3::LevelData<D3::FArrayBox>* >  vectData3d;
  Vector<string> vectNames;
  D3::Box domain3d;
  Real dx;
  Real dt;
  Real time;
  Vector<int> refRatio;
  int numLevels;
#if 0
  status = D3::ReadAMRHierarchyHDF5(fname3d,
                                    vectGrids3d,
                                    vectData3d,
                                    vectNames,
                                    domain3d,
                                    dx,
                                    dt,
                                    time,
                                    refRatio,
                                    numLevels);
#endif

  Vector<std::string> var;

  D3::FArrayBox fab;
  D3::Box box;
  if (procID() == uniqueProc(SerialTask::compute))
      {
#ifdef HAVE_NETCDF
	NCIO::readFAB(in_file,var,fab,dx);
#else
	MayDay::Error("netcdf input requested but netcdf support not built")
#endif
	box = fab.box();
      }// end if serial compute
    broadcast(box,uniqueProc(SerialTask::compute));
    D3::ProblemDomain pd(box);
    D3::Vector<D3::Box> boxes(1,pd.domainBox());
    Vector<int> procAssign(1,uniqueProc(SerialTask::compute));
    D3::DisjointBoxLayout grids(boxes, procAssign, pd);
    D3::LevelData<D3::FArrayBox> levelData(grids,var.size(),D3::IntVect::Zero);

    for (D3::DataIterator dit(grids); dit.ok(); ++dit)
      {
	levelData[dit].copy(fab,0,0,fab.nComp());
      }


  // now, define 2D dataholders
  D2::Vector<D2::DisjointBoxLayout> vectGrids2d(vectGrids3d.size());
  D2::Vector<D2::LevelData<D2::FArrayBox>* >  vectData2d(vectData3d.size(), NULL);
;
  D2::Box domain2d;
  
  // now do the magic slicing from 3d->2d
  D3::SliceSpec slice(2,0);
  sliceBox(domain2d, domain3d, slice);
  for (int lev=0; lev<vectGrids3d.size(); lev++)
    {
      sliceDisjointBoxLayout(vectGrids2d[lev],
                             vectGrids3d[lev],
                             slice);
      vectData2d[lev] = new D2::LevelData<D2::FArrayBox>;
      sliceLevelData(*vectData2d[lev], *vectData3d[lev], slice);
    }

  D2::WriteAMRHierarchyHDF5(fname2d,
                            vectGrids2d,
                            vectData2d,
                            vectNames,
                            domain2d,
                            dx,
                            dt,
                            time,
                            refRatio,
                            numLevels);
  


  return status;
}
