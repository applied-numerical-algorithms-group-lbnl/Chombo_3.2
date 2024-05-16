#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "Averaging.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "AverageF_F.H"
#include "AverageFaceF_F.H"

#include "NamespaceHeader.H"

/** A set of utility functions for averaging data from a fine grid to
    a coarse grid
*/

void 
horizontalAverage(LevelData<FArrayBox>& a_crseData,
                  const LevelData<FArrayBox>& a_fineData,
                  int a_nRef)
{
  const DisjointBoxLayout& crseGrids = a_crseData.getBoxes();
  const DisjointBoxLayout& fineGrids = a_fineData.getBoxes();
  
  // if crseGrids is a simple coarsening of the fine grids, then 
  // both can use the same DataIterator...
  if (crseGrids.compatible(fineGrids))
    {
      DataIterator dit = fineGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& fineFab = a_fineData[dit];
          FArrayBox& crseFab = a_crseData[dit];
          // figure out how much we can average
          Box crseBox = fineFab.box();
          crseBox.coarsen(a_nRef);
          crseBox &= crseFab.box();

          // now call single-fab version
          horizontalAverage(crseFab, fineFab, crseBox, a_nRef);
        }
    } // end if crse and fine grids compatible
  else
    {
      // bail for now -- implement this later if we have to
      MayDay::Error("LevelData<FArrayBox> horizontalAverage not implemented for incompatible crse and fine grids");
    }                
}

void
horizontalAverageFace(LevelData<FluxBox>& a_crseData,
                      const LevelData<FluxBox>& a_fineData,
                      int a_nRef)
{
  const DisjointBoxLayout& crseGrids = a_crseData.getBoxes();
  const DisjointBoxLayout& fineGrids = a_fineData.getBoxes();
  
  // if crseGrids is a simple coarsening of the fine grids, then 
  // both can use the same DataIterator...
  if (crseGrids.compatible(fineGrids))
    {
      DataIterator dit = fineGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FluxBox& fineFlux = a_fineData[dit];
          FluxBox& crseFlux = a_crseData[dit];
          // figure out how much we can average
          Box crseBox = fineFlux.box();
          crseBox.coarsen(a_nRef);
          crseBox &= crseFlux.box();

          // now call single-fab version
          horizontalAverageFace(crseFlux, fineFlux, crseBox, a_nRef);
        }
    } // end if crse and fine grids compatible
  else
    {
      // bail for now -- implement this later if we have to
      MayDay::Error("LevelData<FluxBox> horizontalAverageFace not implemented for incompatible crse and fine grids");
    }                

}

void
averageAllDim(LevelData<FArrayBox>& a_crseData, 
              const LevelData<FArrayBox>& a_fineData,
              int a_nRef)
{
  const DisjointBoxLayout& crseGrids = a_crseData.getBoxes();
  const DisjointBoxLayout& fineGrids = a_fineData.getBoxes();
  
  // if crseGrids is a simple coarsening of the fine grids, then 
  // both can use the same DataIterator...
  if (crseGrids.compatible(fineGrids))
    {
      DataIterator dit = fineGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& fineFab = a_fineData[dit];
          FArrayBox& crseFab = a_crseData[dit];
          // figure out how much we can average
          Box crseBox = fineFab.box();
          crseBox.coarsen(a_nRef);
          crseBox &= crseFab.box();

          // now call single-fab version
          averageAllDim(crseFab, fineFab, crseBox, a_nRef);
        }
    } // end if crse and fine grids compatible
  else
    {
      // bail for now -- implement this later if we have to
      MayDay::Error("LevelData<FArrayBox> averageAllDim not implemented for incompatible crse and fine grids");
    }                

}


void
averageAllDimFace(LevelData<FluxBox>& a_crseData, 
                  const LevelData<FluxBox>& a_fineData,
                  int a_nRef)
{
  const DisjointBoxLayout& crseGrids = a_crseData.getBoxes();
  const DisjointBoxLayout& fineGrids = a_fineData.getBoxes();
  
  // if crseGrids is a simple coarsening of the fine grids, then 
  // both can use the same DataIterator...
  if (crseGrids.compatible(fineGrids))
    {
      DataIterator dit = fineGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FluxBox& fineFlux = a_fineData[dit];
          FluxBox& crseFlux = a_crseData[dit];
          // figure out how much we can average
          Box crseBox = fineFlux.box();
          crseBox.coarsen(a_nRef);
          crseBox &= crseFlux.box();

          // now call single-fab version
          averageAllDimFace(crseFlux, fineFlux, crseBox, a_nRef);
        }
    } // end if crse and fine grids compatible
  else
    {
      // bail for now -- implement this later if we have to
      MayDay::Error("LevelData<FluxBox> averageAllDimFace not implemented for incompatible crse and fine grids");
    }                



}

// --------------------------
//  Single-box versions
// --------------------------


void
horizontalAverage(FArrayBox& a_crseData,
                  const FArrayBox& a_fineData,
                  const Box& a_crseBox,
                  int a_nRef)
{
  if (SpaceDim == 2)
    {
      // in 2D, this is just the normal averaging one does every day...
      Box refbox(IntVect::Zero,
                 (a_nRef-1)*IntVect::Unit);
      // assume arithmetic averaging here
      
      FORT_AVERAGE( CHF_FRA(a_crseData),
                    CHF_CONST_FRA(a_fineData),
                    CHF_BOX(a_crseBox),
                    CHF_CONST_INT(a_nRef),
                    CHF_BOX(refbox)
                    );
      // defer 3D for now
    }
  else if (SpaceDim == 3)
    {
      MayDay::Error("horizontalAverage not defined for 3d");
    }
}

void
horizontalAverageFace(FluxBox& a_crseData,
                      const FluxBox& a_fineData,
                      const Box& a_crseCellBox,
                      int a_nRef)
{
  if (SpaceDim == 2)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& coarseFab = a_crseData[dir];
          const FArrayBox& fineFab = a_fineData[dir];
          
          Box crseFaceBox(a_crseCellBox);
          crseFaceBox.surroundingNodes(dir);

          // set up refinement box
          int boxHi = a_nRef-1;
          IntVect hiVect(D_DECL6(boxHi,boxHi,boxHi,
                                 boxHi,boxHi,boxHi));
          // don't want to index at all in dir direction --
          // instead, want to just march along face.
          hiVect.setVal(dir,0);
          IntVect loVect(D_DECL6(0,0,0,0,0,0));
          Box refBox(loVect, hiVect);
          
          // assume arithmetic averaging here
          FORT_AVERAGEFACE( CHF_FRA(coarseFab),
                            CHF_CONST_FRA(fineFab),
                            CHF_BOX(crseFaceBox),
                            CHF_CONST_INT(dir),
                            CHF_CONST_INT(a_nRef),
                            CHF_CONST_INT(a_nRef),
                            CHF_BOX(refBox));
        }
    } 
  // defer 3D for now
  else if (SpaceDim == 3)
    {
      MayDay::Error("horizontalAverageFace not defined for 3D");
    }
}

void
averageAllDim(FArrayBox& a_crseData, 
              const FArrayBox& a_fineData,
              const Box& a_crseBox,
              int a_nRef)
{
  Box refbox(IntVect::Zero,
             (a_nRef-1)*IntVect::Unit);
  // assume arithmetic averaging here
  
  FORT_AVERAGE( CHF_FRA(a_crseData),
                CHF_CONST_FRA(a_fineData),
                CHF_BOX(a_crseBox),
                CHF_CONST_INT(a_nRef),
                CHF_BOX(refbox)
                );
  
}


void
averageAllDimFace(FluxBox& a_crseData, 
                  const FluxBox& a_fineData,
                  const Box& a_crseCellBox,
                  int a_nRef)
{
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& coarseFab = a_crseData[dir];
      const FArrayBox& fineFab = a_fineData[dir];
      
      Box crseFaceBox(a_crseCellBox);
      crseFaceBox.surroundingNodes(dir);
      
      // set up refinement box
      int boxHi = a_nRef-1;
      IntVect hiVect(D_DECL6(boxHi,boxHi,boxHi,
                             boxHi,boxHi,boxHi));
      // don't want to index at all in dir direction --
      // instead, want to just march along face.
      hiVect.setVal(dir,0);
      IntVect loVect(D_DECL6(0,0,0,0,0,0));
      Box refBox(loVect, hiVect);
      
      // assume arithmetic averaging here
      FORT_AVERAGEFACE( CHF_FRA(coarseFab),
                        CHF_CONST_FRA(fineFab),
                        CHF_BOX(crseFaceBox),
                        CHF_CONST_INT(dir),
                        CHF_CONST_INT(a_nRef),
                        CHF_CONST_INT(a_nRef),
                        CHF_BOX(refBox));
    }
  
}


#include "NamespaceFooter.H"
