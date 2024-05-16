#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LevelMappedDerivatives.H"
#include "CellToEdge.H"
#include "DerivativesF_F.H"

#include "NamespaceHeader.H"

/** A set of utility functions for computing derivatives in the mapped
    sigma coordinate system used by ice-sheet models.
    if we represent (x,y,z) as (sigma, xtilde, ytilde) in 3d, where
    xtilde = x, ytilde = y, and sigma = (z_surface - z)/H, where H =
    ice-sheet thickness, then we have the following geometric
    relations: 

d/(dx) = d/(d xtilde) + (1/H)Delta_x*(d/d sigma)
d/(dy) = d/(d ytilde) + (1/H)Delta_y*(d/d sigma)
d/(d sigma) = -(1/H)d/(d sigma)

where Delta_x and Delta_y are geometric factors computed in the SigmaCS.

These utility functions implement this differencing so that it's
convenient to use. 

*/

/// compute cell-centered derivatives
/** it's assumed that any ghost-cell boundary conditions have already
    been set . optionally (a_mask = true), switch to one-sided differences  
    next to fluid free cells
*/
void 
computeCCDerivatives(LevelData<FArrayBox>& a_ccDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_derivGhost, 
		     const bool a_maskOneSide)
{
  const RealVect& dx = a_coordSys.dx();
  const LevelData<FArrayBox>& levelDeltas = a_coordSys.deltaFactors();
  const DisjointBoxLayout& grids = levelDeltas.getBoxes();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box derivBox = grids[dit];
      derivBox.grow(a_derivGhost);
      const FArrayBox& deltas = levelDeltas[dit];
      const FArrayBox& thisCCdata = a_ccData[dit];      

      CH_assert(deltas.box().contains(derivBox));
      
      for (int derivDir=a_derivDirections.begin(); derivDir<= a_derivDirections.end(); derivDir++)
        {
          // check to be sure that ccData is big enough
          Box growBox(derivBox);
          growBox.grow(derivDir,1);
          // also account for d/d(sigma) term here
          // (z-direction is component 0 in 3d)
          if ((CH_SPACEDIM == 3) && (derivDir != 0)) growBox.grow(0,1);
          CH_assert(thisCCdata.box().contains(growBox));
          
          for (int srcComp = a_comps.begin(); srcComp<=a_comps.end(); srcComp++)
            {
              int derivComp = derivComponent(derivDir, srcComp);
	      if (!a_maskOneSide)
		{
		  FORT_CCDERIV(CHF_FRA1(a_ccDeriv[dit],derivComp),
			       CHF_CONST_FRA1(thisCCdata, srcComp),
			       CHF_BOX(derivBox),
			       CHF_CONST_REAL(dx[derivDir]),
			       CHF_INT(derivDir));
		}
	      else
		{
		  FORT_CCDERIVMASK(CHF_FRA1(a_ccDeriv[dit],derivComp),
				   CHF_CONST_FRA1(thisCCdata, srcComp),
				   CHF_CONST_FRA1(a_coordSys.getH()[dit],0),
				   CHF_BOX(derivBox),
				   CHF_CONST_REAL(dx[derivDir]),
				   CHF_INT(derivDir));
		}

              if (CH_SPACEDIM == 3)
                {
                  // take d/d(sigma) terms into account)
                  MayDay::Warning("3D derivatives not implemented yet");
                }
            }      
        }
    }
}

/// compute face-centered derivatives of cell-centered data on all faces
void
computeFCDerivatives(LevelData<FluxBox>& a_fcDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_derivGhost)
{
  const RealVect& dx = a_coordSys.dx();
  const LevelData<FluxBox>& levelDeltas = a_coordSys.faceDeltaFactors();
  const DisjointBoxLayout& grids = levelDeltas.getBoxes();
  
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      const FluxBox& thisDeltas = levelDeltas[dit];
      const FArrayBox& thisCCdata = a_ccData[dit];
      FluxBox& thisFCderiv = a_fcDeriv[dit];

      Box ccDerivBox = gridBox;
      ccDerivBox.grow(a_derivGhost);

      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          FArrayBox& fcDeriv = thisFCderiv[faceDir];
          const FArrayBox& deltas = thisDeltas[faceDir];
          Box faceDerivBox(ccDerivBox);
          faceDerivBox.surroundingNodes(faceDir);
          CH_assert(deltas.box().contains(faceDerivBox));
  
          for (int derivDir=a_derivDirections.begin(); derivDir<= a_derivDirections.end(); derivDir++)
            {
              // check to be sure that ccData is big enough
              Box growBox(faceDerivBox);
              growBox.enclosedCells();
              growBox.grow(derivDir,1);
              // if transverse, also need to grow in faceDir
              if (derivDir != faceDir)
                {
                  growBox.grow(faceDir,1);
                }
              
              // also account for d/d(sigma) term here
              // (z-direction is component 0 in 3d)
              if ((CH_SPACEDIM == 3) && (derivDir != 0) && (faceDir != 0)) growBox.grow(0,1);
              CH_assert(thisCCdata.box().contains(growBox));
              
              for (int srcComp = a_comps.begin(); srcComp<=a_comps.end(); srcComp++)
                {
                  int derivComp = derivComponent(derivDir, srcComp);          
                  FORT_FACEDERIV(CHF_FRA1(fcDeriv,derivComp),
                                 CHF_CONST_FRA1(thisCCdata, srcComp),
                                 CHF_BOX(faceDerivBox),
                                 CHF_CONST_REAL(dx[derivDir]),
                                 CHF_INT(derivDir),
                                 CHF_INT(faceDir));
                  
                  if (CH_SPACEDIM == 3)
                    {
                      // take d/d(sigma) terms into account)
                      MayDay::Warning("3D derivatives not implemented yet");
                    }
                  
                } // end loop over components
            } // end loop over derivative direction
          
        } // end loop over face directions
    } //end loop over boxes
}

/// compute face-centered derivatives of cell-centered data
/** This version takes cell-centered derivatives in the mapped
    coordinates and averages them to faces and applies the correction
    required due to the mapping. This is essentially designed to
    support the horizontal ice-velocity elliptic operator.

    fcDeriv should be a face-centered FArrayBox

    ccData is the cell-centered data (needed on ghost cells normal 
           to the face in faceDir)

    ccDeriv contains the cell-centered derivatives used to compute
            tangential derivatives. These should be indexed in the
            same way as the cell-centered derivatives returned by the
            computeCCDerivatives function, and need to be defined on  
            one ghost cell in the faceDir direction.

    faceDerivBox should be a face-centered box corresponding to the 
                 region in fcDeriv where derivatives are needed.     
     
*/
void
computeFCDerivatives(LevelData<FluxBox>& a_fcDeriv,
                     const LevelData<FArrayBox>& a_ccData,
                     const LevelData<FArrayBox>& a_ccDeriv,
                     const LevelSigmaCS& a_coordSys,
                     const Interval& a_comps,
                     const Interval& a_derivDirections,
                     const IntVect& a_ghostVect)
{

  const RealVect& dx = a_coordSys.dx();
  const LevelData<FluxBox>& levelDeltas = a_coordSys.faceDeltaFactors();
  const DisjointBoxLayout& grids = levelDeltas.getBoxes();
  
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& thisCCdata = a_ccData[dit];
      const FArrayBox& thisCCderiv = a_ccDeriv[dit];
      FluxBox& thisFCderiv = a_fcDeriv[dit];
      const FluxBox& thisDeltas = levelDeltas[dit];
      Box cellDerivBox = grids[dit];
      cellDerivBox.grow(a_ghostVect);

      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {

          const FArrayBox& deltas = thisDeltas[faceDir];
          FArrayBox& fcDeriv = thisFCderiv[faceDir];
          Box faceDerivBox = cellDerivBox;
          faceDerivBox.surroundingNodes(faceDir);

          CH_assert(deltas.box().contains(faceDerivBox));
  
          for (int derivDir=a_derivDirections.begin(); derivDir<= a_derivDirections.end(); derivDir++)
            {
              // check to be sure that ccData and ccDerivData are big enough
              Box growBox(faceDerivBox);
              growBox.enclosedCells();
              growBox.grow(derivDir,1);
              
              // also account for d/d(sigma) term here
              // (z-direction is component 0 in 3d)
              if ((CH_SPACEDIM == 3) && (derivDir != 0) && (faceDir != 0)) growBox.grow(0,1);
              
              if (derivDir == faceDir)
                {
                  CH_assert(thisCCdata.box().contains(growBox));
                }
              else 
                {
                  CH_assert(thisCCderiv.box().contains(growBox));
                }
              
              for (int srcComp = a_comps.begin(); srcComp<=a_comps.end(); srcComp++)
                {
                  int derivComp = derivComponent(derivDir, srcComp);          
                  
                  // normal direction is the same as the other function
                  if (derivDir == faceDir)
                    {              
                      FORT_FACEDERIV(CHF_FRA1(fcDeriv,derivComp),
                                     CHF_CONST_FRA1(thisCCdata, srcComp),
                                     CHF_BOX(faceDerivBox),
                                     CHF_CONST_REAL(dx[derivDir]),
                                     CHF_INT(derivDir),
                                     CHF_INT(faceDir));
                      
                      if (CH_SPACEDIM == 3)
                        {
                          // take d/d(sigma) terms into account)
                          MayDay::Warning("3D derivatives not implemented yet");
                        }
                    }                 
                  else
                    {
                      // transverse direction derivative at this point is just
                      // a cell->face average
                      CellToEdge(thisCCderiv, derivComp,
                                 fcDeriv, derivComp,
                                 faceDir);
              
                      if (CH_SPACEDIM == 3)
                        {
                          // take d/d(sigma) terms into account)
                          MayDay::Warning("3D derivatives not implemented yet");
                        }
                    } // end if transverse direction
                } // end loop over components
            } // end loop over derivative direction
          
        } // end loop over face directions
    } // end loop over boxes

}
  


#include "NamespaceFooter.H"

