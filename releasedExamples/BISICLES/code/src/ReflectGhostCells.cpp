#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ReflectGhostCells.H"
#include "NamespaceHeader.H"

void 
ReflectGhostCells(LevelData<FArrayBox>& a_phi,
		  const ProblemDomain& a_domain,
		  const int a_dir,
		  const Side::LoHiSide a_side)
{
  if (!a_domain.isPeriodic(a_dir))
    {
      DataIterator dit = a_phi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  ReflectGhostCells(a_phi[dit], a_domain, a_phi.ghostVect(), a_dir, a_side);
	}
    }
}

void 
ReflectGhostCells(FArrayBox& a_phi,
		  const ProblemDomain& a_domain,
		  const IntVect& a_ghostVect,
		  const int a_dir, 
		  const Side::LoHiSide a_side)
{

  const Box& domainBox = a_domain.domainBox();
  int rad = a_ghostVect[a_dir]; 
  Box  stripGhost = adjCellBox(domainBox,  a_dir, a_side, 1);
  // do this to try to catch corner cells
  stripGhost.grow(rad);
  stripGhost.grow(a_dir,-rad);
  
  for (int i =0; i < rad; ++i)
    {
      Box toBox = stripGhost;
      toBox &= a_phi.box();
      Box fromBox = toBox;
      fromBox.shift(a_dir,-sign(a_side) * (2*i + 1));
      a_phi.copy(a_phi,fromBox,0,toBox,0,a_phi.nComp());
      stripGhost.shift(a_dir,sign(a_side));
    }
}
   

#include "NamespaceFooter.H"

