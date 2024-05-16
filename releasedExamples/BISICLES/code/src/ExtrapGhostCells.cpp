#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ExtrapGhostCells.H"
#include "ExtrapBCF_F.H"

#include "NamespaceHeader.H"

void 
ExtrapGhostCells(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain)
{
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      ExtrapGhostCells(a_phi[dit], a_domain, a_phi.ghostVect());
    }
}

void 
ExtrapGhostCells(FArrayBox& a_phi,
                 const ProblemDomain& a_domain,
                 const IntVect& a_ghostVect)
{

  const Box& domainBox = a_domain.domainBox();

  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (!a_domain.isPeriodic(dir))
        {
          int rad = a_ghostVect[dir];
          
          // lo-side 
           int hiLo = 0;
	  // slc :: we need to do the low side one strip at 
	  // at time so that strip n is filled before we
          // extrapolate from it into strip n -1. 
          Box  stripLo = adjCellLo(domainBox,  dir, 1);
          // do this to try to catch corner cells
	  stripLo.grow(rad);
          stripLo.grow(dir,-rad);
	  for (int i =0; i < rad; ++i)
	    {
	      Box ghostBoxLo = stripLo;
	      ghostBoxLo &= a_phi.box();
	      if (!ghostBoxLo.isEmpty())
		{
		  FORT_SIMPLEEXTRAPBC(CHF_FRA(a_phi),
				      CHF_BOX(ghostBoxLo),
				      CHF_INT(dir), 
				      CHF_INT(hiLo));
		}
	      stripLo.shift(dir,-1);
	    }
          
          // hi-side
          hiLo = 1;
          Box ghostBoxHi = adjCellHi(domainBox, 
                                     dir, rad);
          // do this to try to catch corner cells
          ghostBoxHi.grow(1);
          ghostBoxHi.grow(dir,-1);
          ghostBoxHi &= a_phi.box();
          if(!ghostBoxHi.isEmpty())
            {
              FORT_SIMPLEEXTRAPBC(CHF_FRA(a_phi),
                                  CHF_BOX(ghostBoxHi),
                                  CHF_INT(dir), 
                                  CHF_INT(hiLo));
            }
          
        } // end if not periodic in this direction
    } // end loop over directions
}

#include "NamespaceFooter.H"


