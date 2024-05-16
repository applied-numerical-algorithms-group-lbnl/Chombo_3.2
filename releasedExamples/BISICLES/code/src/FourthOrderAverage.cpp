#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FourthOrderAverage.H"
#include "FourthOrderAverageF_F.H"

#include "NamespaceHeader.H"

/// -----------------------------------------------------------
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        int a_sgn)
{
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];

      fourthOrderAverageCell(thisPhi, a_sgn);
    } // end loop over boxes
}


/// -----------------------------------------------------------
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        const ProblemDomain&  a_domain,
                        int a_sgn)
{
  const DisjointBoxLayout& layout = a_phi.disjointBoxLayout();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      const Box& bx = layout[dit];

      fourthOrderAverageCell(thisPhi, a_domain, bx, a_sgn);
    } // end loop over boxes
}


/// -----------------------------------------------------------
void fourthOrderAverageCell(FArrayBox& a_phi,
                            int a_sgn)
{
  // derivBox shrunk by 1
  Box derivBox(a_phi.box());
  derivBox.grow(-1);

  FArrayBox tempLap(derivBox, a_phi.nComp());
  tempLap.setVal(0.0);

  Real factor = a_sgn * (1.0/24.0);
  for (int dir=0; dir <SpaceDim; dir++)
    {
      FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                              CHF_CONST_FRA(a_phi),
                              CHF_BOX(derivBox),
                              CHF_CONST_INT(dir),
                              CHF_CONST_REAL(factor));
    } // end loop over directions
  a_phi.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
}


/// -----------------------------------------------------------
void fourthOrderAverageCell(FArrayBox& a_phi,
                            const ProblemDomain&  a_domain,
                            const Box&            a_bx,
                            int a_sgn)
{
  CH_assert(a_phi.box().contains(a_bx));
  Box derivBox(a_bx);
  derivBox.grow(1);
  derivBox &= a_domain;
  CH_assert(a_phi.box().contains(derivBox));

  FArrayBox tempLap(a_bx, a_phi.nComp());
  tempLap.setVal(0.0);

  Real factor = a_sgn * (1.0/24.0);
  for (int dir=0; dir <SpaceDim; dir++)
    {
      Box derivDirBox(a_bx);
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().smallEnd(dir) == a_bx.smallEnd(dir)))
        {
          // one-sided derivative, low end
          Box loBox(a_bx);
          loBox.setRange(dir, a_bx.smallEnd(dir));
          FORT_INCREMENTLOSIDELAPLACIAN(CHF_FRA(tempLap),
                                        CHF_FRA(a_phi),
                                        CHF_BOX(loBox),
                                        CHF_INT(dir),
                                        CHF_REAL(factor));
          derivDirBox.growLo(dir, -1);
        }
      if (!a_domain.isPeriodic(dir) &&
          (a_domain.domainBox().bigEnd(dir) == a_bx.bigEnd(dir)))
        {
          // one-sided derivative, high end
          Box hiBox(a_bx);
          hiBox.setRange(dir, a_bx.bigEnd(dir));
          FORT_INCREMENTHISIDELAPLACIAN(CHF_FRA(tempLap),
                                        CHF_FRA(a_phi),
                                        CHF_BOX(hiBox),
                                        CHF_INT(dir),
                                        CHF_REAL(factor));
          derivDirBox.growHi(dir, -1);
        }
      FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                              CHF_CONST_FRA(a_phi),
                              CHF_BOX(derivDirBox),
                              CHF_CONST_INT(dir),
                              CHF_CONST_REAL(factor));
    } // end loop over directions
  a_phi.plus(tempLap, a_bx, 0, 0, tempLap.nComp());
}


#include "NamespaceFooter.H"
