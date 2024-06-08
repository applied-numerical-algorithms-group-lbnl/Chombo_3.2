#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "ExternalObj.H"
#include "ExternalObjF_F.H"

#include "NamespaceHeader.H"

ExternalObj::ExternalObj()
{}

ExternalObj::~ExternalObj()
{}

void ExternalObj::NonLinear_level(LevelData<FArrayBox>&        a_NL,
                                  LevelData<FArrayBox>&        a_dNL,
                                  const LevelData<FArrayBox>&  a_u)
{
  CH_TIME("setNL_Level");

  //ParmParse pp("solver");

  Real gamma = 100.0;
  //pp.query("gamma", gamma);

  DataIterator levelDit        = a_NL.dataIterator();
  const DisjointBoxLayout& dbl = a_NL.disjointBoxLayout();
  for (levelDit.begin(); levelDit.ok(); ++levelDit) {

      const Box& region = dbl[levelDit];

      FORT_COMPUTENONLINEARTERMS( CHF_FRA(a_u[levelDit]),
            CHF_CONST_REAL(gamma),
            CHF_BOX(region),
            CHF_FRA(a_NL[levelDit]),
            CHF_FRA(a_dNL[levelDit]) );
  
  } // end loop over grids on this level
}


void ExternalObj::NonLinear_piece(Vector<LevelData<FArrayBox>* > a_NL,
                                  Vector<LevelData<FArrayBox>* > a_dNL,
                                  Vector<LevelData<FArrayBox>* > a_u,
                                  int a_finestLevel)
{
  CH_TIME("setNL_piece");

  //ParmParse pp("solver");

  Real gamma = 100.0;
  //pp.query("gamma", gamma);

  for (int lev=0; lev<=a_finestLevel; lev++) {

    LevelData<FArrayBox>& levelNL   = *(a_NL[lev]);
    LevelData<FArrayBox>& leveldNL  = *(a_dNL[lev]);
    LevelData<FArrayBox>& levelU    = *(a_u[lev]);

    DataIterator levelDit = levelNL.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit) {

        FArrayBox& thisNL    = levelNL[levelDit];
        FArrayBox& thisdNL   = leveldNL[levelDit];
        FArrayBox& thisU     = levelU[levelDit];

        BoxIterator bit(thisNL.box());
        for (bit.begin(); bit.ok(); ++bit) {

            IntVect iv = bit();

            thisNL(iv, 0)  = ( gamma * thisU(iv, 0) * exp( thisU(iv, 0) ) );
            thisdNL(iv, 0) = ( gamma * (1 + thisU(iv, 0)) * exp( thisU(iv, 0) ) );
        }
    } // end loop over grids on this level
  } // end loop over levels
}

#include "NamespaceFooter.H" 
