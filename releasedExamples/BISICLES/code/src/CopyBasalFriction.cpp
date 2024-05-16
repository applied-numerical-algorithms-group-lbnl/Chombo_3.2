#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "FIBCF_F.H"
#include "FillFromReference.H"
#include "CopyBasalFriction.H"
CopyBasalFriction::CopyBasalFriction()
{
  m_verbose = false;
}


CopyBasalFriction::CopyBasalFriction(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_vectBeta,
                                     Vector<ProblemDomain>& a_vectDdomains)
{
  define(a_vectBeta, a_vectDdomains);
}


void
CopyBasalFriction::define(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_vectBeta,
                          Vector<ProblemDomain>& a_vectDomains)
{
  m_vectBeta = a_vectBeta;
  m_vectDomains = a_vectDomains;
}

BasalFriction* CopyBasalFriction::new_basalFriction() const
{
  
  CopyBasalFriction* ptr = new CopyBasalFriction();
  ptr->m_vectBeta = m_vectBeta;
  ptr->m_vectDomains = m_vectDomains;
  ptr->m_verbose = m_verbose;
  return static_cast<BasalFriction*>(ptr);
}


void
CopyBasalFriction::setBasalFriction(LevelData<FArrayBox>& a_C,
                                    LevelSigmaCS& a_coordSys,
                                    Real a_time,
                                    Real a_dt)
{
  const DisjointBoxLayout& levelGrids = a_C.getBoxes();
  const ProblemDomain& levelDomain = levelGrids.physDomain();

  bool foundLevel = false;
  for (int lev=0; lev<m_vectBeta.size(); lev++)
    {
      if (levelDomain.domainBox() == m_vectDomains[lev].domainBox())
        {
          CH_assert(foundLevel == false);
          foundLevel = true;
          if (levelGrids == m_vectBeta[lev]->getBoxes())
            {
              // do this fab-by-fab to avoid comminication and to grab 
              // any relevant ghost cells
              const LevelData<FArrayBox>& srcBeta = *m_vectBeta[lev];
              DataIterator dit = levelGrids.dataIterator();
              for (dit.begin(); dit.ok(); ++dit)
                {
                  a_C[dit].copy(srcBeta[dit]);
                }
            }
          else
            {
              // do this with a copyTo
              m_vectBeta[lev]->copyTo(a_C);
            }
        }

    }
  if (!foundLevel)
    {
      MayDay::Error("CopyBasalFriction -- no matching level found");
    }

}
