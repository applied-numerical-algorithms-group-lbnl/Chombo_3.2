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
#include "FortranInterfaceBasalFriction.H"
FortranInterfaceBasalFriction::FortranInterfaceBasalFriction()
{
  m_verbose = false;
}

BasalFriction* FortranInterfaceBasalFriction::new_basalFriction() const
{
  
  FortranInterfaceBasalFriction* ptr = new FortranInterfaceBasalFriction();
  ptr->m_ghost = m_ghost;
  ptr->m_dx = m_dx;
  ptr->m_fab.define(m_fab.box(),m_fab.nComp());
  ptr->m_fab.copy(m_fab);
  return static_cast<BasalFriction*>(ptr);
}

void
FortranInterfaceBasalFriction::setReferenceFAB(Real* a_data_ptr,
					       const int* a_dimInfo,
					       const RealVect& a_dx,
					       const IntVect& a_ghost,
					       const bool a_nodal)
{
  m_ghost = a_ghost;
  m_dx = a_dx;
  
  IntVect hiVect(D_DECL(a_dimInfo[2]-1,a_dimInfo[3]-1, a_dimInfo[1]-1));
 
  if (a_nodal)
    {
      Box nodeBox (IntVect::Zero, hiVect);
      nodeBox.shift(-a_ghost);
      nodeBox.shift(IntVect::Unit);
      FArrayBox nodeFAB;
      nodeFAB.define(nodeBox,1,a_data_ptr);
      Box fabBox(IntVect::Zero, hiVect - IntVect::Unit);
      m_fab.define(fabBox, 1);
      FORT_NODETOCELL(CHF_CONST_FRA1(nodeFAB,0),
		      CHF_FRA1(m_fab,0),
		      CHF_BOX(fabBox));
    }
  else 
    {
       Box fabBox(IntVect::Zero, hiVect);
       fabBox.shift(-a_ghost);
      m_fab.define(fabBox, 1, a_data_ptr);
    } 
}

void
FortranInterfaceBasalFriction::setBasalFriction(LevelData<FArrayBox>& a_C,
						LevelSigmaCS& a_coordSys,
						Real a_time,
						Real a_dt)
{
  FillFromReference(a_C, m_fab, a_coordSys.dx(), m_dx, m_ghost, m_verbose);
}
