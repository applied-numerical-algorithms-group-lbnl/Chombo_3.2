#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MaskedCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "NamespaceHeader.H"

// default constructor
MaskedCalvingModel::MaskedCalvingModel()
  :m_calvingMask(NULL), m_minThickness(-10.0), m_calvingVal(0.1)
{
}

// constructor with calving mask in the form of a surfaceFlux
MaskedCalvingModel::MaskedCalvingModel(SurfaceFlux* a_calvingMaskPtr,
                                       Real a_minThickness)
{
  define(a_calvingMaskPtr, a_minThickness);
}

MaskedCalvingModel::~MaskedCalvingModel()
{
  if (m_calvingMask != NULL)
    {
      delete m_calvingMask;
      m_calvingMask = NULL;
    }
}

// define
void 
MaskedCalvingModel::define(SurfaceFlux* a_calvingMaskPtr,
                           Real  a_minThickness)
{
  m_minThickness = a_minThickness;
  m_calvingMask = a_calvingMaskPtr->new_surfaceFlux();
  // set a default calving value to be 0.1
  m_calvingVal = 0.1;
}


//alter the thickness field at the end of a time step  
void 
MaskedCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness, 
                                   LevelData<FArrayBox>& a_calvedIce,
                                   LevelData<FArrayBox>& a_addedIce,
                                   LevelData<FArrayBox>& a_removedIce,
                                   LevelData<FArrayBox>& a_iceFrac, 
                                   const AmrIce& a_amrIce,
                                   int a_level,
                                   Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  // for now, don't calve in ghost cells
  IntVect ghostVect = IntVect::Zero;
  const DisjointBoxLayout& grids = levelCoords.grids();
  LevelData<FArrayBox> calvingMask(grids, 1, ghostVect);
  Real dt = 1.0;
  m_calvingMask->evaluate(calvingMask, a_amrIce, a_level, dt);


  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& thisCalvingMask = calvingMask[dit];
      Box b = grids[dit];

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if ((mask(iv) == FLOATINGMASKVAL) && (thisCalvingMask(iv,0) >= m_calvingVal))
	    {
	      thck(iv) = m_minThickness;
	    }
	      
	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }

}



CalvingModel*
MaskedCalvingModel::new_CalvingModel()
{
  MaskedCalvingModel* newCM = new MaskedCalvingModel(m_calvingMask,
                                                     m_minThickness);

  return static_cast<CalvingModel*>(newCM);
}

#include "NamespaceFooter.H"
