
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _MASKEDCALVINGMODEL_H_
#define _MASKEDCALVINGMODEL_H_

#include "CalvingModel.H"
#include "LevelSigmaCS.H"
#include "ParmParse.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

/// Calving model which uses a masked input to enforce shelf collapse
class MaskedCalvingModel : public CalvingModel
{

public:

// default constructor
MaskedCalvingModel();

// constructor with calving mask in the form of a surfaceFlux
MaskedCalvingModel(SurfaceFlux* a_calvingMaskPtr,
                   Real m_minThickness);

// destructor
virtual ~MaskedCalvingModel();

// define
void define(SurfaceFlux* a_calvingMaskPtr,
            Real m_minThickness);


   //alter the thickness field at the end of a time step  
   virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
                               LevelData<FArrayBox>& a_calvedIce,
                               LevelData<FArrayBox>& a_addedIce,
                               LevelData<FArrayBox>& a_removedIce,
			       LevelData<FArrayBox>& a_iceFrac, 
			       const AmrIce& a_amrIce,
			       int a_level,
			       Stage a_stage);


virtual CalvingModel* new_CalvingModel();

/// make this public
Real m_calvingVal;

protected:

SurfaceFlux* m_calvingMask;

Real m_minThickness;

private: 

};



#include "NamespaceFooter.H"
#endif
