#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CREVASSECALVINGMODEL_H_
#define _CREVASSECALVINGMODEL_H_

#include "LevelSigmaCS.H"
#include "ParmParse.H"
#include "CalvingModel.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

/// Calving models based on crevase depth computed from viscous stress 
/// CrevasseCalvingModel implements the methods common to all of these, 
/// and requires its subclasses to compute a remnnat given the stress.
/// The remnant could be the remaining thickness of ice 
///  remnant = (thickness - crevasse depth)/thickness
/// but could be any quantity where remnant = 0 indiciates full fracture and
/// remnant > 0 is a partial fracture.
class CrevasseCalvingModel : public CalvingModel
{
public:

  enum StressMeasure {FirstPrincipalStress,AlongFlowNormalStress,Trace,MAX_STRESS_MEASURE} ;

  CrevasseCalvingModel(ParmParse& a_pp);

  virtual ~CrevasseCalvingModel();

  virtual void computeRemnant(LevelData<FArrayBox>& a_remnant,
			      const LevelData<FArrayBox>& a_stress,
			      const LevelData<FArrayBox>& a_thck,
			      const LevelData<FArrayBox>& a_usrf,
			      const LevelData<FArrayBox>& a_hab, 
			      const LevelData<FArrayBox>& a_waterDepth, 
			      const LevelSigmaCS& a_coords) = 0;
			      

  // CalvingModel implementation
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness,
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,  
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);

  virtual void getWaterDepth(LevelData<FArrayBox>& a_waterDepth, 
			      const AmrIce& a_amrIce,int a_level);

  CrevasseCalvingModel(const CrevasseCalvingModel& a) 
  {
    m_waterDepth = a.m_waterDepth->new_surfaceFlux();
    m_includeBasalCrevasses = a.m_includeBasalCrevasses;
    m_calvingZoneLength = a.m_calvingZoneLength;
    m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*a.m_domainEdgeCalvingModel);
    m_stressMeasure = a.m_stressMeasure;
  }

protected:
  //compute a scalar stress measure from the viscous stress tensor
  void computeStressMeasure(LevelData<FArrayBox>& a_stressMeasure,
			    const AmrIce& a_amrIce,int a_level);


  SurfaceFlux* m_waterDepth;
  bool m_includeBasalCrevasses;
  Real m_calvingZoneLength;
  DomainEdgeCalvingModel* m_domainEdgeCalvingModel;
  StressMeasure m_stressMeasure; 

};

/// see (for 1D) Benn (2007), Earth-Science Reviews vol 82 p 143 doi:10.1016/j.earscirev.2007.02.002
/// Nick et al (2011) J. Glaciol vol 56 p 781 doi:10.3189/002214310794457344
/// author : Andrew Taylor, University of Bristol
class BennCalvingModel : public CrevasseCalvingModel
{

public:  
  
  BennCalvingModel(ParmParse& a_pp);

  ~BennCalvingModel();

  virtual void computeRemnant(LevelData<FArrayBox>& a_remnant,
			      const LevelData<FArrayBox>& a_stress,
			      const LevelData<FArrayBox>& a_thck,
			      const LevelData<FArrayBox>& a_usrf,
			      const LevelData<FArrayBox>& a_hab,
			      const LevelData<FArrayBox>& a_waterDepth,
			      const LevelSigmaCS& a_coords);

  BennCalvingModel(const BennCalvingModel& a) : CrevasseCalvingModel(a)
  {
  }


  /// factory method
  virtual CalvingModel* new_CalvingModel()
  {
    return new BennCalvingModel(*this);
  }
};

/// Calving model based on crevasse opening
/// see Van der Veen (1998) , Cold Regions Sci. Tech. vol 27 p 231
class VdVCalvingModel : public CrevasseCalvingModel
{

public:  
  
  VdVCalvingModel(ParmParse& a_pp);

  ~VdVCalvingModel();

  virtual void computeRemnant(LevelData<FArrayBox>& a_remnant,
			      const LevelData<FArrayBox>& a_stress,
			      const LevelData<FArrayBox>& a_thck,
			      const LevelData<FArrayBox>& a_usrf,
			      const LevelData<FArrayBox>& a_hab, 
			      const LevelData<FArrayBox>& a_waterDepth,
			      const LevelSigmaCS& a_ccords);

  VdVCalvingModel(const VdVCalvingModel& a) : CrevasseCalvingModel(a)
  {
  }

  /// factory method
  virtual CalvingModel* new_CalvingModel()
  {
    return new VdVCalvingModel(*this);
  }
};
  

#include "NamespaceFooter.H"

#endif
