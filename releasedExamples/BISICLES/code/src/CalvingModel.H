
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CALVINGMODEL_H_
#define _CALVINGMODEL_H_

#include "LevelSigmaCS.H"
#include "ParmParse.H"
#include "SurfaceFlux.H"
#include "NamespaceHeader.H"

class AmrIce;

///virtual base for calving models
class CalvingModel
{

public:

  enum Stage { Initialization, PostVelocitySolve, PostThicknessAdvection, PostRegrid, MAX_STAGE};

  static void updateCalvedIce(const Real& a_thck, const Real a_prevThck, const int a_mask, Real& a_added, Real& a_calved, Real& a_removed);


  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness,
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
 			      LevelData<FArrayBox>& a_removedIce, 
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage) = 0;


  virtual ~CalvingModel(){};

  static CalvingModel* parseCalvingModel(const char* a_prefix);

  /// factory method
  virtual CalvingModel* new_CalvingModel() = 0;

  /// frontal ablation rate.
  virtual void getCalvingRate(LevelData<FArrayBox>& a_calvingRate, 
			      const AmrIce& a_amrIce,int a_level);

  /// water depth for crevasse calving.
  virtual void getWaterDepth(LevelData<FArrayBox>& a_waterDepth, 
			      const AmrIce& a_amrIce,int a_level);
 
};

///the simplest calving model : don't do anything
class NoCalvingModel : public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage)
  {

  }

 virtual CalvingModel* new_CalvingModel()
  {
    return new NoCalvingModel( *this);
  }
  

};


// A Deglaciation-only calving model which 
// i. Makes calving irreversible : if a cell mask is open sea / land, set thickness
//    there to zero. No need to worry about the velocity fields in cells
//    which were previously calved.
// ii. Sets any floating ice below a given thickness and in a given depth of
//     ocean to zero thickness
// iii. prevents thickness from dropping below a given thickness
class DeglaciationCalvingModelA : public CalvingModel
{

public:


  //alter the thickness field at the end of a time step
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);



  DeglaciationCalvingModelA(const Real& a_calvingThickness,
			    const Real& a_calvingOceanDepth,
			    const Real& a_minThickness, 
			    const Real& a_startTime, 
			    const Real& a_endTime)
    : m_calvingThickness(a_calvingThickness),
      m_calvingOceanDepth(a_calvingOceanDepth),
      m_minThickness(a_minThickness),
      m_startTime(a_startTime),
      m_endTime(a_endTime)
  {
    ;
  }

  virtual CalvingModel* new_CalvingModel()
  {
    return new DeglaciationCalvingModelA( *this);
  }
  

  
protected:
  Real m_calvingThickness, m_calvingOceanDepth, m_minThickness, m_startTime, m_endTime;
  
private:
  DeglaciationCalvingModelA();
};

/// Ensure open ocean on one or more domain edges
class DomainEdgeCalvingModel :  public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);

  DomainEdgeCalvingModel(const Vector<int>& a_frontLo,
			 const Vector<int>& a_frontHi,
			 bool a_preserveSea, 
			 bool a_preserveLand)
    
    :m_frontLo(a_frontLo),m_frontHi(a_frontHi),
     m_preserveSea(a_preserveSea),m_preserveLand(a_preserveLand)
  {
    ;
  }


  virtual CalvingModel* new_CalvingModel()
  {
    return new  DomainEdgeCalvingModel( *this);
  }
  

  
protected:
  Vector<int> m_frontLo;
  Vector<int> m_frontHi;
  bool m_preserveSea, m_preserveLand;
private:
  DomainEdgeCalvingModel();
};


///Calve ice based on proximity to the grounding line
class ProximityCalvingModel : public CalvingModel
{

  Real m_proximity, m_velocity,  m_startTime, m_endTime;

public:
  ProximityCalvingModel(const Real& a_proximity, 
			const Real& a_velocity,
			const Real& a_startTime, 
			const Real& a_endTime)
    : m_proximity(a_proximity), 
      m_velocity(a_velocity),
      m_startTime(a_startTime), 
      m_endTime(a_endTime)
  {
  }
 
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);
 


  virtual CalvingModel* new_CalvingModel()
  {
    return new ProximityCalvingModel( *this);
  }
  
};



// A Deglacation-only calving model which 
// i. Makes calving irreversible : if a cell mask is open sea / land, set thickness
//    there to zero. No need to worry about the velocity fields in cells
//    which were previously calved.
// ii. Sets any floating ice below a given thickness to the min thickness in (iii) (this is where it differs from ModelA, which doesn't actually use the CalvingThickness)
// iii. prevents thickness from dropping below a given thickness
class DeglaciationCalvingModelB : public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);


  DeglaciationCalvingModelB(const Real& a_calvingThickness,
			    const Real& a_calvingOceanDepth,
			    const Real& a_minThickness, 
			    const Real& a_startTime, 
			    const Real& a_endTime)
    : m_calvingThickness(a_calvingThickness),
      m_calvingOceanDepth(a_calvingOceanDepth),
      m_minThickness(a_minThickness),
      m_startTime(a_startTime),
      m_endTime(a_endTime)
  {
    ;
  }
  virtual CalvingModel* new_CalvingModel()
  {
    return new DeglaciationCalvingModelB (*this); 
  }

protected:
  Real m_calvingThickness, m_calvingOceanDepth, m_minThickness, m_startTime, m_endTime;
  
private:
  DeglaciationCalvingModelB();
};



/// A  thickness-based calving model which should permit advance and retreat
/** (i). Sets any floating ice below a given thickness and in a given depth of
         ocean to zero thickness
   (ii). uses real-valued ice fraction in AmrIce to account for partial cells
         using the scheme described in the "Notes on MISOMIP Thickness-driven
         Calving Law" in the doc directory
*/
class ThicknessCalvingModel : public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);


  ThicknessCalvingModel(const Real& a_calvingThickness,
                        const Real& a_calvingOceanDepth,
                        const Real& a_minThickness, 
                        const Real& a_startTime, 
                        const Real& a_endTime,
			const bool& a_factorMuCoef)
    : m_calvingThickness(a_calvingThickness),
      m_calvingOceanDepth(a_calvingOceanDepth),
      m_minThickness(a_minThickness),
      m_startTime(a_startTime),
      m_endTime(a_endTime),
      m_factorMuCoef(a_factorMuCoef)
  {
    ;
  }
  virtual CalvingModel* new_CalvingModel()
  {
    return new ThicknessCalvingModel( *this);
  }
protected:
  Real m_calvingThickness, m_calvingOceanDepth, m_minThickness, m_startTime, m_endTime;
  bool m_factorMuCoef; // test muCoef*thickness rather than thickness?
private:
  ThicknessCalvingModel();
};


/// A calving model which enforces a maximum extent of the ice shelf
class MaximumExtentCalvingModel : public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);
  
  MaximumExtentCalvingModel(const RealVect& a_highLoc,
                            const RealVect& a_lowLoc,
                            const Real& a_startTime, 
                            const Real& a_endTime)
    : m_startTime(a_startTime),
      m_endTime(a_endTime),
      m_highLoc(a_highLoc),
      m_lowLoc(a_lowLoc)      
  {
    ;
  }

  virtual CalvingModel* new_CalvingModel()
  {
    return new MaximumExtentCalvingModel ( *this);
  }

  
protected:
  Real m_startTime, m_endTime;
  RealVect m_highLoc, m_lowLoc;
  
private:
  MaximumExtentCalvingModel();
};



/// Calving model which lets us compose multiple calving models
/** Allows serial application of  multiple calving models (for example, FixedFrontCalvingModel plus ThicknessCalvingModel). The ordering here is likely important, so always perform operations starting with calving model 0, then 1, etc.
 */
class CompositeCalvingModel : public CalvingModel
{

public:

  CompositeCalvingModel(Vector<CalvingModel*> a_vectModels)
    : m_vectModels(a_vectModels)
  {
    m_vectModels.resize(a_vectModels.size());
    for (int i =0; i < m_vectModels.size(); i++)
      {
	m_vectModels[i] = a_vectModels[i]->new_CalvingModel();
      }   
  }
   
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);


  virtual ~CompositeCalvingModel();
  
  // this is public for simplicity
  Vector<CalvingModel*> m_vectModels;

  virtual CalvingModel* new_CalvingModel()
  {
    return new CompositeCalvingModel (*this); 
  }

  CompositeCalvingModel(const CompositeCalvingModel& a)
  {
    m_vectModels.resize(a.m_vectModels.size());
    for (int i =0; i < a.m_vectModels.size(); i++)
      {
	m_vectModels[i] = a.m_vectModels[i]->new_CalvingModel();
      }
  }
  
  
};


/// calving model which removes all floating ice. 
class FlotationCalvingModel : public CalvingModel
{

public:

  
  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);

  FlotationCalvingModel (const Vector<int>& a_frontLo,
			 const Vector<int>& a_frontHi,
			 bool a_preserveSea,
			 bool a_preserveLand) 
    : m_domainEdgeCalvingModel(a_frontLo,a_frontHi, 
			       a_preserveSea, a_preserveLand)
  {
    
  }

  virtual CalvingModel* new_CalvingModel()
  {
    return new  FlotationCalvingModel (*this); 
  }

  
private:
  DomainEdgeCalvingModel m_domainEdgeCalvingModel ;

  
};

/// A calving model which implements the Pollard and DeConto cliff-collapse mechanism using the ice mask
/** 
   uses real-valued ice mask in AmrIce to account for partial cells
   using the scheme described in the "Notes on MISOMIP Thickness-driven
   Calving Law" in the doc directory -- recession due to cliff collapse is implemented by removing 
   cell-averaged thickness along with cell area
*/
class CliffCollapseCalvingModel : public CalvingModel
{

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);


  CliffCollapseCalvingModel(const Real& a_maxCliffThickness,
			    const Real& a_recessionRate,
			    const Real& a_startTime, 
			    const Real& a_endTime)
    : m_maxCliffThickness(a_maxCliffThickness),
      m_recessionRate(a_recessionRate),
      m_startTime(a_startTime),
      m_endTime(a_endTime)
  {
    ;
  }

  virtual CalvingModel* new_CalvingModel()
  {
    return new CliffCollapseCalvingModel (*this); 
  }

  
protected:
  Real m_maxCliffThickness, m_recessionRate, m_startTime, m_endTime;

private:
  CliffCollapseCalvingModel();
  
};


/// A Calving model that provides a simple rate
/**
   The rate is cell-centred (which fits with the 
   current implementation in AmrIce::advectIceFrac),
   and is specified by a SurfaceFlux.
 */
class VariableRateCalvingModel : public CalvingModel
{

  Real m_startTime, m_endTime;
  DomainEdgeCalvingModel* m_domainEdgeCalvingModel;
  SurfaceFlux* m_calvingRate;

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
 			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);
  
  VariableRateCalvingModel (ParmParse& a_pp);

  virtual ~VariableRateCalvingModel();
  
  virtual CalvingModel* new_CalvingModel();
  
  virtual void getCalvingRate(LevelData<FArrayBox>& a_calvingRate, 
			      const AmrIce& a_amrIce,int a_level);

};


/// A Calving model that provides a rate proportional to speed: the constant of proportionality can vary in x,y,t
/**
   The rate is cell-centred (which fits with the  current implementation in AmrIce::advectIceFrac),
   The constant is specified by a SurfaceData.
 */
class RateProportionalToSpeedCalvingModel : public CalvingModel
{

  Real m_startTime, m_endTime;
  DomainEdgeCalvingModel* m_domainEdgeCalvingModel;
  SurfaceFlux* m_proportion;

public:

  virtual void applyCriterion(LevelData<FArrayBox>& a_thickness, 
 			      LevelData<FArrayBox>& a_calvedIce,
			      LevelData<FArrayBox>& a_addedIce,
			      LevelData<FArrayBox>& a_removedIce,
			      LevelData<FArrayBox>& a_iceFrac, 
			      const AmrIce& a_amrIce,
			      int a_level,
			      Stage a_stage);
  
  RateProportionalToSpeedCalvingModel (ParmParse& a_pp);

  virtual ~RateProportionalToSpeedCalvingModel();
  
  virtual CalvingModel* new_CalvingModel();
  
  virtual void getCalvingRate(LevelData<FArrayBox>& a_calvingRate, 
			      const AmrIce& a_amrIce,int a_level);

};


#include "NamespaceFooter.H"

#endif
