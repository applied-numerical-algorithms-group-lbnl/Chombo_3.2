#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Virtual base class encapsulating coefficient of effective viscosity 
// for ice-sheet model. NB : this is an additional coefficient a, so
// that the viscous tensor is (a * mu * rate-of-strain)

#ifndef _MUCOEFFICIENT_H_
#define _MUCOEFFICIENT_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"


class MuCoefficient
{
public:

  ///assemble a MuCoefficient  object from ParmParse input, return pointer
  static MuCoefficient* parseMuCoefficient(const char* a_prefix);

  /// constructor
  MuCoefficient() {};

  /// destructor
  virtual ~MuCoefficient() {};

  /// factory method
  /** return a pointer to a new MuCoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const = 0;
 /// define cell and face-centered coefficient
  virtual void setMuCoefficient(LevelData<FArrayBox>& a_muCoef,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt) = 0;

 
private:
  // Disallowed for all the usual reasons
  void operator=(const MuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MuCoefficient(const MuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


// simple version which sets mu coefficient to one
class UnitMuCoefficient : public MuCoefficient
{
public:
  /// constructor
  UnitMuCoefficient() {};

  /// destructor
  virtual ~UnitMuCoefficient() {};

  /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  /// define cell and face-centered coefficient
  /** time and dt are included in case this is time-dependent.
  */
  virtual void setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);


private:
  // Disallowed for all the usual reasons
  void operator=(const UnitMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  UnitMuCoefficient(const UnitMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


//class which contains two muCoefficients (x and y)
//and two scalars (a and b) and provides the
//coefficient a x + b y
class AxbyMuCoefficient : public MuCoefficient
{

  MuCoefficient* m_x;
  MuCoefficient* m_y;
  Real m_a,m_b;

public:
  /// constructor
  AxbyMuCoefficient(const Real& a_a, MuCoefficient* a_x, 
		    const Real& a_b, MuCoefficient* a_y);

  /// destructor
  virtual ~AxbyMuCoefficient();

  /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  /// define cell and face-centered coefficient
  /** time and dt are included in case this is time-dependent.
  */
  virtual void setMuCoefficient(LevelData<FArrayBox>& a_cellMuCoef,
                                LevelSigmaCS& a_coordSys,
                                Real a_time,
                                Real a_dt);


private:
  // Disallowed for all the usual reasons
  void operator=(const AxbyMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AxbyMuCoefficient(const AxbyMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

/// MuCoeffcient that computes \f$ \phi \f$ from data on a uniform grid
/**
   
   Run time configuration is carried out in BasalFriction::parse(A, ... )
   If the option `A.type = LevelData` is found, then the 
   expected form is:
   
   - `inputLevelData.frictionFile = <string>` 
   - `inputLevelData.frictionName = <string>` 
   
   for a single file or
   
   - `inputLevelData.frictionFileFormat = <string>` 
   - `inputLevelData.frictionName = <string>` 
   - `inputLevelData.frictionFileStartTime = <int>`
   - `inputLevelData.frictionFileTimeStep = <int>`
   
   for a sequence of files
   
*/
class LevelDataMuCoefficient : public MuCoefficient
{

public:
  //constructor
  LevelDataMuCoefficient
  (RefCountedPtr<LevelData<FArrayBox> > a_muCoef,
   const RealVect& a_dx)
    :m_timeFileMap(NULL),m_startTime(-1.2345678e+300),m_endTime(1.2345678e+300),
    m_muCoef(a_muCoef),m_endMuCoef(a_muCoef),m_dx(a_dx),m_verbose(true),m_name("")
  {
  }

   /// construct from a sequence of files
  /**
     @param a_timeFileMap maps model times to filenames
     @param a_name  name of the field in all hdf5 files
     
   */
  LevelDataMuCoefficient(RefCountedPtr<std::map<Real, std::string> >& a_timeFileMap, 
		       const std::string& a_name)
    :m_timeFileMap(a_timeFileMap),m_startTime(1.2345678e+299),m_endTime(-1.2345678e+299),
     m_muCoef(new LevelData<FArrayBox>()),m_endMuCoef(new LevelData<FArrayBox>())
    ,m_verbose(true),m_name(a_name)
  {
    
  }
  
    /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  virtual void setMuCoefficient
  (LevelData<FArrayBox>& a_muCoef,
   LevelSigmaCS& a_coordSys,
   Real a_time,
   Real a_dt) ;

protected:
  RefCountedPtr<std::map<Real,std::string> > m_timeFileMap;
  Real m_startTime, m_endTime;
  RefCountedPtr< LevelData<FArrayBox> > m_muCoef; // P at t = m_startTime
  RefCountedPtr< LevelData<FArrayBox> > m_endMuCoef; // P at t = m_endTime
  RealVect m_dx;
  bool m_verbose;
  std::string m_name;

private:
  void operator=(const LevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  LevelDataMuCoefficient(const LevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

class MultiLevelDataMuCoefficient : public MuCoefficient
{

public:
  //constructor
  MultiLevelDataMuCoefficient
  (const Vector< RefCountedPtr<LevelData<FArrayBox> > >& a_muCoef,
   const RealVect& a_dxCrse, const Vector<int> & a_ratio)
  :m_muCoef(a_muCoef),m_dxCrse(a_dxCrse),m_ratio(a_ratio),m_verbose(true)
  {
  }
    /// factory method
  /** return a pointer to a new Mucoefficient object
   */
  virtual MuCoefficient* new_muCoefficient() const;

  virtual void setMuCoefficient
  (LevelData<FArrayBox>& a_muCoef,
   RealVect a_dx,
   Real a_time,
   Real a_dt) ;

  virtual void setMuCoefficient
  (LevelData<FArrayBox>& a_muCoef,
   LevelSigmaCS& a_coordSys,
   Real a_time,
   Real a_dt) ;

protected:
  Vector< RefCountedPtr< LevelData<FArrayBox> > > m_muCoef;
  RealVect m_dxCrse;
  Vector<int> m_ratio;
  bool m_verbose;

private:
  void operator=(const MultiLevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

  MultiLevelDataMuCoefficient(const MultiLevelDataMuCoefficient& a_input)
  {
    MayDay::Error("invalid operator");
  }

};



#include "NamespaceFooter.H"


#endif // multiple-include preventer
