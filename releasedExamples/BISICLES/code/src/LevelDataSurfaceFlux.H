#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes determined  
// by copying, coarsening, or interpolating a LevelData<FArrayBox>
// covering an entire domain. Objects can be defined either
// by specifying a RefCountedPtr<LevelData<FArrayBox> > , or by specifying
// a std::map<Real,string> mapping time ti to an hdf5 file f. In the
// latter case, the flux at time t is found by linear interploation
// between time max(ti <= t) and min(ti > t) 

// \todo replace the std::map<Real,string> mechanism with
// a suitable abstraction

#ifndef _LEVELDATASURFACEFLUX_H_
#define _LEVELDATASURFACEFLUX_H_

#include "SurfaceFlux.H"
#include <map>
#include "NamespaceHeader.H"

class LevelDataSurfaceFlux : public SurfaceFlux
{
public:

   LevelDataSurfaceFlux(RefCountedPtr<LevelData<FArrayBox> > a_flux,
			const RealVect& a_dx, const Real& a_defaultValue=0.0)
     :m_timeFileMap(NULL), m_startTime(-1.2345678e+300), m_endTime(1.2345678e+300),
      m_defaultValue(a_defaultValue),
      m_startFlux(a_flux), m_endFlux(a_flux),m_dx(a_dx),m_verbose(true),m_name(""),
      m_linearInterp(true)
     
   {
   }
  
  LevelDataSurfaceFlux(RefCountedPtr<std::map<Real,
                       std::string> >& a_timeFileMap, 
		       const std::string& a_name,
		       const bool& a_linearInterp,
		       const Real& a_defaultValue=0.0)
    :m_timeFileMap(a_timeFileMap), m_startTime(1.2345678e+300), m_endTime(-1.2345678e+300),
     m_defaultValue(a_defaultValue),
     m_startFlux(new LevelData<FArrayBox>()),m_endFlux(new LevelData<FArrayBox>()),
     m_verbose(true),m_name(a_name),m_linearInterp(a_linearInterp)
  {
  }

 


  /// destructor
  virtual ~LevelDataSurfaceFlux() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux() 
  {
    LevelDataSurfaceFlux* ptr;
    if (m_timeFileMap != NULL)
      {
	ptr = new LevelDataSurfaceFlux(m_timeFileMap,m_name,m_linearInterp, m_defaultValue);
      }
    else
      {
     	ptr = new LevelDataSurfaceFlux(m_startFlux,m_dx,m_defaultValue );
      }
    
    return static_cast<SurfaceFlux*>(ptr);
  }

  /// define source term for thickness evolution and place it in flux
  /** 
  */

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  

protected:

  RefCountedPtr<std::map<Real,std::string> > m_timeFileMap;
  Real m_startTime, m_endTime, m_defaultValue;
  RefCountedPtr< LevelData<FArrayBox> > m_startFlux;
  RefCountedPtr< LevelData<FArrayBox> > m_endFlux;
  RealVect m_dx;
  bool m_verbose;
  std::string m_name;
  bool m_linearInterp;
  
private:
  // Disallowed for all the usual reasons
  void operator=(const LevelDataSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
  
  // public interfaces disallowed for all the usual reasons
  LevelDataSurfaceFlux(const LevelDataSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

};


class MultiLevelDataSurfaceFlux : public SurfaceFlux
{
public:

 
  MultiLevelDataSurfaceFlux(RefCountedPtr<std::map<Real,std::string> >& a_timeFileMap, 
			    const std::string& a_name,
			    const bool& a_linearInterp,
                            const Real& a_defaultValue=0.0)
    :m_timeFileMap(a_timeFileMap), m_startTime(1.2345678e+300), m_endTime(-1.2345678e+300), m_defaultValue(a_defaultValue),
     m_verbose(true),m_name(a_name), m_linearInterp(a_linearInterp)
  {
  }

 


  /// destructor
  virtual ~MultiLevelDataSurfaceFlux() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux()
  {
    MultiLevelDataSurfaceFlux* ptr;
    if (m_timeFileMap != NULL)
      {
	ptr = new MultiLevelDataSurfaceFlux(m_timeFileMap,m_name,m_linearInterp,
                                            m_defaultValue);
      }
    // else
    //   {
    // 	ptr = new MultiLevelDataSurfaceFlux(m_startFlux,m_dx);
    //   }

    return static_cast<SurfaceFlux*>(ptr);
  }

  /// define source term for thickness evolution and place it in flux
  /** 
  */

  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				    const AmrIceBase& a_amrIce, 
				    int a_level, Real a_dt);
  

protected:

  RefCountedPtr<std::map<Real,std::string> > m_timeFileMap;
  Real m_startTime, m_endTime, m_defaultValue;
  Vector<RefCountedPtr< LevelData<FArrayBox> > > m_startFlux;
  Vector<RefCountedPtr< LevelData<FArrayBox> > > m_endFlux;
  RealVect m_dxCrse;
  Vector<int> m_ratio;
  bool m_verbose;
  std::string m_name;
  bool m_linearInterp;

private:
  // Disallowed for all the usual reasons
  void operator=(const MultiLevelDataSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }
  
  // public interfaces disallowed for all the usual reasons
  MultiLevelDataSurfaceFlux(const MultiLevelDataSurfaceFlux& a_input)
  {
    MayDay::Error("invalid operator");
  }

};


#include "NamespaceFooter.H"

#endif
