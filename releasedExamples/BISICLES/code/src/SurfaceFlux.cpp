#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SurfaceFlux.H"
#include "ComplexSurfaceFlux.H"
#include "LevelDataSurfaceFlux.H"
#include "ISMIP6OceanForcing.H"
#include "GroundingLineLocalizedFlux.H"
#include "HotspotFlux.H"
#include <map>
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "IceConstants.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "BisiclesF_F.H"
#include "ParmParse.H"
#include "AmrIceBase.H"
#include "FortranInterfaceIBC.H"
#include "FillFromReference.H"

#include "NamespaceHeader.H"

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */

SurfaceFlux* 
zeroFlux::new_surfaceFlux()
{
  zeroFlux* newPtr = new zeroFlux;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
zeroFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
			       const AmrIceBase& a_amrIce, 
			       int a_level, Real a_dt)
{
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}


constantFlux::constantFlux() : m_isValSet(false)
{
}

SurfaceFlux* 
constantFlux::new_surfaceFlux()
{
  constantFlux* newPtr = new constantFlux;
  newPtr->m_fluxVal = m_fluxVal;
  newPtr->m_isValSet = m_isValSet;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
constantFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				   const AmrIceBase& a_amrIce, 
				   int a_level, Real a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(m_fluxVal);
    }
}


///
void
constantFlux::setFluxVal(const Real& a_fluxVal) 
{
  m_fluxVal = a_fluxVal; 
  // input value is in meters/year divide by secondsperyear 
  // to get flux in meters/second
  // slc: switch to flux in m/a
  //m_fluxVal/= secondsperyear;
  
  m_isValSet = true;
}


SurfaceFlux* SurfaceFlux::parse(const char* a_prefix)
{
  
  SurfaceFlux* ptr = NULL;
  std::string type = "";
  
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "zeroFlux")
    {
      ptr = new zeroFlux;
    }
  else if (type == "constantFlux")
    {
      constantFlux* constFluxPtr = new constantFlux;
      Real fluxVal;
      pp.get("flux_value", fluxVal);
      constFluxPtr->setFluxVal(fluxVal);
      ptr = static_cast<SurfaceFlux*>(constFluxPtr);
    }
  else if (type == "hotspotFlux")
    {
      HotspotFlux* hotspotFluxPtr = new HotspotFlux;
      Real fluxVal;
      pp.get("flux_value", fluxVal);
      hotspotFluxPtr->setFluxVal(fluxVal);
      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("radius",vect,0,SpaceDim);
      RealVect radius(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("center",vect,0,SpaceDim);
      RealVect center(D_DECL(vect[0], vect[1],vect[2]));

      hotspotFluxPtr->setSpotLoc(radius, center);
      
      Real startTime = -1.2345e300;
      Real stopTime = 1.2345e300;
      pp.query("start_time", startTime);
      pp.query("stop_time", stopTime);
      hotspotFluxPtr->setSpotTimes(startTime, stopTime);
      
      ptr = static_cast<SurfaceFlux*>(hotspotFluxPtr);
    }
  else if (type == "LevelData")
    {
      std::string fileFormat;
      pp.get("fileFormat",fileFormat);
      int n;
      pp.get("n",n);
      int offset = 0;
      pp.query("offset",offset);

      Real startTime = 0.0, timeStep = 1.0;
      pp.query("startTime", startTime);
      pp.query("timeStep", timeStep);
      std::string name = "flux";
      pp.query("name", name);
      bool linearInterp = true;
      pp.query("linearInterp", linearInterp);

      Real defaultValue = 0.0;
      pp.query("defaultValue", defaultValue);
      
      RefCountedPtr<std::map<Real,std::string> > tf
	(new std::map<Real,std::string>);
      
      for (int i =0; i < n; i++)
	{
	  char* file = new char[fileFormat.length()+32];
	  sprintf(file, fileFormat.c_str(),i + offset);
	  tf->insert(make_pair(startTime + Real(i)*timeStep, file));
	  delete[] file;
	}
      
      LevelDataSurfaceFlux* ldptr = new LevelDataSurfaceFlux(tf,name,linearInterp, defaultValue);
      ptr = static_cast<SurfaceFlux*>(ldptr);
    }
  else if (type == "MultiLevelData")
    {
      std::string fileFormat;
      pp.get("fileFormat",fileFormat);
      int n;
      pp.get("n",n);
      int offset = 0;
      pp.query("offset",offset);

      Real startTime = 0.0, timeStep = 1.0;
      pp.query("startTime", startTime);
      pp.query("timeStep", timeStep);
      std::string name = "flux";
      pp.query("name", name);
      bool linearInterp = true;
      pp.query("linearInterp", linearInterp);

      Real defaultValue = 0.0;
      pp.query("defaultValue", defaultValue);
      
      RefCountedPtr<std::map<Real,std::string> > tf
	(new std::map<Real,std::string>);
      
      for (int i =0; i < n; i++)
	{
	  char* file = new char[fileFormat.length()+32];
	  sprintf(file, fileFormat.c_str(),i + offset);
	  tf->insert(make_pair(startTime + Real(i)*timeStep, file));
	  delete[] file;
	}
      
      MultiLevelDataSurfaceFlux* ldptr = new MultiLevelDataSurfaceFlux(tf,name,linearInterp, defaultValue);
      ptr = static_cast<SurfaceFlux*>(ldptr);
    }
  else if (type == "fortran")
    {
      // don't have the context here to actually set values, but
      // we can at least allocate the object here and return it 
      fortranInterfaceFlux* fifptr = new fortranInterfaceFlux;
      ptr = static_cast<SurfaceFlux*>(fifptr);
    }
  else if (type == "piecewiseLinearFlux")
    {
      int n = 1;  
      pp.query("n",n);
      Vector<Real> vabs(n,0.0);
      Vector<Real> vord(n,0.0);
      pp.getarr("abscissae",vabs,0,n);
      pp.getarr("ordinates",vord,0,n);
      
      Real dmin = -1.0;
      pp.query("minWaterDepth",dmin);
      
      PiecewiseLinearFlux* pptr = new PiecewiseLinearFlux(vabs,vord,dmin);
      ptr = static_cast<SurfaceFlux*>(pptr);
    }
  else if (type == "productFlux")
    {
      std::string flux1Prefix(a_prefix);
      flux1Prefix += ".flux1";
      SurfaceFlux* flux1Ptr = parse(flux1Prefix.c_str());
      if (flux1Ptr == NULL)
	{
	  MayDay::Error("undefined flux1 in productFlux");
	}

      std::string flux2Prefix(a_prefix);
      flux2Prefix += ".flux2";
      SurfaceFlux* flux2Ptr = parse(flux2Prefix.c_str());
      if (flux2Ptr == NULL)
	{
	  MayDay::Error("undefined flux2 in productFlux");
	}
 

      ptr = static_cast<SurfaceFlux*>
	(new ProductSurfaceFlux(flux1Ptr->new_surfaceFlux(),
				flux2Ptr->new_surfaceFlux()));

      
      delete flux1Ptr;
      delete flux2Ptr;
    }
  else if (type == "maskedFlux")
    {
      std::string groundedPrefix(a_prefix);
      groundedPrefix += ".grounded";
      SurfaceFlux* groundedPtr = parse(groundedPrefix.c_str());
      if (groundedPtr == NULL)
	{
	  groundedPtr = new zeroFlux;
	}
 
      std::string floatingPrefix(a_prefix);
      floatingPrefix += ".floating";
      SurfaceFlux* floatingPtr = parse(floatingPrefix.c_str());
      if (floatingPtr == NULL)
	{
	  floatingPtr = new zeroFlux;
	}

      std::string openLandPrefix(a_prefix);
      openLandPrefix += ".openLand";
      SurfaceFlux* openLandPtr = parse(openLandPrefix.c_str());
      if (openLandPtr == NULL)
	{
	  openLandPtr = groundedPtr->new_surfaceFlux();
	}

      
      std::string openSeaPrefix(a_prefix);
      openSeaPrefix += ".openSea";
      SurfaceFlux* openSeaPtr = parse(openSeaPrefix.c_str());
      if (openSeaPtr == NULL)
	{
	  openSeaPtr = floatingPtr->new_surfaceFlux();
	}

      ptr = static_cast<SurfaceFlux*>
	(new MaskedFlux(groundedPtr->new_surfaceFlux(),
			floatingPtr->new_surfaceFlux(),
			openSeaPtr->new_surfaceFlux(),
			openLandPtr->new_surfaceFlux()));
      
      delete groundedPtr;
      delete floatingPtr;
      delete openSeaPtr;
      delete openLandPtr;
    }
  else if (type == "boxBoundedFlux")
    {
      Vector<Real> tmp(SpaceDim); 
      pp.getarr("lo",tmp,0,SpaceDim);
      RealVect lo (D_DECL(tmp[0], tmp[1],tmp[2]));
      pp.getarr("hi",tmp,0,SpaceDim);
      RealVect hi (D_DECL(tmp[0], tmp[1],tmp[2]));

      Vector<Real> time(2);
      time[0] = -1.2345678e+300;
      time[1] = 1.2345678e+300;
      pp.queryarr("time",time,0,2);
           
      std::string prefix(a_prefix);
      prefix += ".flux";
      SurfaceFlux* fluxPtr = parse(prefix.c_str());
      CH_assert(fluxPtr != NULL);
      BoxBoundedFlux bbf(lo,hi,time[0],time[1],fluxPtr);
      ptr = static_cast<SurfaceFlux*>(bbf.new_surfaceFlux());

    }
  else if (type == "axbyFlux")
   {
     Real a; 
     pp.get("a",a);
     
     std::string xpre(a_prefix);
     xpre += ".x";
     SurfaceFlux* x = parse(xpre.c_str());
     
     Real b; 
     pp.get("b",b);
     
     std::string ypre(a_prefix);
     ypre += ".y";
     SurfaceFlux* y = parse(ypre.c_str());
    
     AxbyFlux axbyFlux(a,x,b,y);
     ptr = static_cast<SurfaceFlux*>(axbyFlux.new_surfaceFlux());

   }
  else if (type == "compositeFlux")
   {
     
     
     int nElements;
     pp.get("nElements",nElements);
     
     std::string elementPrefix(a_prefix);
     elementPrefix += ".element";

     Vector<SurfaceFlux*> elements(nElements);
     for (int i = 0; i < nElements; i++)
       {
	 std::string prefix(elementPrefix);
	 char s[32];
	 sprintf(s,"%i",i);
	 prefix += s;
	 ParmParse pe(prefix.c_str());
	 elements[i] = parse(prefix.c_str());
	 CH_assert(elements[i] != NULL);
	 
       }
     CompositeFlux compositeFlux(elements);
     ptr = static_cast<SurfaceFlux*>(compositeFlux.new_surfaceFlux());
   
   }
  else if (type == "normalizedFlux")
   {
     
     Real amplitude;
     pp.get("amplitude",amplitude);
     std::string prefix(a_prefix);
     prefix += ".direction";
     SurfaceFlux* direction = parse(prefix.c_str());
     if (direction == NULL)
       {
	 pout() << " no flux defined in " << prefix  << std::endl;
	 MayDay::Error("no flux defined in <normalized flux>.direction");
       }
       
     NormalizedFlux flux(direction, amplitude);
     ptr = static_cast<SurfaceFlux*>(flux.new_surfaceFlux());
   
   }
  else if (type == "targetThicknessFlux")
   {
     
     Real timescale = 1.0e+4;
     pp.get("timescale",timescale);
     std::string prefix(a_prefix);
     prefix += ".target";
     SurfaceFlux* target = parse(prefix.c_str());
     if (target == NULL)
       {
	 pout() << " no flux defined in " << prefix  << std::endl;
	 MayDay::Error("no flux defined in <target flux>.target");
	 
       }
     
     TargetThicknessFlux flux(target, timescale);
     ptr = static_cast<SurfaceFlux*>(flux.new_surfaceFlux());
   
   }
  else if (type == "groundingLineLocalizedFlux")
    {
      Real powerOfThickness = 0.0;
      pp.query("powerOfThickness",powerOfThickness);

      std::string glPrefix(a_prefix);
      glPrefix += ".groundingLine";
      SurfaceFlux* glPtr = parse(glPrefix.c_str());
      if (glPtr == NULL)
	{
	  glPtr = new zeroFlux;
	}
       
      std::string ambientPrefix(a_prefix);
      ambientPrefix += ".ambient";
      SurfaceFlux* ambientPtr = parse(ambientPrefix.c_str());
      if (ambientPtr == NULL)
	{
	  ambientPtr = new zeroFlux;
	}
      
      ptr = static_cast<SurfaceFlux*>
	(new GroundingLineLocalizedFlux(glPtr->new_surfaceFlux(),
					ambientPtr->new_surfaceFlux(),
					powerOfThickness ));
	 
      delete glPtr;
      delete ambientPtr;

    }
  else if (type == "FloatingDivUHLocalizedFlux")
    {
      std::string prefix(a_prefix);
      prefix +=  ".flux";
      SurfaceFlux* fluxptr = parse(prefix.c_str());

      //ParmParse pp2(prefix.c_str());
      Real dx;
      pp.get("mesh_spacing",dx);
      
      ptr = static_cast<SurfaceFlux*>
	(new  FloatingDivUHLocalizedFlux(fluxptr->new_surfaceFlux(), dx));
      delete fluxptr;
      
    }
  else if ((type == "IMSIP6OceanForcing") || (type == "ISMIP6OceanForcing"))
    {
      ptr = new ISMIP6OceanForcing(pp);
    }
#ifdef HAVE_PYTHON
  else if (type == "pythonFlux") {
    
    std::string module;
    pp.get("module",module);
    std::string function;
    pp.get("function",function);

    int nkwargs = 0;
    pp.query("n_kwargs",nkwargs);
    Vector<std::string> kwargName(nkwargs);
    if (nkwargs > 0)
      {
	pp.queryarr("kwargs",kwargName,0,nkwargs);
      }
    std::map<std::string, Real> kwarg;
    for (int i = 0; i < nkwargs; i++)
      {
	kwarg[kwargName[i]] = 0.0;
      }
    PythonInterface::PythonSurfaceFlux pythonFlux(module, function, kwarg);
    ptr = static_cast<SurfaceFlux*>(pythonFlux.new_surfaceFlux());

  }
#endif
  else if (type == "")
    {
      ptr = NULL; // return a NULL and leave it up to the caller to care
    }
  else
    {
      // a type was specified but it made no sense...
      pout() << "unknown flux type " << type << std::endl;
      MayDay::Error("unknown flux type");
    }
  return ptr;
  
}




#ifdef HAVE_PYTHON
#include "signal.h"


#endif
#include "NamespaceFooter.H"
