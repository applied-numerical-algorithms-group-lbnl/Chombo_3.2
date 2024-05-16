#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifdef HAVE_PYTHON

// PythonInterface.cpp
// ============
// BISICLES embedded Python interface. 

#include "PythonInterface.H"
#include "ReflectGhostCells.H"
#include "FillFromReference.H"
#include "CoarseAverage.H"
#include "IceConstants.H"
#include "AmrIceBase.H"
#include "ParmParse.H"
#include "signal.h"
#include "VelocityBC.H"
#include "IceThermodynamics.H"
#include "NamespaceHeader.H"

void PythonInterface::PythonEval(PyObject* a_pFunc, 
				 Vector<Real>& a_value, 
				 Vector<Real>& a_arg, 
				 std::map<std::string, Real>* a_kwarg)
{
  PyObject *args, *kwargs, *value;
  args = PyTuple_New(a_arg.size());

  //contruct positional args
  for (int j = 0; j < a_arg.size(); ++j)
    {
      value = PyFloat_FromDouble(a_arg[j]);
      if (!value)
	{
	  pout() << "PythonInterface::PythonEval failed to construct Python args" << endl;
	  CH_assert(value != NULL);
	  MayDay::Error("PythonInterface::PythonEval failed to construct Python args");
	}
      PyTuple_SetItem(args, j, value);
    }
  
  //contruct keyword args
  kwargs = NULL;
  if (a_kwarg)
    {
      std::map<std::string,Real>& kwarg = *a_kwarg;
      if (kwarg.size() > 0)
	{
	  kwargs = PyDict_New();
	  for (std::map<std::string, Real>::const_iterator j = kwarg.begin(); j!= kwarg.end(); ++j)
	    {
	      value = PyFloat_FromDouble(j->second);
	      if (!value)
		{
		  pout() << "PythonInterface::PythonEval failed to construct Python kwargs" << endl;
		  CH_assert(value != NULL);
		  MayDay::Error("PythonInterface::PythonEval failed to construct Python kwargs");
		}
	      PyDict_SetItemString(kwargs,j->first.c_str(),value);
	    }
	}
    }
  
  value = PyObject_Call(a_pFunc, args, kwargs);
  Py_DECREF(args);
  if (kwargs)
    {Py_DECREF(kwargs);}
  if (!value)
    {
      pout() << "PythonInterface::PythonEval  python call failed" << endl;
      PyErr_Print();
      CH_assert(value != NULL);
      MayDay::Error("PythonInterface::PythonEval python call failed");
    }

  bool value_ok = false;
  if (PyFloat_CheckExact(value))
    {
      a_value.resize(1,0.0);
      a_value[0] =  PyFloat_AS_DOUBLE(value);
      value_ok = true;
    }
  else if (PyTuple_CheckExact(value))
    {
      Py_ssize_t n = PyTuple_Size(value);
      a_value.resize(n);
      value_ok = true;
      for (int i = 0; i < n; i++)
	{
	  PyObject* t = PyTuple_GetItem(value,i);
	  if (PyFloat_CheckExact(t))
	    {
	      a_value[i] =  PyFloat_AS_DOUBLE(t);
	    }
	  else
	    {
	      value_ok = false;
	    }
	}
    }

  Py_DECREF(value);
  
  if (!value_ok)
    {
      MayDay::Error("PythonInterface::PythonEval python return value is not a float or tuple of floats");
    }

}

void PythonInterface::InitializePythonModule(PyObject **a_pModule,
					     const std::string& a_pyModuleName)
{
  if ( !Py_IsInitialized() )
    {
      //initialize Python without overriding SIGINT  
      PyOS_sighandler_t sigint = PyOS_getsig(SIGINT);
      Py_Initialize();
      sigint =  PyOS_setsig(SIGINT, sigint);
    }
  
  PyObject* pName = PyUnicode_FromString(a_pyModuleName.c_str());
  *a_pModule =  PyImport_Import(pName);

  if (*a_pModule == NULL)
    {
      PyErr_Print();
      pout() << "failed to import Python module " << a_pyModuleName << endl;
      CH_assert(a_pModule != NULL);
      MayDay::Error("failed to import Python module ");
    }

  Py_DECREF(pName);

}

void PythonInterface::InitializePythonFunction(  PyObject **a_pFunc,
						 PyObject*  a_pModule,
						 const std::string& a_pyFuncName)
{

 
  pout() << "PythonInterface::InitializePythonFunction initializing " 
	 << a_pyFuncName << endl;
  
  if ( a_pModule == NULL)
    {
      pout() << "null python module " << a_pyFuncName << endl;
      CH_assert(a_pModule != NULL);
      MayDay::Error("null python module");
    }

  char* t = const_cast<char*>(a_pyFuncName.c_str());
  *a_pFunc = PyObject_GetAttrString(a_pModule,t);
  
  if (a_pFunc == NULL)
    {
      pout() << "failed to initialize Python function " << a_pyFuncName << endl;
      CH_assert(*a_pFunc != NULL);
      MayDay::Error("failed to initialize Python function");
    }

  if (!PyCallable_Check(*a_pFunc))
    {
      pout() << "Python function " << a_pyFuncName << " not callable " <<  endl;
      CH_assert(PyCallable_Check(*a_pFunc));
      MayDay::Error("Python function not callable");
    }


}



PythonInterface::PythonIBC::PythonIBC(const std::string& a_pyModuleName,
				      const std::string& a_pyFuncThckName,
				      const std::string& a_pyFuncTopgName, 
				      const std::string& a_pyFuncRhsName,
				      const std::string& a_pyFuncFaceVelName)
  :m_isBCsetUp(false),m_verbose(true)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFuncThck, m_pModule,  a_pyFuncThckName);
  InitializePythonFunction(&m_pFuncTopg, m_pModule,  a_pyFuncTopgName);
  if (a_pyFuncRhsName != "")
    {
      InitializePythonFunction(&m_pFuncRhs, m_pModule,  a_pyFuncRhsName);
    }
  else
    {
      m_pFuncRhs = NULL;
    }
  if (a_pyFuncFaceVelName != "")
    {
      InitializePythonFunction(&m_pFuncFaceVel, m_pModule,  a_pyFuncFaceVelName);
    }
  else
    {
      m_pFuncFaceVel = NULL;
    }


}


PythonInterface::PythonIBC::PythonIBC()
  :m_isBCsetUp(false),m_verbose(true)
{

}

PythonInterface::PythonIBC::~PythonIBC()
{
  
  
  Py_DECREF(m_pFuncThck);
  Py_DECREF(m_pFuncTopg);
  if (m_pFuncRhs != NULL)
    {
      Py_DECREF(m_pFuncRhs);
    }
  if (m_pFuncRhs != NULL)
    {
      Py_DECREF(m_pFuncFaceVel);
    }
  Py_DECREF(m_pModule);
}

void PythonInterface::PythonIBC::define(const ProblemDomain& a_domain,
		       const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
  //m_dx = a_dx * RealVect::Unit;
}

void PythonInterface::PythonIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("PythonIBC::initialize not implemented");
}

// set boundary fluxes
void PythonInterface::PythonIBC::primBC(FArrayBox&            a_WGdnv,
			   const FArrayBox&      a_Wextrap,
			   const FArrayBox&      a_W,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const Real&           a_time)
{
 if (!m_domain.isPeriodic(a_dir))
   {
     int lohisign;
     Box tmp = a_WGdnv.box();
     // Determine which side and thus shifting directions
     lohisign = sign(a_side);
     tmp.shiftHalf(a_dir,lohisign);
     // (DFM - 5/28/10) this little dance with the ghostBox is a bit 
     // of a kluge to handle the case where a_WGdnv has more than one layer   
     // of ghosting, in which case just testing agains tmp isn't 
     // sufficient to determine whether you're up against the domain edge
     Box ghostBox = (a_side == Side::Lo)?
       adjCellLo(m_domain.domainBox(),a_dir, 1):
       adjCellHi(m_domain.domainBox(),a_dir, 1);
     ghostBox &= tmp;
     // Is there a domain boundary next to this grid
     if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
       {
	 tmp &= m_domain;
	 Box boundaryBox = (a_side == Side::Lo)?
	   bdryLo(tmp,a_dir):bdryHi(tmp,a_dir);
	 BoxIterator bit(boundaryBox);
	 for (bit.begin(); bit.ok(); ++bit){
	   const IntVect& i = bit();
	   a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	 }
       }   

   }
}

void  PythonInterface::PythonIBC::setBdrySlopes(FArrayBox&       a_dW,
				  const FArrayBox& a_W,
				  const int&       a_dir,
				  const Real&      a_time)
{
  // one-sided differences sounds fine with me, so do nothing...
}


void PythonInterface::PythonIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here
  MayDay::Error("PythonIBC::artViscBC not implemented");
}

/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC> PythonInterface::PythonIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  //return m_velBCs; 

  return m_velBC;
}

void  iceDivideBC_PyBC(FArrayBox& a_vel,
					const Box& a_valid,
					const ProblemDomain& a_domain,
					Real a_dx,
					bool a_homogeneous)
{

  const IntVect ghostVect = IntVect::Unit;
   for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Lo);
	  ReflectGhostCells(a_vel, a_domain,ghostVect, dir, Side::Hi);

	  Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      ghostBoxLo &= a_vel.box();
	      //a_vel.mult(-1.0,ghostBoxLo,dir,1);
		  
	    }
	  Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      ghostBoxHi &= a_vel.box();
	      //a_vel.mult(-1.0,ghostBoxHi,dir,1);
	    }

	}
    }

}



void  PythonInterface::PythonIBC::setupBCs()
{
  Vector<RefCountedPtr<CompGridVTOBC> > loBC;
  Vector<RefCountedPtr<CompGridVTOBC> > hiBC;

  loBC.resize(SpaceDim);
  hiBC.resize(SpaceDim);

  ParmParse pp("PythonIBC");
  Vector<int> loBCType(SpaceDim),hiBCType(SpaceDim) ;

  //default to refelection BC's on all boundaries
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      loBCType[dir] = IceDivide;
      hiBCType[dir] = IceDivide;
    }
  pp.queryarr("bc_lo",loBCType,0,SpaceDim);
  pp.queryarr("bc_hi",hiBCType,0,SpaceDim);

  for (int dir = 0; dir < SpaceDim; dir++)
    {
      loBC[dir] = BCFactory(BCType(loBCType[dir]),dir,Side::Lo);
      hiBC[dir] = BCFactory(BCType(hiBCType[dir]),dir,Side::Hi);
    }
  m_velBC = RefCountedPtr<CompGridVTOBC>(new PerEdgeCompGridVTOBC(loBC, hiBC));
  m_isBCsetUp = true;

}

/// set non-periodic ghost cells for surface height z_s. 
void  PythonInterface::PythonIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
				       LevelSigmaCS& a_coords,
				       const ProblemDomain& a_domain,
				       const RealVect& a_dx,
				       Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
	}
    }
  a_zSurface.exchange();
}


/// set non-periodic ghost cells for ice fraction
void  PythonInterface::PythonIBC::setIceFractionBCs(LevelData<FArrayBox>& a_fraction,
						    const ProblemDomain& a_domain,
						    const RealVect& a_dx,
						    Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir))){
	ReflectGhostCells(a_fraction, a_domain, dir, Side::Lo);
	ReflectGhostCells(a_fraction, a_domain, dir, Side::Hi);	
      }
    }
  a_fraction.exchange();
}


/// set non-periodic ghost cells for thickness & topography
void  PythonInterface::PythonIBC::setGeometryBCs(LevelSigmaCS& a_coords,
						 const ProblemDomain& a_domain,
						 const RealVect& a_dx,
						 Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(a_domain.isPeriodic(dir))){
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
	ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
      }
    }
  a_coords.getTopography().exchange();
  a_coords.getH().exchange();
}

void  PythonInterface::PythonIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
					 const RealVect& a_dx,
					 const RealVect& a_domainSize,
					 const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{
  CH_TIME("PythonInterface::PythonIBC::initializeIceGeometry");
  if (m_verbose)
    {
      pout() << " PythonIBC::initializeIceGeometry" << endl;
    }

  
  Vector<Real> args(SpaceDim);
  Vector<Real> rval;
  const DisjointBoxLayout& grids = a_coords.grids();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      FArrayBox& topg = a_coords.getTopography()[dit];
      FArrayBox& thck = a_coords.getH()[dit];

      Box b = topg.box();
      b &= thck.box();
      b &= a_coords.grids().physDomain();
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  IntVect iv = bit();
	  int i = 0;
	   for (int dir=0; dir < SpaceDim; dir++)
	     {
	       args[i] = (Real(iv[dir]) + 0.5)*a_dx[dir];
	       i++;
	     }
	   PythonEval(m_pFuncThck, rval,  args);
	   thck(iv) = rval[0];
	   PythonEval(m_pFuncTopg, rval,  args);
	   topg(iv) = rval[0]; 	   
	}
    }
  a_coords.getTopography().exchange();
  a_coords.getH().exchange();
  Real dt = 0.0;
  setGeometryBCs(a_coords, grids.physDomain(), a_dx, a_time, dt);  
}

bool  PythonInterface::PythonIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				     const RealVect& a_dx,
				     const RealVect& a_domainSize,
				     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio)
{
  CH_TIME("PythonInterface::PythonIBC::regridIceGeometry");
  Vector<Real> args(SpaceDim);
  Vector<Real> rval;
  const DisjointBoxLayout& grids = a_coords.grids();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
     
      FArrayBox& topg = a_coords.getTopography()[dit];
      Box b = grids[dit];
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	   IntVect iv = bit();
	   int i = 0;
	   for (int dir=0; dir < SpaceDim; dir++)
	     {
	       args[i] = (Real(iv[dir]) + 0.5)*a_dx[dir];
	       i++;
	     }
	   
	   PythonEval(m_pFuncTopg, rval,  args);
	   topg(iv) = rval[0]; 
	}
    }

  a_coords.getTopography().exchange();
  a_coords.getH().exchange();

  // pythonIBC always regrids topography
  return true;
}

void PythonInterface::PythonIBC::modifyFaceVelocity
(LevelData<FluxBox>& a_faceVel,
 const LevelSigmaCS& a_coords,
 const ProblemDomain& a_domain) const
{

   CH_TIME("PythonInterface::PythonIBC::modifyFaceVelocity");
   if (m_pFuncFaceVel != NULL)
     {
       //args are : x,[y,],dx[dir],dir,facevel[dir]
       Vector<Real> args(SpaceDim + 3);
       Vector<Real> rval;
       const DisjointBoxLayout& grids = a_coords.grids();
       for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      FArrayBox& fvel = a_faceVel[dit][dir];
	      Box b = fvel.box();
	      for (BoxIterator bit(b); bit.ok(); ++bit)
		{
		  IntVect iv = bit();
		  int i = 0;
		 
		  for (int idir=0; idir < SpaceDim; idir++)
		    {
		      args[i] = (Real(iv[idir]) + ( (idir==dir)?0.0:0.5 ) )*a_coords.dx()[idir];
		      i++;
		    }
		  
		  args[i] = a_coords.dx()[dir] ; i++;
		  args[i] = dir; i++;
		  args[i] = fvel(iv); i++;

		  PythonEval(m_pFuncFaceVel, rval,  args);
		  fvel(iv) = rval[0];

		  if ( (dir == 1) && (a_coords.dx()[0] < 1500.0))
		    {
		      CH_assert(0 == 0);
		      int dbg = 0; 
		      dbg++;
		    }

		}
	      
	    }
	}

     }
}


/// if appropriate, modify velocity solve RHS in a problem-dependent way. 
  /** default is to do nothing
   */
void  PythonInterface::PythonIBC::modifyVelocityRHS
(LevelData<FArrayBox>& a_rhs,
 LevelSigmaCS& a_coords,
 const ProblemDomain& a_domain,
 Real a_time, Real a_dt)
{
  CH_TIME("PythonInterface::PythonIBC::modifyVelocityRHS");
  if (m_pFuncRhs != NULL)
    {
      //args are : x,[y,],dx[dir],dir,rhs,rho*g*H,s,sleft,sright
      Vector<Real> args(SpaceDim + 7);  
      Vector<Real> rval;
      const DisjointBoxLayout& grids = a_coords.grids();
      Real rhog = a_coords.iceDensity() * a_coords.gravity();

      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
	  const FArrayBox& usrf = a_coords.getSurfaceHeight()[dit];
	  const FArrayBox& thck = a_coords.getH()[dit];
	  FArrayBox& rhs = a_rhs[dit];
	  Box b = grids[dit];
	  
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      for (BoxIterator bit(b);bit.ok();++bit)
		{
		  IntVect iv = bit();
		  int i = 0;
		  // x , y
		  for (int idir=0; idir < SpaceDim; idir++)
		    {
		      args[i] = (Real(iv[idir]) + 0.5)*a_coords.dx()[idir];
		      i++;
		    }
		  // dx,
		  args[i] = a_coords.dx()[dir] ; i++;
		  //dir
		  args[i] = dir; i++;
		  //rhs
		  args[i] = rhs(iv,dir); i++;
		  // rho * g * H
		  args[i] = thck(iv) * rhog; i++;
		  // s, sleft, sright
		  args[i] = usrf(iv); i++;
		  args[i] = usrf(iv - BASISV(dir)); i++;
		  args[i] = usrf(iv + BASISV(dir) ); i++;

		  PythonEval(m_pFuncRhs, rval,  args);
		  rhs(iv,dir) = rval[0];
		  // if (dir == 1)
		  //   {
		  //     CH_assert(std::abs( rhs(iv,1) < 1.0e+6));
		  //   }
	   
		}
	    }
	}
    }
}

IceThicknessIBC*  PythonInterface::PythonIBC::new_thicknessIBC()
{
  PythonIBC* retval = new PythonIBC( m_pModule, m_pFuncThck, m_pFuncTopg, m_pFuncRhs, m_pFuncFaceVel);
  return static_cast<IceThicknessIBC*>(retval);
}


PythonInterface::PythonSurfaceFlux::PythonSurfaceFlux(const std::string& a_pyModuleName,
						      const std::string& a_pyFuncName, 
						      std::map<std::string,Real>& a_kwarg)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
  m_kwarg = a_kwarg;
}

PythonInterface::PythonSurfaceFlux::~PythonSurfaceFlux()
{
  Py_DECREF(m_pFunc);
  Py_DECREF(m_pModule);
  
}

SurfaceFlux* PythonInterface::PythonSurfaceFlux::new_surfaceFlux()
{
  PythonSurfaceFlux* ptr = new PythonSurfaceFlux(m_pModule, m_pFunc, m_kwarg);
  return static_cast<SurfaceFlux*>( ptr);
}

/**
   Map names to cell centered data given (a_amrIce,a_level,a_dit,a_iv).
   Silently sets a_value[i] to 0.0 if a_name[i] is not known
*/
void  PythonInterface::FillKwargs(std::map<std::string,Real>& a_kwarg,
				       const AmrIceBase& a_amrIce, int a_level, 
				       const DataIterator& a_dit, const IntVect& a_iv)
{
  //\todo : clumsy, since it does a series of string matches on a cell-by-cell basis, 
  // Might want to carry this on a FAB-by-FAB basis,
  // but on a grander scale a standard 'look-up-by-name' addition to AmrIceBase itself
  // might be useful. For now, this function is not widely used
  
  for (std::map<std::string,Real>::iterator i = a_kwarg.begin(); i!= a_kwarg.end(); ++i)
    {
      if (i->first == "gl_proximity")
	{
	  i->second = (*a_amrIce.groundingLineProximity(a_level))[a_dit](a_iv);
	} 
      else if (i->first == "gl_proximity_scale")
	{
	  i->second = a_amrIce.groundingLineProximityScale();
	}
      else
	{
	  i->second = 0.0;
	}
    }
}


void PythonInterface::PythonSurfaceFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
							      const AmrIceBase& a_amrIce, 
							      int a_level, Real a_dt)
{
   CH_TIME("PythonInterface::PythonSurfaceFlux::surfaceThicknessFlux");
 
  Vector<Real> args(SpaceDim + 3);
  Vector<Real> rval;
  const RefCountedPtr<LevelSigmaCS> levelCS = a_amrIce.geometry(a_level);

  const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      a_flux[dit].setVal(0.0);
      const Box& b = a_flux[dit].box();//grids[dit];
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  int i = 0;
	  for (int dir=0; dir < SpaceDim; dir++)
	    {
	      args[i] = (Real(iv[dir]) + 0.5)*levelCS->dx()[dir];
	      i++;
	    }
	  args[i] = a_amrIce.time();i++;
	  args[i] = levelCS->getH()[dit](iv);i++;
	  args[i] = levelCS->getTopography()[dit](iv);i++;
	  FillKwargs(m_kwarg, a_amrIce, a_level,  dit, iv);
	  PythonEval(m_pFunc, rval,  args, &m_kwarg);
	  a_flux[dit](iv) =  rval[0];

	}
    }
}


PythonInterface::PythonIceTemperatureIBC::PythonIceTemperatureIBC
(const std::string& a_pyModuleName,const std::string& a_pyFuncName)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
}

PythonInterface::PythonIceTemperatureIBC::~PythonIceTemperatureIBC()
{
  Py_DECREF(m_pFunc);
  Py_DECREF(m_pModule);
  
}

PythonInterface::PythonIceTemperatureIBC* 
PythonInterface::PythonIceTemperatureIBC::new_internalEnergyIBC()
{
  PythonInterface::PythonIceTemperatureIBC* ptr = 
    new PythonInterface::PythonIceTemperatureIBC(m_pModule, m_pFunc);
  return ptr;
}

void
PythonInterface::PythonIceTemperatureIBC::basalHeatFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  MayDay::Error("PythonInterface::PythonIceTemperatureIBC::basalHeatFlux not implemented");
}

#if BISICLES_Z == BISICLES_LAYERED
void PythonInterface::PythonIceTemperatureIBC::initializeIceInternalEnergy
(LevelData<FArrayBox>& a_E,  LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceE, 
 LevelData<FArrayBox>& a_basalE,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
 CH_TIME("PythonInterface::PythonIceTemperatureIBC::initializeIceTemperature");
 //function args are (x,[y,],thickness,topography,sigma)
 Vector<Real> args(SpaceDim + 3);
 Vector<Real> rval;
  
 const DisjointBoxLayout& grids = a_E.disjointBoxLayout();
 const LevelSigmaCS& coordSys = *a_amrIce.geometry(a_level); 
 for (DataIterator dit(grids);dit.ok();++dit)
   {
     
     FArrayBox T( a_E[dit].box(), a_E[dit].nComp()); //temperature
     FArrayBox w( a_E[dit].box(), a_E[dit].nComp()); //water fraction, set to zero for now
     w.setVal(0.0);
     FArrayBox sT( a_surfaceE[dit].box(), 1); //surface temperature
     FArrayBox sw( a_surfaceE[dit].box(), 1); //surface water fraction, set to zero for now
     sw.setVal(0.0);
     FArrayBox bT( a_basalE[dit].box(), 1); //basal temperature
     FArrayBox bw( a_basalE[dit].box(), 1); //basal water fraction, set to zero for now
     bw.setVal(0.0);
     const Box& b = a_E[dit].box();//grids[dit];
     for (BoxIterator bit(b);bit.ok();++bit)
       {
	 const IntVect& iv = bit();
	 int i = 0;
	 for (int dir=0; dir < SpaceDim; dir++)
	   {
	     args[i] = (Real(iv[dir]) + 0.5)*coordSys.dx()[dir];
	     i++;
	   }
	 args[i] = coordSys.getH()[dit](iv);i++;
	 args[i] = coordSys.getTopography()[dit](iv);i++;
	 //layer midpoint temperatures
	 for (int layer = 0; layer < T.nComp(); layer++)
	   {
	     args[i] = coordSys.getSigma()[layer];
	     PythonEval(m_pFunc, rval,  args);
	     T(iv,layer) =  rval[0];
	   }
	 //surface temperature
	 args[i] = 0.0;
	 PythonEval(m_pFunc, rval,  args);
	 sT(iv) =  rval[0];
	 //basal temperature
	 args[i] = 1.0;
	 PythonEval(m_pFunc, rval,  args);
	 bT(iv) =  rval[0];
       }
     IceThermodynamics::composeInternalEnergy(a_E[dit],T,w,a_E[dit].box());
     IceThermodynamics::composeInternalEnergy(a_surfaceE[dit],sT,w,a_E[dit].box());
     IceThermodynamics::composeInternalEnergy(a_basalE[dit],bT,w,a_E[dit].box());
     a_tillWaterDepth[dit].setVal(0.0);
   }
 
}

void PythonInterface::PythonIceTemperatureIBC::setIceInternalEnergyBC
(LevelData<FArrayBox>& a_E, 
 LevelData<FArrayBox>& a_tillWaterDepth,
 LevelData<FArrayBox>& a_surfaceT, 
 LevelData<FArrayBox>& a_basalT,
 const LevelSigmaCS& a_coordSys)
{
  const ProblemDomain& domain = a_coordSys.grids().physDomain();

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir)))
	{
	  ReflectGhostCells(a_E, domain, dir, Side::Lo);
	  ReflectGhostCells(a_E, domain, dir, Side::Hi);
	  ReflectGhostCells(a_tillWaterDepth, domain, dir, Side::Lo);
	  ReflectGhostCells(a_tillWaterDepth, domain, dir, Side::Hi);
	  ReflectGhostCells(a_surfaceT, domain, dir, Side::Lo);
	  ReflectGhostCells(a_surfaceT, domain, dir, Side::Hi);
	  ReflectGhostCells(a_basalT, domain, dir, Side::Lo);
	  ReflectGhostCells(a_basalT, domain, dir, Side::Hi);
	}
    }
  //\todo make an option to call python function
}

#elif BISICLES_Z == BISICLES_FULLZ
void PythonInterface::PythonIceTemperatureIBC::initializeIceTemperature
(LevelData<FArrayBox>& a_T,
 const LevelSigmaCS& a_coordSys)
{

  MayDay::Error("BISICLES_Z == BISICLES_FULLZ PythonIceTemperatureIBC::initializeIceTemperature not yet implemented");
}


void PythonInterface::PythonIceTemperatureIBC::setIceTemperatureBC
(LevelData<FArrayBox>& a_T,
 const LevelSigmaCS& a_coordSys)
{

  MayDay::Error("BISICLES_Z == BISICLES_FULLZ PythonIceTemperatureIBC::setIceTemperatureBC not yet implemented");
}

#endif



PythonInterface::PythonBasalFriction::PythonBasalFriction
(const std::string& a_pyModuleName,
 const std::string& a_pyFuncName)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
}

PythonInterface::PythonBasalFriction::~PythonBasalFriction()
{
  Py_DECREF(m_pFunc);
  Py_DECREF(m_pModule);
  
}

BasalFriction* PythonInterface::PythonBasalFriction::new_basalFriction() const
{
  PythonBasalFriction* ptr = new PythonBasalFriction(m_pModule, m_pFunc);
  return static_cast<BasalFriction*>(ptr);
}


void PythonInterface::PythonBasalFriction::setBasalFriction
(LevelData<FArrayBox>& a_C,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{
  CH_TIME("PythonInterface::PythonBasalFriction::setBasalFriction");
  Vector<Real> args(SpaceDim + 3);
  Vector<Real> rval;
  
  const DisjointBoxLayout& grids = a_C.disjointBoxLayout();
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      a_C[dit].setVal(0.0);
      const Box& b = a_C[dit].box();//grids[dit];
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  int i = 0;
	  for (int dir=0; dir < SpaceDim; dir++)
	    {
	      args[i] = (Real(iv[dir]) + 0.5)*a_coordSys.dx()[dir];
	      i++;
	    }
	  args[i] = a_time;i++;
	  args[i] = a_coordSys.getH()[dit](iv);i++;
	  args[i] = a_coordSys.getTopography()[dit](iv);i++;
	  PythonEval(m_pFunc, rval,  args);
	  a_C[dit](iv) =  rval[0];
	  
	}
    }
}



PythonInterface::PythonMuCoefficient::PythonMuCoefficient
(const std::string& a_pyModuleName,
 const std::string& a_pyFuncName)
{
  InitializePythonModule(&m_pModule,  a_pyModuleName);
  InitializePythonFunction(&m_pFunc, m_pModule,  a_pyFuncName);
}




PythonInterface::PythonMuCoefficient::~PythonMuCoefficient()
{
  Py_DECREF(m_pFunc);
  Py_DECREF(m_pModule);
 
}

MuCoefficient* PythonInterface::PythonMuCoefficient::new_muCoefficient() const
{
  PythonMuCoefficient* ptr = new PythonMuCoefficient(m_pModule, m_pFunc);
  return static_cast<MuCoefficient*>(ptr);
}


void PythonInterface::PythonMuCoefficient::setMuCoefficient
(LevelData<FArrayBox>& a_muCoef,
 LevelSigmaCS& a_coordSys,
 Real a_time,
 Real a_dt)
{
  CH_TIME("PythonInterface::PythonMuCoefficient::setMuCoefficient");
  Vector<Real> args(SpaceDim + 3);
  Vector<Real> rval;
  
  const DisjointBoxLayout& grids = a_muCoef.disjointBoxLayout();
  for (DataIterator dit(grids);dit.ok();++dit)
    {
      a_muCoef[dit].setVal(0.0);
      const Box& b = a_muCoef[dit].box();
      for (BoxIterator bit(b);bit.ok();++bit)
	{
	  const IntVect& iv = bit();
	  int i = 0;
	  for (int dir=0; dir < SpaceDim; dir++)
	    {
	      args[i] = (Real(iv[dir]) + 0.5)*a_coordSys.dx()[dir];
	      i++;
	    }
	  args[i] = a_time;i++;
	  args[i] = a_coordSys.getH()[dit](iv);i++;
	  args[i] = a_coordSys.getTopography()[dit](iv);i++;
	  PythonEval(m_pFunc, rval,  args);
	  a_muCoef[dit](iv) =  rval[0];
	  
	}
    }
}

void PythonInterface::PythonVelocitySolver::define
(const ProblemDomain& a_coarseDomain,
 ConstitutiveRelation* a_constRel,
 BasalFrictionRelation* a_basalFrictionRel,
 const Vector<DisjointBoxLayout>& a_vectGrids,
 const Vector<int>& a_vectRefRatio,
 const RealVect& a_dxCrse,
 IceThicknessIBC* a_bc,
 int a_numLevels)
{
  ParmParse pp("PythonVelocitySolver");
  std::string module;
  pp.get("module",module);
  InitializePythonModule(&m_pModule,  module);
  std::string function;
  pp.get("function",function);
  InitializePythonFunction(&m_pFunc, m_pModule,  function);
}

int PythonInterface::PythonVelocitySolver::solve
(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
 Vector<LevelData<FArrayBox>* >& a_calvedIce,
 Vector<LevelData<FArrayBox>* >& a_addedIce,
 Vector<LevelData<FArrayBox>* >& a_removedIce,
 Real& a_initialResidualNorm, Real& a_finalResidualNorm,
 const Real a_convergenceMetric,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 const Vector<LevelData<FArrayBox>* >& a_beta,
 const Vector<LevelData<FArrayBox>* >& a_beta0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 const Vector<LevelData<FArrayBox>* >& a_muCoef,
 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 Real a_time,
 int a_lbase, int a_maxLevel)
{
  CH_TIME("PythonInterface::PythonVelocitySolver::solve");

  for (int lev =0; lev <= a_maxLevel; lev++)
    { 
      LevelData<FArrayBox>& vel = *a_horizontalVel[lev];
      const DisjointBoxLayout& grids = vel.disjointBoxLayout();
      const LevelSigmaCS& coordSys = *a_coordSys[lev];
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  Vector<Real> args(SpaceDim + 3);
	  Vector<Real> rval;
	  const Box& b = vel[dit].box();
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      int i = 0;
	      for (int dir=0; dir < SpaceDim; dir++)
		{
		  args[i] = (Real(iv[dir]) + 0.5)*coordSys.dx()[dir];
		  i++;
		}
	      args[i] = a_time;i++;
	      args[i] = coordSys.getH()[dit](iv);i++;
	      args[i] = coordSys.getTopography()[dit](iv);i++;
	      PythonEval(m_pFunc, rval,  args);
	      CH_assert(rval.size() >= SpaceDim);
	      vel[dit](iv,0) =  rval[0];
	      vel[dit](iv,1) =  rval[1];
	    }
	  
	}
      
    }

  return 0;
}

PythonInterface::PythonVelocitySolver::~PythonVelocitySolver()
{
  Py_DECREF(m_pFunc);
  Py_DECREF(m_pModule);
}

#include "NamespaceFooter.H"

#endif
