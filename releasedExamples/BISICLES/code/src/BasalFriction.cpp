#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BasalFriction.H"
#include "IceConstants.H"
#include "ParmParse.H"
#include "twistyStreamFriction.H"
#include "singularStreamFriction.H"
#include "GaussianBumpFriction.H"
#include "LevelDataBasalFriction.H"
#include "ReadLevelData.H"
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "CONSTANTS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// -----------------------------------------------------------------
//   zeroFriction
// -----------------------------------------------------------------
BasalFriction* 
zeroFriction::new_basalFriction() const
{
  zeroFriction* newPtr = new zeroFriction;
  return static_cast<BasalFriction*>(newPtr);
}


void zeroFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
                               LevelSigmaCS& a_coordSys,
                               Real a_time,
                               Real a_dt)
{
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_betaSqr[dit].setVal(0.0);
    }
}

// -----------------------------------------------------------------
//   constantFriction
// -----------------------------------------------------------------
constantFriction::constantFriction() : m_isValSet(false)
{
}

constantFriction::constantFriction(Real a_value)
  : m_frictionVal(a_value), m_isValSet(true)
{

}

BasalFriction* 
constantFriction::new_basalFriction() const
{
  constantFriction* newPtr = new constantFriction;
  newPtr->m_frictionVal = m_frictionVal;
  newPtr->m_isValSet = m_isValSet;
  return static_cast<BasalFriction*>(newPtr);
}

/// define source term for thickness evolution and place it in flux
/** dt is included in case one needs integrals or averages over a
    timestep
*/
void
constantFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
				   LevelSigmaCS& a_coordSys,
				   Real a_time,
				   Real a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_betaSqr[dit].setVal(m_frictionVal);
    }
}


///
void
constantFriction::setFrictionVal(const Real& a_betaSqr) 
{
  m_frictionVal = a_betaSqr; 
  m_isValSet = true;
}


// -----------------------------------------------------------------
//   sinusoidalFriction
// -----------------------------------------------------------------



/// constructor
sinusoidalFriction::sinusoidalFriction()
{
}

sinusoidalFriction::sinusoidalFriction(const Real& a_betaVal,
                                       const RealVect& a_omega,
                                       const Real& a_eps,
                                       const RealVect& a_domainSize) 
  : m_betaVal(a_betaVal), m_omega(a_omega), m_eps(a_eps), m_domainSize(a_domainSize)
{

}
                                                 

/// destructor
sinusoidalFriction::~sinusoidalFriction()
{
}

/// factory method
/** return a pointer to a new BasalFriction object
 */
BasalFriction* 
sinusoidalFriction::new_basalFriction() const
{
  sinusoidalFriction* newPtr = new sinusoidalFriction;
  newPtr->m_betaVal = m_betaVal;
  newPtr->m_omega = m_omega;
  newPtr->m_eps = m_eps;
  newPtr->m_domainSize = m_domainSize;
  return static_cast<BasalFriction*>(newPtr);
}

/// define source term for thickness evolution and place it in friction
/** dt is included in case one needs integrals or averages over a
    timestep. friction should be defined in meters/second in the current 
    implementation. 
*/
void
sinusoidalFriction::setBasalFriction(LevelData<FArrayBox>& a_betaSqr,
				     LevelSigmaCS& a_coordSys,
				     Real a_time,
				     Real a_dt)
{
  RealVect twoPiOmega(m_omega);
  twoPiOmega /= m_domainSize;
  twoPiOmega *= 2*Pi;

  RealVect dx = a_coordSys.dx();

  DataIterator dit = a_betaSqr.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisBeta = a_betaSqr[dit];

      BoxIterator bit(thisBeta.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;
          
          Real betaVal = m_betaVal*(1.0+m_eps+D_TERM(cos(twoPiOmega[0]*x[0]),
                                                     *cos(twoPiOmega[1]*x[1]),
                                                     *cos(twoPiOmega[2]*x[2])));
          thisBeta(iv,0) = betaVal;
        } // end loop over cells
    } // end loop over grids

}

/// set friction value in Pa*a/m)
void
sinusoidalFriction::setSinParameters(const Real& a_betaVal, 
                                     const RealVect& a_omega,
                                     const Real& a_eps,
                                     const RealVect& a_domainSize)
{
  m_betaVal = a_betaVal;
  m_omega = a_omega;
  m_eps = a_eps;
  m_domainSize = a_domainSize;
}



BasalFriction* BasalFriction::parse(const char* a_prefix, const RealVect& a_domainSize)
{
  BasalFriction* basalFrictionPtr = NULL;
  std::string type = "";

  ParmParse pp(a_prefix);

  pp.query("beta_type",type);
    

  if (type == "constantBeta")
    {
      Real betaVal;
      pp.get("betaValue", betaVal);
      basalFrictionPtr = static_cast<BasalFriction*>(new constantFriction(betaVal));
    }
  else if (type == "sinusoidalBeta")
    {
      Real betaVal, eps;
      RealVect omega(RealVect::Unit);
      Vector<Real> omegaVect(SpaceDim);
      pp.get("betaValue", betaVal);
      if (pp.contains("omega"))
	{
	  pp.getarr("omega", omegaVect, 0, SpaceDim);
	  omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
	}
      pp.get("betaEps", eps);
      basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
									    omega, 
									    eps,
									    a_domainSize));
    }
  // keep this one around for backward compatibility, even if it
  // is a special case of sinusoidalBeta
  else if (type == "sinusoidalBetay")
    {
      Real betaVal, eps, omegaVal;
      RealVect omega(RealVect::Zero);
      omega[1] = 1;
        
      pp.get("betaValue", betaVal);
      if (pp.contains("omega"))
	{
	  pp.get("omega", omegaVal);
	  omega[1] = omegaVal;
	}
      pp.get("betaEps", eps);
      basalFrictionPtr = static_cast<BasalFriction*>(new sinusoidalFriction(betaVal, 
									    omega, 
									    eps,
									    a_domainSize));

    }
  else if (type == "twistyStreamx")
    {
      Real betaVal, eps, magOffset;
      magOffset = 0.25;
      RealVect omega(RealVect::Unit);
      Vector<Real> omegaVect(SpaceDim);
      pp.get("betaValue", betaVal);
      if (pp.contains("omega"))
	{
	  pp.getarr("omega", omegaVect, 0, SpaceDim);
	  omega = RealVect(D_DECL(omegaVect[0], omegaVect[1], omegaVect[2]));
	}
      pp.query("magOffset", magOffset);
      pp.get("betaEps", eps);
      basalFrictionPtr = static_cast<BasalFriction*>(new twistyStreamFriction(betaVal, 
									      omega, 
									      magOffset, 
									      eps,
									      a_domainSize));          
    }
  else if (type == "singularStream")
    {
      Real slippyC,stickyC,width,twistNumber,twistAmplitude;
      pp.get("slippyC",slippyC);
      pp.get("stickyC",stickyC);
      pp.get("width",width);
      pp.get("twistNumber",twistNumber);
      pp.get("twistAmplitude",twistAmplitude);
      basalFrictionPtr = static_cast<BasalFriction*>
	(new singularStreamFriction
	 (slippyC,stickyC,width,twistNumber,twistAmplitude,a_domainSize));

    }
  else if (type == "gaussianBump")
    {
      int nt;
      pp.get("gaussianBump_nt", nt);
      Vector<Real> t(nt-1);
      Vector<Real> C0(nt),a(nt);
      Vector<RealVect> b(nt),c(nt);

      pp.getarr("gaussianBump_t", t, 0, nt-1);
      pp.getarr("gaussianBump_C", C0, 0, nt);
      pp.getarr("gaussianBump_a", a, 0, nt);
     
#if CH_SPACEDIM == 1
      Vector<Real> xb(nt),xc(nt);
      pp.getarr("gaussianBump_xb", xb, 0, nt);
      pp.getarr("gaussianBump_xc", xc, 0, nt);
      for (int i = 0; i < nt; ++i)
	{
	  b[i][0] = xb[i];
	  c[i][0] = xc[i];
	}
#elif CH_SPACEDIM == 2
      Vector<Real> xb(nt),yb(nt),xc(nt),yc(nt);
      pp.getarr("gaussianBump_xb", xb, 0, nt);
      pp.getarr("gaussianBump_xc", xc, 0, nt);
      pp.getarr("gaussianBump_yb", yb, 0, nt);
      pp.getarr("gaussianBump_yc", yc, 0, nt);
      for (int i = 0; i < nt; ++i)
	{
	  b[i][0] = xb[i];
	  b[i][1] = yb[i];
	  c[i][0] = xc[i];
	  c[i][1] = yc[i];
	}
#else
      MayDay::Error("type = gaussianBump not implemented for CH_SPACEDIM > 2")
#endif
        basalFrictionPtr = static_cast<BasalFriction*>
	(new GaussianBumpFriction(t, C0, a, b, c));
    }
  else if (type == "LevelData")
    {
      //read a one level beta^2 from an AMR Hierarchy, and  store it in a LevelDataBasalFriction
      ParmParse ildPP("inputLevelData");
      std::string frictionName = "btrc";
      ildPP.query("frictionName",frictionName);
      
      std::string infileFormat = "";
      ildPP.query("frictionFileFormat",infileFormat);
     
      if (infileFormat == "")
	{
	  //basic case : one file is specified. this might be unnecessary since the multiple file case
	  //includes the one file case, but for now we are holding onto the older code
	  std::string infile = "";
	  ildPP.get("frictionFile",infile);
	  
	  RefCountedPtr<LevelData<FArrayBox> > levelC (new LevelData<FArrayBox>());
	  Real dx;
	  
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectC;
	  vectC.push_back(levelC);
	  
	  Vector<std::string> names(1);
	  names[0] = frictionName;
	  
	  readLevelData(vectC,dx,infile,names,1);
	  
	  RealVect levelDx = RealVect::Unit * dx;
	  basalFrictionPtr = static_cast<BasalFriction*>
	    (new LevelDataBasalFriction(levelC,levelDx));
	}
      else
	{
	  int n = 1;
	  ildPP.query("frictionFileSteps",n);
	  
	  int startTime = 0, timeStep = 1;
	  ildPP.query("frictionFileStartTime", startTime);
	  ildPP.query("frictionFileTimeStep", timeStep);
	  RefCountedPtr<std::map<Real,std::string> > tf
	    (new std::map<Real,std::string>);
	  for (int i =0; i < n; i++)
	    {
	      char* file = new char[infileFormat.length()+32];
	      int j = startTime + i * timeStep;
	      sprintf(file, infileFormat.c_str(),j);
	      tf->insert(make_pair(Real(j), file));
	      delete file;
	    }
	  
	  basalFrictionPtr = static_cast<BasalFriction*>
	    (new LevelDataBasalFriction(tf,frictionName));
	}
    }
  else if (type == "MultiLevelData")
    {
      //read a multi level beta^2 from an AMR Hierarchy, and  store it in a MultiLevelDataBasalFriction
      ParmParse ildPP("inputLevelData");
      std::string infile;
      ildPP.get("frictionFile",infile);
      std::string frictionName = "btrc";
      ildPP.query("frictionName",frictionName);
      Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > vectC;
      Real dx;
      Vector<int> ratio;
      Vector<std::string> names(1);
      names[0] = frictionName;
      readMultiLevelData(vectC,dx,ratio,infile,names,1);
      RealVect dxCrse = RealVect::Unit * dx;
      basalFrictionPtr = static_cast<BasalFriction*>
	(new MultiLevelDataBasalFriction(vectC[0],dxCrse,ratio));
    }
#ifdef HAVE_PYTHON
  else if (type == "Python")
    {
      ParmParse pyPP("PythonBasalFriction");
      std::string module;
      pyPP.get("module",module);
      std::string funcName = "friction";
      pyPP.query("function",funcName);
      basalFrictionPtr = static_cast<BasalFriction*>
	(new PythonInterface::PythonBasalFriction(module, funcName));

    }
#endif
  return basalFrictionPtr;
}
#include "NamespaceFooter.H"
