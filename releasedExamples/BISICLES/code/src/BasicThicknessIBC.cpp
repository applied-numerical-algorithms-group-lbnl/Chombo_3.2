#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "BasicThicknessIBC.H"
#include "defineLevelSigmaCS.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"


void zeroBCValue(Real* pos,
                 int* dir,
                 Side::LoHiSide* side,
                 Real* a_values)
{
  a_values[0]=0.0;
}


void iceNeumannBC(FArrayBox& a_state,
                  const Box& a_valid,
                  const ProblemDomain& a_domain,
                  Real a_dx,
                  bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}


void iceDirichletBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Lo);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  //Real bcVal = 0.0;
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}




// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
BasicThicknessIBC::BasicThicknessIBC() : m_boundaryThickness(0.0)
{
  m_isBCsetUp = false;
  m_isDefined = false;
}

BasicThicknessIBC::~BasicThicknessIBC()
{
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
BasicThicknessIBC::define(const ProblemDomain& a_domain,
                          const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
BasicThicknessIBC::new_thicknessIBC()
{
  BasicThicknessIBC* retval = new BasicThicknessIBC();

  retval->m_boundaryThickness = m_boundaryThickness;

  return static_cast<IceThicknessIBC*>(retval);
}

/// Set up initial conditions
/**
 */
void
BasicThicknessIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("BasicThicknessIBC::initialize not implemented");
}

/// Set boundary fluxes
/**
 */
void 
BasicThicknessIBC::primBC(FArrayBox&            a_WGdnv,
                          const FArrayBox&      a_Wextrap,
                          const FArrayBox&      a_W,
                          const int&            a_dir,
                          const Side::LoHiSide& a_side,
                          const Real&           a_time)
{
  // do nothing in periodic case
  if (!m_domain.isPeriodic(a_dir))
    {    
      int lohisign;
      Box tmp = a_WGdnv.box();
      
      // Determine which side and thus shifting directions
      lohisign = sign(a_side);
      tmp.shiftHalf(a_dir,lohisign);
      
      // Is there a domain boundary next to this grid
      if (!m_domain.contains(tmp))
        {
          tmp &= m_domain;
          
          Box boundaryBox;
          
          if (a_side == Side::Lo)
            {
              boundaryBox = bdryLo(tmp,a_dir);
            }
          else
            {
              boundaryBox = bdryHi(tmp,a_dir);
            }
          
          // Set the boundary values
          a_WGdnv.setVal(m_boundaryThickness, boundaryBox, 0, 1);
        }
    }
  
}

/// Set boundary slopes
/**
   The boundary slopes in a_dW are already set to one sided difference
   approximations.  If this function doesn't change them they will be
   used for the slopes at the boundaries.
*/

void 
BasicThicknessIBC::setBdrySlopes(FArrayBox&       a_dW,
                                 const FArrayBox& a_W,
                                 const int&       a_dir,
                                 const Real&      a_time)
{
  // one-sided differences sounds fine with me, so do nothing...
}

/// Adjust boundary fluxes to account for artificial viscosity
/**
 */
void 
BasicThicknessIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("BasicThicknessIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
BasicThicknessIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs;
}



  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
      by calling the defineLevelSigmaCS function in the util directory
   */
void
BasicThicknessIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                         const RealVect& a_dx,
                                         const RealVect& a_domainSize,
                                         const Real& a_time, 
					 const LevelSigmaCS* a_crseCoords,
					 const int a_refRatio)
{
  // set some boring default values
  int basal_type = constantZb;
  int thickness_type = constantThickness;

  ParmParse geomPP("geometry");

  std::string string_in;
  geomPP.get("basal_type", string_in);
  // there's probably a more clever way to do this...
  if (string_in == "constantZb")
    {
      basal_type = constantZb;
    }
  else if (string_in == "xInclineZb")
    {
      basal_type = xInclineZb;
    }
  else if (string_in == "yInclineZb")
    {
      basal_type = yInclineZb;
    }
  else if (string_in == "sinusoidalYZb")
    {
      basal_type = sinusoidalYZb;
    }
  else if (string_in == "pattynAZb")
    {
      basal_type = pattynAZb;
    }
    else if (string_in == "pattynBZb")
    {
      basal_type = pattynBZb;
    }
  else 
    {
      MayDay::Error("unknown basal type");
    }

  std::string thickness_in;
  geomPP.get("thickness_type", thickness_in);

  if (thickness_in == "constantThickness")
    {
      thickness_type = constantThickness;
    }
  else if (thickness_in == "constantThickness1km")
    {
      thickness_type = constantThickness1km;
    }
  else if (thickness_in == "constantZs1km")
    {
      thickness_type = constantZs1km;
    }
  else if (thickness_in == "doubleZb")
    {
      thickness_type = doubleZb;
    }                   
  else if (thickness_in == "sinusoidalH")
    {
      thickness_type = sinusoidalH;
    }                                      
  else if (thickness_in == "sinusoidalHx")
    {
      thickness_type = sinusoidalHx;
    }                                      
  else if (thickness_in == "singleSinBump")
    {
      thickness_type = singleSinBump;
    }                                      
  else if (thickness_in == "circle")
    {
      thickness_type = circle;
    }  
  else if (thickness_in == "pattynAH")
    {
      thickness_type = pattynAH;
    } 
  else if (thickness_in == "pattynBH")
    {
      thickness_type = pattynBH;
    } 
  else 
    {
      MayDay::Error("unknown thickness type");
    }
                   
  // get scaling factor for ice thickness
  Real thicknessScale=1.0;
  if (thickness_in == "pattynAH" || thickness_in == "pattynBH" )
    {
      thicknessScale = 500.0;
    }
  geomPP.query("thickness_scale", thicknessScale);
  
  RealVect basalSlope(RealVect::Zero);
  
  if (geomPP.contains("basalSlope"))
    {
      Vector<Real> t(SpaceDim, 0.0);
      geomPP.getarr("basalSlope", t, 0, SpaceDim);
       D_TERM(basalSlope[0] = t[0];,
	      basalSlope[1] = t[1];,
	      basalSlope[2] = t[2];);
    }
  

  defineLevelSigmaCS(a_coords,
		     a_domainSize,
		     thickness_type,
		     basal_type,
		     basalSlope,
		     thicknessScale);

}


void 
BasicThicknessIBC::setupBCs()
{
  ParmParse ppBC("bc");

  // get boundary conditions 
  Vector<int> loBCvect(SpaceDim), hiBCvect(SpaceDim);
  ppBC.getarr("lo_bc", loBCvect, 0, SpaceDim);
  ppBC.getarr("hi_bc", hiBCvect, 0, SpaceDim);

  // this is a placeholder until I can get a BCHolder to work...
  // require all directions to have the same BC for now
  CH_assert(loBCvect[0] == loBCvect[1]);
  CH_assert(hiBCvect[0] == hiBCvect[1]);
  CH_assert(loBCvect[0] == hiBCvect[0]);

  if (loBCvect[0] == 0)
    {
      //BCFuncWrapper* newBCPtr = new BCFuncWrapper(iceDirichletBC);
      m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceDirichletBC));
    }
  else if (loBCvect[0] == 1)
    {
      m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceNeumannBC));
    }
  else
    {
      MayDay::Error("bad BC type");
    }
  m_isBCsetUp = true;
}
      


#include "NamespaceFooter.H"
