#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "VieliPayneIBC.H"
#include "ParmParse.H"
#include "BoxIterator.H"
#include "ExtrapBCF_F.H"
#include "IceConstants.H"

#include "NamespaceHeader.H"

//RealVect VieliPayneIBC::s_edgeThickness(RealVect::Zero);

void zeroBCValueVP(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  // first set all values to zero
  a_values[0] = 0.0;
  a_values[1] = 0.0;

  // move this to RHS of solve
#if 0 
  // if hi-side, compute pressure bc
  if (*side == Side::Hi)
    {
      // value of A from Vieli&Payne (2005)
      // Pa^{-3}a^{-1} -- convert to Pa^{-3}/s
      Real A = 1e-18/SECONDSPERYEAR; 
      Real defaultgravity = 9.81;
      int n = 3;
      // bcVal from Vieli&Payne
      Real bcVal = (0.25*defaultgravity*defaulticedensity*(1.0 - defaulticedensity/defaultseawaterdensity));
      bcVal *= VieliPayneIBC::s_edgeThickness[*dir];
      bcVal = A*pow(bcVal, n);
      a_values[*dir]= bcVal;
    }
#endif

}


void VieliPayneVelBC(FArrayBox& a_state,
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

              // boundary conditions are 
              // Neumann on high side,
              // Dirichlet on low side,for half-domain, 
              // Neumann on low-side for full-domain
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);


              valid.grow(1);
              valid.grow(dir,-1);

              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxLo.grow(1);
                  ghostBoxLo.grow(dir,-1);
                  if (valid.loVect()[dir] < 0)
                    {
                      // full domain case -- we're at a calving front
                      // normal-component BC is Neumann, 
                      // transverse vel BC is Dirichlet
                      Interval NeumInterval(dir,dir);
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             zeroBCValueVP,
                             dir,
                             Side::Lo,
                             NeumInterval);                      
                      for (int comp=0; comp<a_state.nComp(); comp++)
                        {
                          if (comp != dir)
                            {
                              Interval DiriInterval(comp,comp);
                              DiriBC(a_state,
                                     valid,
                                     a_dx,
                                     a_homogeneous,
                                     zeroBCValueVP,
                                     dir,                                 
                                     Side::Lo,
                                     DiriInterval);
                            } // end if comp != dir
                        } // end loop over components
                    }
                  else
                    {
                      // half-domain - we're at ice divide                    
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             zeroBCValueVP,
                             dir,
                             Side::Lo);
                      
                      for (int comp=0; comp<a_state.nComp(); comp++)
                        {
                          if (comp != dir)
                            {
                              Interval DiriInterval(comp,comp);
                              DiriBC(a_state,
                                     valid,
                                     a_dx,
                                     a_homogeneous,
                                     zeroBCValueVP,
                                     dir,                                 
                                     Side::Lo,
                                     DiriInterval);
                            } // end if comp != dir
                        } // end loop over components 
                    }
                }
              
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxHi.grow(1);
                  ghostBoxHi.grow(dir,-1);

                  // normal-component BC is Neumann, 
                  // transverse vel BC is Dirichlet
                  Interval NeumInterval(dir,dir);
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValueVP,
                         dir,
                         Side::Hi,
                         NeumInterval);
                  
                  for (int comp=0; comp<a_state.nComp(); comp++)
                    {
                      if (comp != dir)
                        {
                          Interval DiriInterval(comp,comp);
                          DiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 zeroBCValueVP,
                                 dir,                                 
                                 Side::Hi,
                                 DiriInterval);
                        } // end if comp != dir
                    } // end loop over components
                }

            } // end if is not periodic in ith direction
        }
    }
}


// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
VieliPayneIBC::VieliPayneIBC() 
{
  m_isBCsetUp = false;
  m_paramsSet = false;
  m_isDefined = false;
}

VieliPayneIBC::~VieliPayneIBC()
{
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
VieliPayneIBC::define(const ProblemDomain& a_domain,
                      const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}


/// set parameters
void
VieliPayneIBC::setParameters(const Real& a_thickness,
                             const RealVect& a_slope,
                             const Real& a_originElevation,
                             const Real& a_seaLevel)
{

  m_thickness = a_thickness;
  m_slope = a_slope;
  m_originElevation = a_originElevation;
  m_seaLevel = a_seaLevel;
  m_paramsSet = true;

}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
VieliPayneIBC::new_thicknessIBC()
{
  VieliPayneIBC* retval = new VieliPayneIBC();

  retval->m_thickness = m_thickness;
  retval->m_slope = m_slope;
  retval->m_originElevation = m_originElevation;
  retval->m_seaLevel = m_seaLevel;
  retval->m_paramsSet = true;
  retval->m_isDefined = m_isDefined;
  retval->m_velBCs = m_velBCs;
  retval->m_isBCsetUp = m_isBCsetUp;

  return static_cast<IceThicknessIBC*>(retval);
}

/// Set up initial conditions
/**
 */
void
VieliPayneIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("VieliPayneIBC::initialize not implemented");
}

/// Set boundary fluxes
/**
 */
void 
VieliPayneIBC::primBC(FArrayBox&            a_WGdnv,
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
      
      // (DFM - 5/28/10) this little dance with the ghostBox is a bit 
      // of a kluge to handle the case where a_WGdnv has more than one layer   
      // of ghosting, in which case just testing agains tmp isn't 
      // sufficient to determine whether you're up against the domain edge
      Box ghostBox;
      if (a_side == Side::Lo)
        {
          ghostBox = adjCellLo(m_domain.domainBox(),a_dir, 1);
        }
      else
        {
          ghostBox = adjCellHi(m_domain.domainBox(),a_dir, 1);
        }
      ghostBox &= tmp;

      // Is there a domain boundary next to this grid
      if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
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
	  //          a_WGdnv.setVal(m_thickness, boundaryBox, 0, 1);
	  BoxIterator bit(boundaryBox);
	  for (bit.begin(); bit.ok(); ++bit){
	    const IntVect& i = bit();
	    a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	  }

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
VieliPayneIBC::setBdrySlopes(FArrayBox&       a_dW,
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
VieliPayneIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("VieliPayneIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
VieliPayneIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs;
}


/// fill ghost cells on velocity
void
VieliPayneIBC::velocityGhostBC(LevelData<FArrayBox>& a_velocity,
                               LevelSigmaCS& a_coords,
                               const ProblemDomain& a_domain,
                               Real a_time)
{
  //BCFunc* velSolveBC = velocitySolveBC();

  Real dx = a_coords.dx()[0];
  const DisjointBoxLayout& grids = a_velocity.getBoxes();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      FArrayBox& thisVel = a_velocity[dit];

      bool homogeneous  = false;


      // at some point, we probably want to handle this in a way
      // that doesn't involve so much cut-n-paste reuse...
      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (!a_domain.isPeriodic(dir))
            {
              // lo-side
              Box ghostBoxLow = adjCellLo(gridBox, dir, 1);
              if (!a_domain.contains(ghostBoxLow))
                {
                  // do lo-side BC's (cut and paste from solverBC)
                  ghostBoxLow.grow(1);
                  ghostBoxLow.grow(dir, -1);
                  DiriBC(thisVel,
                         gridBox,
                         dx,
                         homogeneous,
                         zeroBCValueVP,
                         dir,
                         Side::Lo);                  
                }
              Box ghostBoxHi = adjCellHi(gridBox, dir, 1);
              if (!a_domain.contains(ghostBoxHi))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxHi.grow(1);
                  ghostBoxHi.grow(dir,-1);

                  // normal-component high BC is Extrap, for now 

                  Interval ExtrapInterval(dir,dir);
                  ExtrapolateBC(thisVel,
                                gridBox,
                                dx,
                                dir,
                                Side::Hi,
                                ExtrapInterval);
                  
                  // transverse vel BC is Dirichlet
                  for (int comp=0; comp<thisVel.nComp(); comp++)
                    {
                      if (comp != dir)
                        {
                          Interval DiriInterval(comp,comp);
                          DiriBC(thisVel,
                                 gridBox,
                                 dx,
                                 homogeneous,
                                 zeroBCValueVP,
                                 dir,                                 
                                 Side::Hi,
                                 DiriInterval);
                        } // end if comp != dir
                    } // end loop over components

                }
            

            } // end if not periodic 
        } // end loop over directions
    } // end loop over grids

}


/// if appropriate, modify velocity solve RHS in a problem-dependent way. 
/** 
    add a momemtum source s[dir] = 1/Dx * 0.5 *  ri/rw * (ri - rw) * g * H^2
    to the rhs[dir] at the high dir boundary (1D or 2D cases), or 
    s[dir] = - 1/Dx * ri/rw * (ri - rw) * g * H (3D case) , which is equvalent 
    to setting the momentum flux through the surface (with unit normal n)  
    T.n n = s n
*/
void
VieliPayneIBC::modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                                 LevelSigmaCS& a_coords,
                                 const ProblemDomain& a_domain,
                                 Real a_time, Real a_dt)
{
  // this is no longer necessary
#if 0
  const DisjointBoxLayout& grids = a_coords.getBoxes();

  const RealVect& dxVect = a_coords.dx();
  for (int dir=0; dir<SpaceDim; dir++)
    {
     
      // incorrect for 3D, for now
      CH_assert(SpaceDim < 3);

      if (!a_domain.isPeriodic(dir))
        {
          // row of cells just inside domain
          Box edgeBox = adjCellHi(a_domain.domainBox(), dir, -1);
          DataIterator dit = grids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {              
              Box intersectBox(grids[dit]);
              intersectBox &= edgeBox;
              if (!intersectBox.isEmpty())
                {
                  // for now, use cell-centered values
                  // (probably want to revisit this later)
                  Real Dx = dxVect[dir];
		  
                  const FArrayBox& thisH = a_coords[dit].getH();
		  Real f = 0.5 * a_coords.gravity()*a_coords.iceDensity()/a_coords.waterDensity() * 
		    (a_coords.iceDensity() - a_coords.waterDensity()) / Dx;

		  for (BoxIterator bit(intersectBox); bit.ok(); ++bit){

		    a_rhs[dit](bit(),dir) = a_rhs[dit](bit(),dir) 
                      + f *  std::pow(thisH(bit()),2); 
		  }
                } // end if this box abuts non-periodic high domain boundary
            } // end loop over grids
	} // end if not periodic in this direction  
    } // end loop over directions
#endif
}

/// set non-periodic ghost cells for surface height z_s. 
/** 
    use linear extrapolation (corresponds to one-sided differencing)
 */
void
VieliPayneIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx,
                                   Real a_time, Real a_dt)
{
  Box domainBox = a_domain.domainBox();
  const DisjointBoxLayout& grids = a_zSurface.getBoxes();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      FArrayBox& thisZsurf = a_zSurface[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (!a_domain.isPeriodic(dir))
            {
              Box hiBox = adjCellHi(gridBox,dir,1);
              if (!domainBox.contains(hiBox))
                {
                  // set high bc...
                  int hiLo = 1;
                  FORT_SIMPLEREFLECTBC(CHF_FRA(thisZsurf),
                                      CHF_BOX(hiBox),
                                      CHF_INT(dir),
                                      CHF_INT(hiLo));

                }
              Box loBox = adjCellLo(gridBox, dir, 1);
              if (!domainBox.contains(loBox))
                {
                  // set low bc
                  int hiLo = 0;                                    
                  FORT_SIMPLEREFLECTBC(CHF_FRA(thisZsurf),
                                      CHF_BOX(loBox),
                                      CHF_INT(dir),
                                      CHF_INT(hiLo));
                  
                } 
 
            } // end if not periodic in this direction
        } // end loop over directions
    } // end loop over boxes

}


/// set AMR grid hierarchy (for the BC's which need this)
void
VieliPayneIBC::setGridHierarchy(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                                Vector<ProblemDomain>& a_vectDomain)

{
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  //  VieliPayneBCFunction* vpBCptr = dynamic_cast<VieliPayneBCFunction*>(*m_BCfunction);
  m_BCfunction->define(a_vectCS, a_vectDomain);

}


  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
void
VieliPayneIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time, 
				     const LevelSigmaCS* a_crseCoords,
				     const int a_refRatio) 
{

  CH_assert(m_paramsSet);
  const RealVect& dx = a_coords.dx();
  const LevelData<FArrayBox>& zBref = a_coords.getTopography();
  const DisjointBoxLayout& grids = zBref.getBoxes();
  LevelData<FArrayBox> zBlocal(grids, 1, zBref.ghostVect());
  LevelData<FArrayBox>& H = a_coords.getH();
      
  a_coords.setSeaLevel(m_seaLevel);

  DataIterator dit = grids.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    { 
      FArrayBox& thisZbLocal = zBlocal[dit];
      BoxIterator bit(thisZbLocal.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;
          
          Real baseHeight = D_TERM(abs(x[0])*m_slope[0],+abs(x[1])*m_slope[1],+abs(x[2])*m_slope[2]);
          baseHeight += m_originElevation;

          thisZbLocal(iv,0) = baseHeight;
        }
    
      
      // constant thickness
      H[dit].setVal(m_thickness);
      
    } // end loop over boxes

  a_coords.setTopography(zBlocal);
  //a_coords.recomputeGeometry();
}


void 
VieliPayneIBC::setupBCs()
{
  //RefCountedPtr<VieliPayneBCFunction> thisBC(new VieliPayneBCFunction());
  m_BCfunction = RefCountedPtr<VieliPayneBCFunction>(new VieliPayneBCFunction);
  m_velBCs = m_BCfunction;

  m_isBCsetUp = true;
}
      
// ----------------------------------------------------------------------
// VieliPayneBCFunction
// ----------------------------------------------------------------------


VieliPayneBCFunction::VieliPayneBCFunction(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                                           Vector<ProblemDomain>& a_vectDomain)
{
  define(a_vectCS, a_vectDomain);
}

  
VieliPayneBCFunction::~VieliPayneBCFunction()
{
}

void
VieliPayneBCFunction:: define(Vector<RefCountedPtr<LevelSigmaCS > >& a_vectCS,
                              Vector<ProblemDomain>& a_vectDomain)
{
  m_vectCS = a_vectCS;
  m_vectDomain = a_vectDomain;
}
  

void 
VieliPayneBCFunction::operator()(FArrayBox&           a_state,
                                 const Box&           a_valid,
                                 const ProblemDomain& a_domain,
                                 Real                 a_dx,
                                 bool                 a_homogeneous)
{
  MayDay::Error("VieliPayneBCFunction requires a DataIndex");
}


void 
VieliPayneBCFunction::operator()(FArrayBox&           a_state,
                                 const Box&           a_valid,
                                 const ProblemDomain& a_domain,
                                 Real                 a_dx,
                                 const DataIndex&     a_index,
                                 bool                 a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(dir))
            {
              
              // boundary conditions are Dirichlet on low side,
              // Neumann on high side,
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);


              valid.grow(1);
              valid.grow(dir,-1);

              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxLo.grow(1);
                  ghostBoxLo.grow(dir,-1);
                  if (valid.loVect()[dir] < 0)
                    {
                      // full domain case -- we're at a calving front
                      // normal-component BC is Neumann,
                      // transverse vel BC is Dirichlet
                      // normal component BC is marine boundary condition.
                      // if homogeneous, just use homogeneous Neumann BCFunc
                      if (a_homogeneous)
                        {                      
                          // normal-component BC is Neumann, 
                          // transverse vel BC is Dirichlet
                          Interval NeumInterval(dir,dir);
                          NeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 zeroBCValueVP,
                                 dir,
                                 Side::Lo,
                                 NeumInterval);
                        }
                      else
                        {
                          // this is where we actually do something interesting
                          // and apply the inhomogeneous form of the marine bc
                          
                          // first, figure out which level we're on
                          int thisLev=-1; 
                          for (int lev=0; lev<m_vectDomain.size(); lev++)
                            {
                              if (a_domain == m_vectDomain[lev])
                                {
                                  thisLev = lev;
                                }
                            }
                          // check to be sure we actually found a match
                          if (thisLev == -1)
                            { 
                              MayDay::Error("unable to determine AMR level for inhomogeneous Marine BC");
                            }
                          
                          // now grab the SigmaCS for this level
                          const LevelSigmaCS& thisCS = (*m_vectCS[thisLev]);
                          
                          // get face-centered ice thickness 
                          const FluxBox& thisFaceH = thisCS.getFaceH()[a_index];
                          const FArrayBox& thisFaceHDir = thisFaceH[dir];
                      
                          Real dx = thisCS.dx()[dir];
                          
                          //                      Real A = 1e-18/SECONDSPERYEAR;  
                          //Real A = 1e-18;
                          Real A = 9.2e-18;
                          int n = 3;
                          // bcVal from Vieli&Payne
                          Real f = (0.25*thisCS.gravity()*thisCS.iceDensity()
                                    *(1.0 - thisCS.iceDensity()/thisCS.waterDensity()));
                          // negative shiftIV for low side
                          IntVect shiftIV = -BASISV(dir);
                          
                          BoxIterator bit(ghostBoxLo);
                          for (bit.begin(); bit.ok(); ++bit)
                            {
                              IntVect iv = bit();
                              // H from face on the boundary, to the right of ghost cell
                              Real bcVal = thisFaceHDir(iv-shiftIV,0);
                              bcVal *= f;
                              bcVal = A*pow(bcVal,n);
                              a_state(iv,dir) = a_state(iv-shiftIV,dir) - dx*bcVal;
                            }    
                        }                    
                      
                      for (int comp=0; comp<a_state.nComp(); comp++)
                        {
                          if (comp != dir)
                            {
                              int order = 1;
                              Interval DiriInterval(comp,comp);
                              DiriBC(a_state,
                                     valid,
                                     a_dx,
                                     a_homogeneous,
                                     zeroBCValueVP,
                                     dir,                                 
                                     Side::Lo,
                                     DiriInterval,
                                     order);
                            } // end if comp != dir
                        } // end loop over components
                    }
                  else
                    {
                      // at ice divide
                      int order = 2;
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             zeroBCValueVP,
                             dir,
                             Side::Lo,
                             order);
                    }
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxHi.grow(1);
                  ghostBoxHi.grow(dir,-1);

                  // normal component BC is marine boundary condition.
                  // if homogeneous, just use homogeneous Neumann BCFunc
                  if (a_homogeneous)
                    {                      
                      // normal-component BC is Neumann, 
                      // transverse vel BC is Dirichlet
                      Interval NeumInterval(dir,dir);
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             zeroBCValueVP,
                             dir,
                             Side::Hi,
                             NeumInterval);
                    }
                  else
                    {
                      // this is where we actually do something interesting
                      // and apply the inhomogeneous form of the marine bc

                      // first, figure out which level we're on
                      int thisLev=-1; 
                      for (int lev=0; lev<m_vectDomain.size(); lev++)
                        {
                          if (a_domain == m_vectDomain[lev])
                            {
                              thisLev = lev;
                            }
                        }
                      // check to be sure we actually found a match
                      if (thisLev == -1)
                        { 
                          MayDay::Error("unable to determine AMR level for inhomogeneous Marine BC");
                        }
                      
                      // now grab the SigmaCS for this level
                      const LevelSigmaCS& thisCS = (*m_vectCS[thisLev]);
                      
                      // get face-centered ice thickness 
                      const FluxBox& thisFaceH = thisCS.getFaceH()[a_index];
                      const FArrayBox& thisFaceHDir = thisFaceH[dir];
                      
                      Real dx = thisCS.dx()[dir];
                      
                      //                      Real A = 1e-18/SECONDSPERYEAR;  
                      //Real A = 1e-18;
                      Real A = 9.2e-18;
                      int n = 3;
                      // bcVal from Vieli&Payne
                      Real f = (0.25*thisCS.gravity()*thisCS.iceDensity()
				*(1.0 - thisCS.iceDensity()/thisCS.waterDensity()));
                      IntVect shiftIV = BASISV(dir);
                      
                      BoxIterator bit(ghostBoxHi);
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv = bit();
                          Real bcVal = thisFaceHDir(iv,0);
                          bcVal *= f;
                          bcVal = A*pow(bcVal,n);
                          a_state(iv,dir) = a_state(iv-shiftIV,dir) + dx*bcVal;
                        }    
                    }


                  for (int comp=0; comp<a_state.nComp(); comp++)
                    {
                      if (comp != dir)
                        {
                          int order = 1;
                          Interval DiriInterval(comp,comp);
                          DiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 zeroBCValueVP,
                                 dir,                                 
                                 Side::Hi,
                                 DiriInterval,
                                 order);
                        } // end if comp != dir
                    } // end loop over components
                }

            } // end if is not periodic in ith direction
        }
    }



}

#include "NamespaceFooter.H"
