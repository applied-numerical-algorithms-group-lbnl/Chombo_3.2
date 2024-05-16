#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MarineIBC.H"
#include "ParmParse.H"
#include "BoxIterator.H"
#include "ExtrapBCF_F.H"
#include "PetscCompGridVTO.H"
#include "IceConstants.H"
#include "ReflectGhostCells.H"

#include "NamespaceHeader.H"



void zeroBCValueMarine(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  // first set all values to zero
  a_values[0] = 0.0;
  a_values[1] = 0.0;
}


void MarineVelBC(FArrayBox& a_state,
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

              // boundary conditions are Dirichlet on low side,
              // Neumann on high side for dir == 0, and Neumann for dir == 1
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);


              //valid.grow(1);
              //valid.grow(dir,-1);

              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxLo.grow(1);
                  ghostBoxLo.grow(dir,-1);


		  if (dir == 0)
		    {
		      DiriBC(a_state,
			     valid,
			     a_dx,
			     a_homogeneous,
			     zeroBCValueMarine,
			     dir,
			     Side::Lo);
		    } 
		  else 
		    {
		      Interval NeumInterval(0,SpaceDim-1);
		      NeumBC(a_state,
			     valid,
			     a_dx,
			     a_homogeneous,
			     zeroBCValueMarine,
			     dir,
			     Side::Lo,
			     NeumInterval);

		    }
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
                  // want to set corners if possible, so grow ghostBox in the 
                  // transverse directions
                  ghostBoxHi.grow(1);
                  ghostBoxHi.grow(dir,-1);

		  Interval NeumInterval(0,SpaceDim-1);
		  NeumBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          zeroBCValueMarine,
                          dir,
                          Side::Hi,
                          NeumInterval);

                  // // normal-component BC is Neumann, 
                  // // transverse vel BC is Dirichlet
                  // Interval NeumInterval(dir,dir);
                  // NeumBC(a_state,
                  //        valid,
                  //        a_dx,
                  //        a_homogeneous,
                  //        zeroBCValueMarine,
                  //        dir,
                  //        Side::Hi,
                  //        NeumInterval);
                  // for (int comp=0; comp<a_state.nComp(); comp++)
                  //   {
                  //     if (comp != dir)
                  //       {
                  //         Interval DiriInterval(comp,comp);
                  //         DiriBC(a_state,
                  //                valid,
                  //                a_dx,
                  //                a_homogeneous,
                  //                zeroBCValueMarine,
                  //                dir,                                 
                  //                Side::Hi,
                  //                DiriInterval);
                  //       } // end if comp != dir
                  //   } // end loop over components
                }

            } // end if is not periodic in ith direction
        }
    }
}


// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
MarineIBC::MarineIBC() 
{
  m_isBCsetUp = false;
  m_paramsSet = false;
  m_isDefined = false;
}

MarineIBC::~MarineIBC()
{
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
MarineIBC::define(const ProblemDomain& a_domain,
                      const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

void
MarineIBC::setParameters(RefCountedPtr<RealFunction<RealVect > > a_thicknessFunction,
			 RefCountedPtr<RealFunction<RealVect > > a_bedrockFunction,
			 const Real& a_seaLevel)
{
  Vector< RefCountedPtr<RealFunction<RealVect > > > bedrockFuncVect(1, a_bedrockFunction);
  
  setParameters(a_thicknessFunction, bedrockFuncVect, a_seaLevel);

}

void
MarineIBC::setParameters(RefCountedPtr<RealFunction<RealVect > > a_thicknessFunction,
			 Vector<RefCountedPtr<RealFunction<RealVect > > > a_bedrockFunction,
			 const Real& a_seaLevel)
{
  m_bedrockFunction = a_bedrockFunction;
  m_thicknessFunction = a_thicknessFunction;
  m_seaLevel = a_seaLevel;
  m_paramsSet = true;

}


// /// set parameters : inclined plane geomety
// void
// MarineIBC::setParameters(const Real& a_thickness,
// 	    const Real& a_originElevation,
// 	    const RealVect& a_slope,
// 	    const Real& a_seaLevel)
// {

//   if  (m_bedrockFunction != NULL)
//     delete m_bedrockFunction;

//   m_bedrockFunction = new InclinedPlaneFunction(a_originElevation, a_slope);

//   if  (m_thicknessFunction != NULL)
//     delete m_thicknessFunction;

//   m_thicknessFunction = new ConstantRealFunction<RealVect>(a_thickness);


//   m_seaLevel = a_seaLevel;
//   m_paramsSet = true;

// }


// // set parameters : katz geomety
// void
// MarineIBC::setParametersKatz(const Real& a_thickness,
// 			     const Real& a_domainLength, 
// 			     const Real& a_domainWidth, 
// 			     const Real& a_originElevation, 
// 			     const Real& a_alpha, 
// 			     const Real& a_sigma,
// 			     const Real& a_coeff2, 
// 			     const Real& a_coeff4, 
// 			     const Real& a_coeff6,
// 			     const Real& a_seaLevel)
// {

//   if  (m_bedrockFunction != NULL)
//     delete m_bedrockFunction;

//   m_bedrockFunction = new KatzBasalElevation(a_domainLength, a_domainWidth, 
// 					     a_originElevation,  a_alpha,  a_sigma,
// 					     a_coeff2,  a_coeff4,  a_coeff6);
  
//   if  (m_thicknessFunction != NULL)
//     delete m_thicknessFunction;

//   m_thicknessFunction = new ConstantRealFunction<RealVect>(a_thickness);

//   m_seaLevel = a_seaLevel;
//   m_paramsSet = true;

// }


/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
MarineIBC::new_thicknessIBC()
{
  MarineIBC* retval = new MarineIBC();

  retval->m_thicknessFunction = m_thicknessFunction;
  retval->m_bedrockFunction = m_bedrockFunction;
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
MarineIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("MarineIBC::initialize not implemented");
}

/// Set boundary fluxes
/**
 */
void 
MarineIBC::primBC(FArrayBox&            a_WGdnv,
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
MarineIBC::setBdrySlopes(FArrayBox&       a_dW,
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
MarineIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("MarineIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
MarineIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  //  BCFunction* thisBC = dynamic_cast<BCFunction*>(&(*(m_velBCs)));
  return m_velBCs;
}

#if 0
/// if appropriate, modify velocity solve RHS in a problem-dependent way. 
/** 
    add a momemtum source s[dir] = 1/Dx * 0.5 *  ri/rw * (ri - rw) * g * H^2
    to the rhs[dir] at the high dir boundary (1D or 2D cases), or 
    s[dir] = - 1/Dx * ri/rw * (ri - rw) * g * H (3D case) , which is equvalent 
    to setting the momentum flux through the surface (with unit normal n)  
    T.n n = s n

    currently, dir = 0 only

*/
void
MarineIBC::modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
			     LevelSigmaCS& a_coords,
			     const ProblemDomain& a_domain,
			     Real a_time, Real a_dt)
{
  const DisjointBoxLayout& grids = a_coords.grids();
  const LevelData<FArrayBox>& H = a_coords.getH();    

  for (int dir=0; dir<1; dir++)
    {
     
      // incorrect for 3D, for now
      CH_assert(SpaceDim < 3);

      if (!a_domain.isPeriodic(dir))
        {
          // row of cells just inside domain
          Box edgeBox = adjCellHi(a_domain.domainBox(), dir, -1);
          const Real& Dx = a_coords.dx()[dir];
                  
          DataIterator dit = grids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {              
              Box intersectBox(grids[dit]);
              intersectBox &= edgeBox;
              if (!intersectBox.isEmpty())
                {
                  // for now, use cell-centered values
                  // (probably want to revisit this later)
		  
                  const FArrayBox& thisH = H[dit];
		  Real f = 0.5 * a_coords.gravity()*a_coords.iceDensity()/a_coords.waterDensity() * 
		    (a_coords.iceDensity() - a_coords.waterDensity()) / Dx;

		  for (BoxIterator bit(intersectBox); bit.ok(); ++bit){

		    a_rhs[dit](bit(),dir) = f *  std::pow(thisH(bit()),2); 
		  }
                } // end if this box abuts non-periodic high domain boundary
            } // end loop over grids
	} // end if not periodic in this direction  
    } // end loop over directions
}
#endif
/// set non-periodic ghost cells for thickness & topography, 
/** reflection */
void
MarineIBC::setGeometryBCs(LevelSigmaCS& a_coords,
			  const ProblemDomain& a_domain,
			  const RealVect& a_dx,
			  Real a_time, Real a_dt)
{
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
      ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
      ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
      ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
    }

}

/// set non-periodic ghost cells for surface height z_s. 
/** reflection */
void
MarineIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                   LevelSigmaCS& a_coords,
                                   const ProblemDomain& a_domain,
                                   const RealVect& a_dx,
                                   Real a_time, Real a_dt)
{

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
      ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
    }
}

/// set up topography, etc at regrid time
bool
MarineIBC::regridIceGeometry(LevelSigmaCS& a_coords,
			     const RealVect& a_dx,
			     const RealVect& a_domainSize,
			     const Real& a_time,
			     const LevelSigmaCS* a_crseCoords,
			     const int a_refRatio)
 {
   CH_assert(m_paramsSet);
   const DisjointBoxLayout& grids = a_coords.grids();
   DataIterator dit = grids.dataIterator();
   const RealVect& dx = a_coords.dx();
   LevelData<FArrayBox>& levelZb = a_coords.getTopography();
    for (dit.begin(); dit.ok(); ++dit)
    { 
      FArrayBox& zB = levelZb[dit];
      zB.setVal(0.0);
      BoxIterator bit(zB.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;
          for (int i=0; i<m_bedrockFunction.size(); i++)
            {
              zB(iv,0) += (*(m_bedrockFunction[i]))(x);
            }
        } 
    } // end loop over boxes

    return true;
 }


/// set up initial ice state  
void
MarineIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                 const RealVect& a_dx,
                                 const RealVect& a_domainSize,
                                 const Real& a_time,
				 const LevelSigmaCS* a_crseCoords,
				 const int a_refRatio)
{

  CH_assert(m_paramsSet);
  a_coords.setSeaLevel(m_seaLevel);
      
  const DisjointBoxLayout& grids = a_coords.grids();

  DataIterator dit = grids.dataIterator();
  const RealVect& dx = a_coords.dx();

  const LevelData<FArrayBox>& levelZbRef = a_coords.getTopography();
  LevelData<FArrayBox>& levelH = a_coords.getH();
  LevelData<FArrayBox> levelZb(grids, 1, levelZbRef.ghostVect());

  for (dit.begin(); dit.ok(); ++dit)
    { 
      FArrayBox& zB = levelZb[dit];
      FArrayBox& H = levelH[dit];
      
      zB.setVal(0.0);

      BoxIterator bit(zB.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= dx;

          for (int i=0; i<m_bedrockFunction.size(); i++)
            {              
              zB(iv,0) += (*m_bedrockFunction[i])(x);
            }

	   if (SpaceDim == 2){
	     //a little perturbations to encourage instabilities to happen away 
	     //from the edges only meaningful on periodic (in y) domains
	     //zB(iv,0) += 1.0 * sin(2.0 * M_PI * x[1] / a_domainSize[1]);
	   }

	   H(iv,0) = (*m_thicknessFunction)(x);
        } 
              
    } // end loop over boxes

  a_coords.setTopography(levelZb);
  //a_coords.recomputeGeometry();
}

void 
MarineIBC::setupBCs()
{
  ParmParse ppBC("bc");
  Vector<int> loBCvect(SpaceDim);
  Vector<int> hiBCvect(SpaceDim);
  ppBC.getarr("lo_bc", loBCvect, 0, SpaceDim);
  ppBC.getarr("hi_bc", hiBCvect, 0, SpaceDim);
  bool new_bc = false;
  ppBC.query("new_bc", new_bc);

  if (new_bc)
    {
      m_velBCs = RefCountedPtr<BCFunction>(dynamic_cast<BCFunction*>(new CompGridVTOBC));
      // now specifically set bcs (default is dirichlet)
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // first do low-side...
          // dirichlet case
          if (loBCvect[dir] == 0)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][0][comp] = true;
                }
            }
          // neumann case
          else if (loBCvect[dir] == 1)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][0][comp] = false;
                }
            }
          // slip-wall case -- Dirichlet in normal direction, free-slip otherwise
          else if (loBCvect[dir] == 2)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][0][comp] = false;
                }
              m_velBCs->m_bcDiri[dir][0][dir] = true;
            }
          else
            {
              MayDay::Error("Unknown BC type in MarineIBC::setupBCs");
            }

          // then do high side...
          // dirichlet case
          if (hiBCvect[dir] == 0)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][1][comp] = true;
                }
            }
          // neumann case
          else if (hiBCvect[dir] == 1)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][1][comp] = false;
                }
            }
          // slip-wall case -- Dirichlet in normal direction, free-slip otherwise
          else if (hiBCvect[dir] == 2)
            {
              for (int comp=0; comp<SpaceDim; comp++)
                {
                  m_velBCs->m_bcDiri[dir][1][comp] = false;
                }
              m_velBCs->m_bcDiri[dir][1][dir] = true;
            }
          else
            {
              MayDay::Error("Unknown BC type in MarineIBC::setupBCs");
            }

        } // end loop over directions    
    } // end if new BCs
  else
    {
      m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(MarineVelBC));
    }

  m_isBCsetUp = true;
}
      


#include "NamespaceFooter.H"
