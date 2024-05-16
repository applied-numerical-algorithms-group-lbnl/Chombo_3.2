#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "HumpIBC.H"
#include "ParmParse.H"
#include "BoxIterator.H"
#include "ExtrapBCF_F.H"
#include "AverageF_F.H"
#include "IceConstants.H"
#include "CONSTANTS.H"

#include "FourthOrderAverage.H"

#include "NamespaceHeader.H"

#define IC_ORDER  2

void zeroBCValueHump(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  // first set all values to zero
  a_values[0] = 0.0;
  a_values[1] = 0.0;
}


void HumpVelBC(FArrayBox& a_state,
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
                  DiriBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValueHump,
                         dir,
                         Side::Lo);
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
                          zeroBCValueHump,
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
                  //        zeroBCValueHump,
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
                  //                zeroBCValueHump,
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
HumpIBC::HumpIBC() 
{
  m_isBCsetUp = false;
  m_paramsSet = false;
  m_isDefined = false;
  m_doHOinit = true;
}

HumpIBC::~HumpIBC()
{
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
HumpIBC::define(const ProblemDomain& a_domain,
                const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

void
HumpIBC::setParameters(Real a_maxThickness,
                       Real a_radSqr,
                       Real a_baseElevation,
                       Real a_minThickness,
                       RealVect a_center,
                       const Real& a_seaLevel,
                       const RealVect& a_widthScale)
{

  m_maxThickness = a_maxThickness/pow(a_radSqr,0.5);
  m_radSqr = a_radSqr;
  m_baseElevation = a_baseElevation;
  m_minThickness = a_minThickness;
  m_center = a_center;
  m_widthScale = a_widthScale;
  m_seaLevel = a_seaLevel;

  m_paramsSet = true;

}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
HumpIBC::new_thicknessIBC()
{
  HumpIBC* retval = new HumpIBC();

  retval->m_maxThickness = m_maxThickness;
  retval->m_radSqr = m_radSqr;
  retval->m_baseElevation = m_baseElevation;
  retval->m_minThickness = m_minThickness;
  retval->m_center = m_center;
  retval->m_widthScale = m_widthScale;
  retval->m_seaLevel = m_seaLevel;
  retval->m_doHOinit = m_doHOinit;

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
HumpIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("HumpIBC::initialize not implemented");
}

/// Set boundary fluxes
/**
 */
void 
HumpIBC::primBC(FArrayBox&            a_WGdnv,
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
HumpIBC::setBdrySlopes(FArrayBox&       a_dW,
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
HumpIBC::artViscBC(FArrayBox&       a_F,
                   const FArrayBox& a_U,
                   const FArrayBox& a_divVel,
                   const int&       a_dir,
                   const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("HumpIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
HumpIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs;
}

/// if appropriate, modify velocity solve RHS in a problem-dependent way. 
/** 
    (not necessary for the hump problem)
*/
void
HumpIBC::modifyVelocityRHS(LevelData<FArrayBox>& a_rhs,
                           LevelSigmaCS& a_coords,
                           const ProblemDomain& a_domain,
                           Real a_time, Real a_dt)
{

}

/// set non-periodic ghost cells for surface height z_s. 
/** 
    use linear extrapolation (corresponds to one-sided differencing)
 */
void
HumpIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
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


  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
      by calling the defineSigmaCS function in the util directory
   */
void
HumpIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                               const RealVect& a_dx,
                               const RealVect& a_domainSize,
                               const Real& a_time, 
			       const LevelSigmaCS* a_crseCoords,
			       const int a_refRatio)
{

  CH_assert(m_paramsSet);
  a_coords.setSeaLevel(m_seaLevel);
      
  const DisjointBoxLayout& grids = a_coords.grids();

#ifdef USEFOURTHORDER
  // reality check to make sure we're doing what we think we're doing..."
  pout() << "using 4th-order cell-averaging for initial thickness profile" << endl;
#endif

  DataIterator dit = grids.dataIterator();

  const RealVect& dx = a_coords.dx();
  const LevelData<FArrayBox>& zBref = a_coords.getTopography();
  LevelData<FArrayBox> zBlocal(grids, 1, zBref.ghostVect());
  LevelData<FArrayBox>& levelH = a_coords.getH();

  // generate sufficiently-accurate initial condition for H by 
  // computing a refined solution and then averaging down.
  int initRef = 8;
  RealVect fineDx = dx/initRef;
  Box refbox(IntVect::Zero,
             (initRef-1)*IntVect::Unit);

  // // to switch back to 2nd order, set factor to be zero...

  // Real factor;
  // if (IC_ORDER == 4)
  //   {
  //     factor = 1.0/24.0;
  //   }
  // else
  //   {
  //     factor = 0.0;
  //   }

  for (dit.begin(); dit.ok(); ++dit)
    { 
      FArrayBox& zB = zBlocal[dit];
      FArrayBox& H = levelH[dit];

      zB.setVal(m_baseElevation);
      
      // allocate storage for refined patch
      Box fineBox(H.box());
      fineBox.refine(initRef);
      FArrayBox fineH(fineBox, 1);

      BoxIterator bit(fineH.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += 0.5*RealVect::Unit;
          x *= fineDx;

          x -= m_center;

#define HUMP2D
#ifdef HUMP2D
          Real radSqr = D_TERM(m_widthScale[0]*x[0]*x[0], 
                               +m_widthScale[1]*x[1]*x[1], 
                               +m_widthScale[2]*x[2]*x[2]);
#else
          Real radSqr = m_widthScale[0]*x[0]*x[0];
#endif


          Real thickness;

          //#define PARABOLIC_PROFILE
#ifdef PARABOLIC_PROFILE
          //  parabolic profile
          if (radSqr <= m_radSqr)
            {
              radSqr = m_radSqr - radSqr;
              radSqr /= (a_domainSize[0]*a_domainSize[1]);
              thickness = m_maxThickness*pow(radSqr, 0.5);
            }
#else
          // cosine profile
          Real phi = 0.5*Pi*(radSqr/m_radSqr);
          if (phi < Pi)
            {
              thickness = max(m_maxThickness*cos(phi),0.0);
            }
          else 
            {
              thickness = 0.0;
            }
#endif
          

          thickness += m_minThickness;
          if (initRef > 1)
            {
              fineH(iv,0) = thickness;
            }
          else
            {
              H(iv,0) = thickness;
            }
        } // end loop over refined thickness box

      if (initRef > 1)
        {
#ifdef USEFOURTHORDER
          // convert from point values to 4th-order cell averages
          fourthOrderAverageCell(fineH);
#endif

          // now average down to fill H
          FORT_AVERAGE(CHF_FRA(H),
                       CHF_CONST_FRA(fineH),
                       CHF_BOX(H.box()),
                       CHF_CONST_INT(initRef),
                       CHF_BOX(refbox));
        }
          
    } // end loop over boxes

  a_coords.setTopography(zBlocal);
  //a_coords.recomputeGeometry();
}

void 
HumpIBC::setupBCs()
{
  m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(HumpVelBC));

  m_isBCsetUp = true;
}
      


#include "NamespaceFooter.H"
