#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//#include "NewCoordSys.H"
#include "LevelSigmaCS.H"
#include "SigmaCSF_F.H"
#include "DerivativesF_F.H"
#include "AMRPoissonOpF_F.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "Averaging.H"
#include "BoxIterator.H"
#include "NodeFArrayBox.H"
#include "NodeAMRIO.H"
#include "LevelData.H"
#include "IceConstants.H"
#include "FineInterp.H"
#include "FineInterpFace.H"
#include "IntFineInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "NamespaceHeader.H"

/// default constructor
LevelSigmaCS::LevelSigmaCS()
{
  setDefaultValues();
  m_isDefined = false;
}

/// defining constructor
LevelSigmaCS::LevelSigmaCS(const DisjointBoxLayout& a_grids,
                           const RealVect& a_dx,
                           const IntVect& a_ghostVect)
{  
  setDefaultValues();
  // this just calls the define function
  define(a_grids, a_dx, a_ghostVect);
}

/**
   Destructor.
*/
LevelSigmaCS::~LevelSigmaCS()
{
}

void
LevelSigmaCS::define(const DisjointBoxLayout& a_grids,
                     const RealVect& a_dx,
                     const IntVect& a_ghostVect)
{
  m_grids = a_grids;
  m_dx = a_dx;
  m_ghostVect = a_ghostVect;

  // now define local storage
  // note that we define H and topography on a one-cell-deep
  // box in the z-direction
  // note that we grow topography and H by one in the horizontal
  // directions so we can define derivatives on a_ghostVect
  // in 2d, however, we just do x,y

  // define H grids
  if (SpaceDim == 1)
    {
      m_Hgrids = m_grids;
    }
  else if (SpaceDim == 2)
    {
      m_Hgrids = m_grids;
    }
  else if (SpaceDim == 3)
    {
      // defer this until we're actually doing 3d
      MayDay::Error("3D LevelSigmaCS not implemented yet");
#if 0
      Box HBox(a_box);
      HBox.grow(1);
      if (SpaceDim == 3)
        {
          HBox.setSmall(0,a_box.smallEnd(0));
          HBox.setBig(0,a_box.smallEnd(0));  
        }
#endif
    }
  
  IntVect Hghost = a_ghostVect;
  // Hghost grown by one in the horizontal directions
  Hghost += IntVect::Unit;
  if (SpaceDim == 3)
    {
      Hghost[0] -= 1;
    }
      
  m_topography.define(m_Hgrids, 1, Hghost);
  m_surface.define(m_Hgrids, 1, 2 * IntVect::Unit);
  m_gradSurface.define(m_Hgrids, SpaceDim, IntVect::Unit);
  m_gradSurfaceFace.define(m_Hgrids, SpaceDim, IntVect::Unit);
  m_thicknessOverFlotation.define(m_Hgrids, 1, IntVect::Unit) ;

  m_H.define(m_Hgrids, 1, Hghost);
  m_faceH.define(m_Hgrids, 1, Hghost);

  //CH_assert((SpaceDim == 2) | (SpaceDim ==3));
  m_deltaFactors.define(a_grids, 2, a_ghostVect);
  
  m_floatingMask.define(m_Hgrids, 1, Hghost);
  m_anyFloating.define(m_Hgrids);

  m_faceDeltaFactors.define(a_grids, 2, a_ghostVect);

  m_isDefined = true;
}

/// define as a coarsening of fineCS by nRef
void
LevelSigmaCS::define(const LevelSigmaCS& a_fineCS, 
                     int a_nRef)
{

  // coarsen grids
  DisjointBoxLayout crseGrids;
  coarsen(crseGrids, m_grids, a_nRef);

  // use same ghosting as a_fineCS  
  IntVect ghostVect = m_deltaFactors.ghostVect();
  RealVect crseDx = a_fineCS.m_dx;
  crseDx *= a_nRef;

  m_seaLevel = a_fineCS.seaLevel();

  // now call the other define function to set up storage
  define(crseGrids, crseDx, ghostVect);

  // now need to average-down data
  horizontalAverage(m_topography, a_fineCS.m_topography, a_nRef);

  horizontalAverage(m_gradSurface, a_fineCS.m_gradSurface,  a_nRef);

  horizontalAverageFace(m_gradSurfaceFace, a_fineCS.m_gradSurfaceFace, 
		    a_nRef);

  horizontalAverage(m_H, a_fineCS.m_H, a_nRef);

  horizontalAverageFace(m_faceH, a_fineCS.m_faceH, a_nRef);

  averageAllDim(m_deltaFactors, a_fineCS.m_deltaFactors, a_nRef);
  
  averageAllDimFace(m_faceDeltaFactors, a_fineCS.m_faceDeltaFactors, 
                    a_nRef);

  
#if CH_SPACEDIM == 2
  m_faceSigma = a_fineCS.m_faceSigma;
  m_sigma = a_fineCS.m_sigma;
  m_dSigma = a_fineCS.m_dSigma;
#endif

}

/// define as a coarsening of fineCS 
LevelSigmaCS::LevelSigmaCS( const DisjointBoxLayout& a_grids,
			    const RealVect& a_dx,
			    const LevelSigmaCS& a_fineCS, 
			    int a_nRef
			    )
{
  // now call the other define function to set up storage
  setDefaultValues();
  define( a_grids, a_dx, a_fineCS.ghostVect() );
 
  DisjointBoxLayout crseGrids;
  coarsen( crseGrids, a_fineCS.m_grids, a_nRef );

  // now need to average-down data
  Copier copier( crseGrids, a_grids );

  /// cell-centered topography
  LevelData<FArrayBox> crs_ldfTop( crseGrids, m_topography.nComp(), IntVect::Zero );
  horizontalAverage( crs_ldfTop, a_fineCS.m_topography, a_nRef );
  crs_ldfTop.copyTo( crs_ldfTop.interval(), m_topography, m_topography.interval(), copier );
  m_topography.exchange();

  /// cell-centered gradient of surface elevation
  LevelData<FArrayBox> crs_gradSurface( crseGrids, m_gradSurface.nComp(), IntVect::Zero );
  horizontalAverage( crs_gradSurface, a_fineCS.m_gradSurface,  a_nRef );
  crs_gradSurface.copyTo( crs_gradSurface.interval(), m_gradSurface, m_gradSurface.interval(), copier);
  m_gradSurface.exchange();

  /// cell-centered ice thickness
  LevelData<FArrayBox> crs_ldfH( crseGrids, m_H.nComp(), IntVect::Zero );
  horizontalAverage( crs_ldfH, a_fineCS.m_H, a_nRef );
  crs_ldfH.copyTo( crs_ldfH.interval(), m_H, m_H.interval(), copier);
  m_H.exchange();

  /// cell-centered surface elevation
  LevelData<FArrayBox> crs_ldfsurface( crseGrids, m_surface.nComp(), IntVect::Zero );
  horizontalAverage( crs_ldfsurface, a_fineCS.m_surface, a_nRef );
  crs_ldfsurface.copyTo( crs_ldfsurface.interval(), m_surface, m_surface.interval(), copier );
  m_surface.exchange();

  /// cell-centered 
  LevelData<FArrayBox> crs_deltaFactors( crseGrids, m_deltaFactors.nComp(), IntVect::Zero );
  averageAllDim(crs_deltaFactors, a_fineCS.m_deltaFactors, a_nRef);
  crs_deltaFactors.copyTo( crs_deltaFactors.interval(), m_deltaFactors, m_deltaFactors.interval(), copier);
  m_deltaFactors.exchange();

  /// cell-centered 
  LevelData<FArrayBox> crs_thicknessOverFlotation( crseGrids, m_thicknessOverFlotation.nComp(), IntVect::Zero );
  //horizontalAverage( crs_thicknessOverFlotation, a_fineCS.m_thicknessOverFlotation, a_nRef );
  averageAllDim( crs_thicknessOverFlotation, a_fineCS.m_thicknessOverFlotation, a_nRef );
  crs_thicknessOverFlotation.copyTo(crs_thicknessOverFlotation.interval(), m_thicknessOverFlotation, m_thicknessOverFlotation.interval(), copier);
  m_thicknessOverFlotation.exchange();

  /// face-centered ice thickness
  CellToEdge( m_H, m_faceH );
  m_faceH.exchange();

  /// face-centered gradient of surface elevation
  CellToEdge( m_gradSurface, m_gradSurfaceFace );
  m_gradSurfaceFace.exchange();

  /// face-centered
  CellToEdge( m_deltaFactors, m_faceDeltaFactors );
  m_faceDeltaFactors.exchange();

  // simple copy stuff
  m_seaLevel = a_fineCS.seaLevel();
#if CH_SPACEDIM == 2
  m_faceSigma = a_fineCS.m_faceSigma;
  m_sigma = a_fineCS.m_sigma;
  m_dSigma = a_fineCS.m_dSigma;
#endif
}

/// given coordinate in mapped space, return its location in real
/// space -- this will be a bit slow; probably want to use Fab-based
/// one instead (once it's implemented in Fortran)
/// note that a_Xi is (sigma,x,y) while realCoord returns (x,y,z)
RealVect
LevelSigmaCS::realCoord(const RealVect& a_Xi, const DataIndex& a_index) const
{
  // note that the vertical direction is the first coordinate 
  // in computational space in 3d, but 2d is regular (x,y)
  // first need to find H for this location -- note that we use the
  // 2D location
#if CH_SPACEDIM == 1
  RealVect loc = a_Xi;
#elif CH_SPACEDIM == 2
  //RealVect loc(a_Xi[1],z);
  RealVect loc = a_Xi;
#elif CH_SPACEDIM == 3
  IntVect iv(D_DECL(a_Xi[0]/m_dx[0],
                    a_Xi[1]/m_dx[1], a_Xi[2]/m_dx[2]));
  if (SpaceDim == 3) iv[0] = m_H[a_index].box().smallEnd(0);
  CH_assert(m_topography[a_index].box().contains(iv));
  Real H = m_H[a_index](iv,0);
  Real zB = m_topography[a_index](iv,0);

  Real z = zB + (1.0 - a_Xi[0])*H;
  RealVect loc(a_Xi[1],a_Xi[2], z);
#else
  MayDay::Error("SigmaCoordSys::RealLoc not defined for dim = CH_SPACEDIM");
#endif

  return loc;

}


/// given coordinate in real space, return its location in the mapped space
RealVect LevelSigmaCS::mappedCoord(const RealVect& a_x, 
                                   const DataIndex& a_index) const
{
  // locate where we are in 2D coordinates

  // first need to find H for this location -- note that we use the
  // 2D location
  RealVect loc;
  if (SpaceDim == 3)
    {
      IntVect iv(D_DECL(m_H[a_index].box().smallEnd(0),
                        (int)(a_x[1]/m_dx[1]), a_x[2]/m_dx[2]));
            
      CH_assert(m_topography[a_index].box().contains(iv));
      Real H = m_H[a_index](iv,0);
      Real zB = m_topography[a_index](iv,0);
      Real sigma = (zB + H - a_x[0])/H;
      loc = RealVect(D_DECL(sigma,a_x[0], a_x[1]));
    }
  else if (SpaceDim == 2)
    {
      loc = a_x;
    }
  return loc;
  
}

/// return Cartesian XYZ locations of nodes 
/** nodeCoords should have dimension returned by dimension() 
 */
void
LevelSigmaCS::getNodeRealCoordinates(LevelData<NodeFArrayBox>& a_nodeCoords) const
{
  CH_assert (a_nodeCoords.getBoxes().compatible(m_grids));

  DataIterator dit = a_nodeCoords.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      NodeFArrayBox& thisNodeFab = a_nodeCoords[dit];
      FArrayBox& thisNodeCoords = thisNodeFab.getFab();
      const Box thisBox = thisNodeCoords.box();
      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          // node centering
          RealVect X = m_dx * iv;
          
          // default implementation just calls RealCoord
          RealVect nodeCoord = realCoord(X, dit());
          D_TERM6(thisNodeCoords(iv,0) = nodeCoord[0];,
                  thisNodeCoords(iv,1) = nodeCoord[1];,
                  thisNodeCoords(iv,2) = nodeCoord[2];,
                  thisNodeCoords(iv,3) = nodeCoord[3];,
                  thisNodeCoords(iv,4) = nodeCoord[4];,
                  thisNodeCoords(iv,5) = nodeCoord[5];)
            } // end loop over nodes
    } // end loop over grids
}

/// returns modifiable cell-centered H (ice sheet thickness) for this
LevelData<FArrayBox>& 
LevelSigmaCS::getH()
{
  
  return m_H;
}

/// returns const reference to cell-centered H (ice sheet thickness) 
const LevelData<FArrayBox>&
LevelSigmaCS::getH() const
{
  return m_H;
}

/// returns modifiable cell-centered H (ice sheet thickness) for this
LevelData<FluxBox>& 
LevelSigmaCS::getFaceH()
{
  return m_faceH;
}


/// returns const reference to face-centered H (ice sheet thickness) 
const LevelData<FluxBox>&
LevelSigmaCS::getFaceH() const
{
 
  return m_faceH;
}

/// returns a const-reference to the cell-centered topography 
const LevelData<FArrayBox>& 
LevelSigmaCS::getTopography() const
{
  return m_topography;
}

/// returns a modifiable reference to the cell-centered topography 
LevelData<FArrayBox>& 
LevelSigmaCS::getTopography() 
{
  
  return m_topography;
}



///sets the base height. 
/** In practice, this will probably done by a derived class.*/
void
LevelSigmaCS::setTopography(const LevelData<FArrayBox>& a_topography)
{
  
  if (a_topography.getBoxes().compatible(m_topography.getBoxes()))
    {
      DataIterator dit = a_topography.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_topography[dit].copy(a_topography[dit], 
                                 m_topography[dit].box());
        }
    } 
  else
    {
      a_topography.copyTo(m_topography);
    }

  

}


///sets the cell-centered surface height
void
LevelSigmaCS::setSurfaceHeight(const LevelData<FArrayBox>& a_surfaceHeight)
{
  // first set to zero
  DataIterator dit = m_surface.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_surface[dit].setVal(0.0);
    }

  // do fab-by-fab copy if possible
  if (a_surfaceHeight.getBoxes().compatible(m_surface.getBoxes()))
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_surface[dit].copy(a_surfaceHeight[dit], a_surfaceHeight[dit].box() );
        }
    } // end if compatible
  else 
    {
      // otherwise, just do a copyTo
      a_surfaceHeight.copyTo(m_surface);
    }
}

///sets the cell-centered surface gradient
void
LevelSigmaCS::setGradSurface(const LevelData<FArrayBox>& a_gradSurface)
{
  // first set to zero
  DataIterator dit = m_gradSurface.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_gradSurface[dit].setVal(0.0);
    }

  // do fab-by-fab copy if possible
  if (a_gradSurface.getBoxes().compatible(m_gradSurface.getBoxes()))
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_gradSurface[dit].copy(a_gradSurface[dit], m_gradSurface[dit].box() );
        }
    } // end if compatible
  else 
    {
      // otherwise, just do a copyTo
      a_gradSurface.copyTo(m_gradSurface);
    }
}

///sets the face-centred surface gradient
void
LevelSigmaCS::setGradSurfaceFace(const LevelData<FluxBox>& a_gradSurfaceFace)
{
  // first, set to zero
  DataIterator dit = m_gradSurfaceFace.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_gradSurfaceFace[dit].setVal(0.0);
    }
  
  if (m_gradSurfaceFace.getBoxes().compatible(a_gradSurfaceFace.getBoxes()))
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir = 0; dir < SpaceDim; ++dir)
            {
              m_gradSurfaceFace[dit][dir].copy(a_gradSurfaceFace[dit][dir], 
                                               m_gradSurfaceFace[dit][dir].box());
            } // end loop over dir
        } // end loop over grids
    } // end if compatible
  else
    {
      // otherwise, just do a copyto
      a_gradSurfaceFace.copyTo(m_gradSurfaceFace);
    }
  
}

///returns the cell-centred surface elevation
const LevelData<FArrayBox>& 
LevelSigmaCS::getSurfaceHeight() const
{
  return m_surface;
}

///returns the cell-centred surface gradient
const LevelData<FArrayBox>& 
LevelSigmaCS::getGradSurface() const
{
  return m_gradSurface;
}

///returns the face-centred surface gradient
const LevelData<FluxBox>& 
LevelSigmaCS::getGradSurfaceFace() const
{
  return m_gradSurfaceFace;
}

const LevelData<FArrayBox>&
LevelSigmaCS::getThicknessOverFlotation() const
{
  return m_thicknessOverFlotation;
}


/// returns ice surface height
/** 
    if ice is grounded, then this equals z_base + H
    if ice is floating, then this equals H*(1 - rho_ice/rho_water)
*/
void
LevelSigmaCS::getSurfaceHeight(LevelData<FArrayBox>& a_zSurface) const
{
  CH_assert(a_zSurface.getBoxes().compatible(m_grids));
  CH_TIME("LevelSigmaCS::getSurfaceHeight");

  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      
      FArrayBox& thisZsurf = a_zSurface[dit];
      const Box& thisBox = thisZsurf.box();
      CH_assert(m_H[dit].box().contains(thisBox));
      
      FORT_SURFACEHEIGHT(CHF_FRA1(thisZsurf,0),
                         CHF_CONST_FRA1(m_H[dit],0),
                         CHF_CONST_FRA1(m_topography[dit],0),
                         CHF_CONST_REAL(m_iceDensity),
                         CHF_CONST_REAL(m_waterDensity),
                         CHF_CONST_REAL(m_seaLevel),
                         CHF_BOX(thisBox));
    } //end loop over grids
                              
}


/// recomputes quantities which are dependent on cell-centered H. 
/** Should be called after H is modified to ensure that things 
    don't get out of sync. */
void
LevelSigmaCS::recomputeGeometry(const LevelSigmaCS* a_crseCoords, 
				const int a_refRatio)
{
  CH_TIME("LevelSigmaCS::recomputeGeometry");
  CH_assert(m_isDefined);
  
  // compute face-averaged Him 
  m_H.exchange();
  m_topography.exchange();
  CellToEdge(m_H, m_faceH);    

  // probably will eventually want to do this more efficiently. For
  // now, go in and set faces adjacent to cells with zero H to zero
  // as well. Also, note that in 3d, the vertical faces are meaningless
  // (eventually go back and switch to 2 FArrayBox's vs. a FluxBox

  int xDir = 0;
  if (SpaceDim == 3) xDir = 1;
  int yDir = xDir +1;
  if (SpaceDim == 1) yDir = 0;
  
  DataIterator dit= m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FORT_FIXFACEH(CHF_FRA1(m_faceH[dit][xDir],0),
                    CHF_FRA1(m_faceH[dit][yDir],0),
                    CHF_CONST_FRA1(m_H[dit],0),
                    CHF_BOX(m_faceH[dit].box()),
                    CHF_INT(xDir),
                    CHF_INT(yDir));
      
    }

  if (a_crseCoords)
    {
      //faceH needs to be zero on the fine side interface if it is on the coarse side. 
      //If not, we end up with some open sea / open land cells on the coarse side
      //connected to grounded ice / floating ice cells on the fine side.
      //The catch is, we don't want to tinker with any faceH not at a coarse-fine
      //interface...

      //first off, need to decide where the coarse-fine interfaces are. This
      //is a pretty stupid way, I suppose, but I was in a hurry, and it works.
      LevelData<FArrayBox> tmp(m_grids, 1, IntVect::Unit);
      for (dit.begin(); dit.ok(); ++dit)
	{
	  tmp[dit].setVal(0.0);
	  tmp[dit].setVal(1.0,m_grids[dit],0,1);
	}
      tmp.exchange();
      //et voila, the coarse-fine interfaces have a 0 in the ghost cell adjacent to
      //them, and all other faces have a 1

      FineInterpFace faceInterp(m_grids, 1, a_refRatio , m_grids.physDomain());
      LevelData<FluxBox> tmpFaceH; tmpFaceH.define(m_faceH);
      faceInterp.interpToFine(tmpFaceH, a_crseCoords->m_faceH);
      
      Real tol = 1.0e-18;
      for (dit.begin(); dit.ok(); ++dit)
	{
	  const Box& b = m_grids[dit];
	  
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      const IntVect& lolo = b.smallEnd();
	      const IntVect& hihi = b.bigEnd();

	      IntVect lohi = hihi;
	      lohi[dir] = lolo[dir];
	      Box loBox(lolo,lohi);
	      for (BoxIterator bit(loBox);bit.ok();++bit)
		{
		   const IntVect& iv = bit();
		   const IntVect ivm = iv - BASISV(dir);
		   if (tmp[dit](ivm) < tol && tmpFaceH[dit][dir](iv) < tol)
		     m_faceH[dit][dir](iv) = 0.0;
		}

	      IntVect hilo = lolo;
	      hilo[dir] = hihi[dir];
	      Box hiBox(hilo,hihi);
	      for (BoxIterator bit(hiBox);bit.ok();++bit)
		{
		   const IntVect& iv = bit();
		   const IntVect ivp = iv + BASISV(dir);
		   if (tmp[dit](ivp) < tol && tmpFaceH[dit][dir](ivp) < tol)
		     m_faceH[dit][dir](ivp) = 0.0;
		}
	    } // end loop over direction
	} // end loop over boxes
    } // end if (a_crseCoords)


  computeDeltaFactors();
  computeSurface(a_crseCoords,  a_refRatio);


  
  
}


/// recomputes quantities which are dependent on cell-centered H. 
/** Should be called after H is modified to ensure that things 
    don't get out of sync. */
void
LevelSigmaCS::recomputeGeometryFace(const LevelSigmaCS* a_crseCoords, 
				    const int a_refRatio)
{
  CH_assert(m_isDefined);
 CH_TIME("LevelSigmaCS::recomputeGeometryFace");
  // compute cell-averaged H
  // H = 0.5*(x-average + y_average)
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox temp(m_H[dit].box(), 1);
      m_H[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          EdgeToCell(m_faceH[dit], 0, temp, 0, dir);    
          temp *= 0.5;
          m_H[dit].plus(temp);
        }
    }
  
  // probably will eventually want to do this more efficiently. For
  // now, go in and set faces adjacent to cells with zero H to zero
  // as well. Also, note that in 3d, the vertical faces are meaningless
  // (eventually go back and switch to 2 FArrayBox's vs. a FluxBox

  int xDir = 0;
  if (SpaceDim == 3) xDir = 1;
  int yDir = xDir +1;
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      FORT_FIXFACEH(CHF_FRA1(m_faceH[dit][xDir],0),
                    CHF_FRA1(m_faceH[dit][yDir],0),
                    CHF_CONST_FRA1(m_H[dit],0),
                    CHF_BOX(m_faceH[dit].box()),
                    CHF_INT(xDir),
                    CHF_INT(yDir));
    }

  computeDeltaFactors();
  //computeFloatingMask();
  computeSurface(a_crseCoords,  a_refRatio);
}


void
LevelSigmaCS::setDefaultValues()
{
  m_seaLevel = 0.0;
  /// densities of ice and seawater, in kg/(m^3)
  /// from Pattyn (2003)
  m_iceDensity = 910.0;
  /// from Vieli and Payne (2005)
  m_waterDensity = 1028.0;

  m_gravity = 9.81;

  m_backgroundSlope = RealVect::Zero;

#if CH_SPACEDIM == 2
  //default is a single layer, which makes sense only for isothermal SSA
  m_faceSigma.resize(2);
  m_faceSigma[0] = 0.0; // top of the ice sheet
  m_faceSigma[0] = 1.0; // bottom of the ice sheet
  m_sigma.resize(1);
  m_sigma[0] = 0.5;
  m_dSigma.resize(1);
  m_dSigma[0] = 1.0;
#endif

}

void
LevelSigmaCS::computeDeltaFactors()
{
 
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {      
      // cell-centered DeltaFactors
      FArrayBox& thisH = m_H[dit];
      FArrayBox& thisTopography = m_topography[dit];
      FArrayBox& thisDeltaFactors = m_deltaFactors[dit];
      Box thisBox = thisH.box();
      if (SpaceDim == 2)
        {
          thisBox.grow(-1);
        }
      else if (SpaceDim == 3)
        {
          IntVect growVect = -IntVect::Unit;
          growVect[0] = 0;
          thisBox.grow(growVect);
        }

      Box derivBox(thisH.box());
      IntVect shrinkVect(D_DECL(0,-1,-1));
      derivBox.grow(shrinkVect);
      
      FArrayBox deltaH(derivBox,1);
      FArrayBox deltaZb(derivBox,1);
  
      for (int derivDir=1; derivDir<SpaceDim; derivDir++)
        {
          // use cc-derivative fortran from DerivativesF.ChF
          FORT_CCDERIV(CHF_FRA1(deltaH,0),
                       CHF_CONST_FRA1(thisH,0),
                       CHF_BOX(derivBox),
                       CHF_CONST_REAL(m_dx[derivDir]),
                       CHF_INT(derivDir));
          
          FORT_CCDERIV(CHF_FRA1(deltaZb,0),
                       CHF_CONST_FRA1(thisTopography,0),
                       CHF_BOX(derivBox),
                       CHF_CONST_REAL(m_dx[derivDir]),
                       CHF_INT(derivDir));
          
          
          // define cell-centered deltas
          FORT_DEFINECELLGEOM(CHF_FRA1(thisDeltaFactors,derivDir-1),
                              CHF_FRA1(deltaH,0),
                              CHF_FRA1(deltaZb,0),
                              CHF_REALVECT(m_dx),
                              CHF_BOX(thisBox),
                              CHF_INT(derivDir));
        } // end loop over derivative directions

    } // end loop over grids

} // end context for cell-centered deltaFactors

void
LevelSigmaCS::computeSurface(const LevelSigmaCS* a_crseCoords, 
			     const int a_refRatio)
{
  
 CH_TIME("LevelSigmaCS::computeSurface");
  const LevelSigmaCS* crseCoords = a_crseCoords;
  const int& refRatio = a_refRatio;
  
  //first off, the surface elevation
  getSurfaceHeight(m_surface);

  if (false && a_crseCoords != NULL)
    {
      
      // fill ghost regions of m_surface
      int nGhost = m_surface.ghostVect()[0];
      PiecewiseLinearFillPatch surfaceFiller
	(m_grids, crseCoords->m_grids , 1, 
	 crseCoords->m_grids.physDomain() ,
	 refRatio, nGhost);
      Real time_interp_coeff = 0.0; //we don't need this
      surfaceFiller.fillInterp
	(m_surface,  crseCoords->m_surface, crseCoords->m_surface,
	 time_interp_coeff, 0, 0, 1);
      
    }
  m_surface.exchange();
  //update the mask
  computeFloatingMask(m_surface);

  //thickness over flotation
 for (DataIterator dit(m_grids); dit.ok(); ++dit)
   {
     Real rhoi = iceDensity();
     Real rhoo = waterDensity(); 
     Real sl = seaLevel();
     
     FORT_THICKNESSOVERFLOTATION(CHF_FRA1(m_thicknessOverFlotation[dit],0),
				 CHF_CONST_FRA1(m_H[dit],0),
				 CHF_CONST_FRA1(m_topography[dit],0),
				 CHF_CONST_REAL(rhoi),
				 CHF_CONST_REAL(rhoo),
				 CHF_CONST_REAL(sl),
				 CHF_BOX(m_thicknessOverFlotation[dit].box()));
   }
  //next, the surface elevation gradient : need to put in the clipping code
  for (DataIterator dit(m_grids); dit.ok(); ++dit)
    {
      
      const BaseFab<int>& thisMask = m_floatingMask[dit];
      const FArrayBox& thisTopg = m_topography[dit];
      FArrayBox& thisSurf = m_surface[dit];
      const FArrayBox& thisThck = m_H[dit];
      // set surface elevation for any open sea or land regions
      FORT_SSETOPENSURFACE(CHF_FRA1(thisSurf,0),
			  CHF_CONST_FIA1(thisMask,0),
			  CHF_CONST_FRA1(thisTopg,0),
			  CHF_CONST_REAL(m_seaLevel),
			  CHF_BOX(thisSurf.box()));

      //cell centered gradient including the first ghost cell
      Box gridBoxPlus = m_grids[dit]; gridBoxPlus.grow(1);

      int xDir = 0;
      if (SpaceDim == 3) xDir = 1;
      int yDir = xDir + 1;
      if (SpaceDim == 1) yDir = 0;

      for (int dir = xDir; dir <= yDir; ++dir)
	{
	  
	  Box faceBox = gridBoxPlus;
	  faceBox.growHi(dir,1);
	  // FORT_SGLGRADS needs a workspace at cell faces
	  FArrayBox faceS(faceBox,1);
	  Real ratio = iceDensity() / waterDensity() ; 
	  FORT_SGLGRADS(CHF_FRA1(m_gradSurface[dit],dir),
	  	       CHF_FRA1(faceS,0),
	  	       CHF_CONST_FRA1(thisThck,0),
	  	       CHF_CONST_FRA1(thisSurf,0),
	  	       CHF_CONST_FRA1(thisTopg,0),
	  	       CHF_CONST_FIA1(thisMask,0),
	  	       CHF_CONST_REAL(ratio),
	  	       CHF_CONST_REAL(m_seaLevel),
	  	       CHF_CONST_REAL(m_dx[dir]),
	  	       CHF_CONST_INT(dir),
	  	       CHF_BOX(gridBoxPlus),
	  	       CHF_BOX(faceBox));

	  m_gradSurface[dit].plus(m_backgroundSlope[dir],dir);

	} //end loop over derivative directions
      
      //face centered gradient
      for (int faceDir = xDir; faceDir <= yDir; ++faceDir)
	{
	  for (int derivDir = xDir; derivDir <= yDir; ++derivDir)
	    {
	      FArrayBox& faceGrad = m_gradSurfaceFace[dit][faceDir];
	      const Box& faceBox = faceGrad.box();
	      Real dx = m_dx[faceDir];
	      FORT_FACEDERIV(CHF_FRA1(faceGrad,derivDir),
			     CHF_CONST_FRA1(m_surface[dit],0),
			     CHF_BOX(faceBox),
			     CHF_CONST_REAL(dx),
			     CHF_INT(derivDir),
			     CHF_INT(faceDir));
	      faceGrad.plus(m_backgroundSlope[derivDir],derivDir);

	    } //end loop over derivative directions
	} //end loop over face directions
    } //end loop over boxes


}


void
LevelSigmaCS::computeFloatingMask(const LevelData<FArrayBox>& a_surface)
{
  CH_TIME("LevelSigmaCS::computeFloatingMask");

  // first, need to create temporary surface height (I don't think 
  // we actually need to store this anywhere, however)
  //LevelData<FArrayBox> surfaceHeight(m_Hgrids, 1, m_H.ghostVect());
  //getSurfaceHeight(surfaceHeight);

  //assume surface has been computed...
  const LevelData<FArrayBox>& surfaceHeight = a_surface;
  // now set up mask of where ice is floating
  DataIterator dit= m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      BaseFab<int>& thisMask = m_floatingMask[dit];
      const FArrayBox& thisZsurf = surfaceHeight[dit];
      const FArrayBox& zBase = m_topography[dit];
      const FArrayBox& thisH = m_H[dit];

      // set to "not floating"
      thisMask.setVal(GROUNDEDMASKVAL);
      int localAnyFloating = 0;
      
      if (!m_grids.physDomain().isPeriodic() || 
	  m_backgroundSlope == RealVect::Zero)
	{
	  Box box = thisZsurf.box();
	  
	  FORT_SETFLOATINGMASK(CHF_FIA1(thisMask,0),
			       CHF_CONST_FRA1(thisZsurf,0),
			       CHF_CONST_FRA1(zBase,0),
			       CHF_CONST_FRA1(thisH,0),
			       CHF_INT(localAnyFloating),
			       CHF_REAL(m_iceDensity),
			       CHF_REAL(m_waterDensity),
			       CHF_REAL(m_seaLevel),
			       CHF_BOX(box));
	}
      if (localAnyFloating) 
	{
	  
	  m_anyFloating[dit] = true;
	}
      else 
	{
	  m_anyFloating[dit] = false;
        }

    } // end loop over grids

}



/// return cell-centered $\Delta_{\tilde{x}}$ and $\Delta_{\tilde{y}}$,
/** stored as components 0 and 1, respectively. */
const LevelData<FArrayBox>& 
LevelSigmaCS::deltaFactors() const
{
  return m_deltaFactors;
}

/// return face-centered $\Delta_{\tilde{x}}$ and $\Delta_{\tilde{y}}$
/** stored as components 0 and 1, respectively. */
const LevelData<FluxBox>& 
LevelSigmaCS::faceDeltaFactors() const
{
  return m_faceDeltaFactors;
}

void WriteSigmaMappedUGHDF5(const string&               a_fileRoot,
                            const DisjointBoxLayout&    a_grids,
                            const LevelData<FArrayBox>& a_data,
                            const LevelSigmaCS&  a_CoordSys,
                            const Box& a_domainBox)
{
  // make up default names and call the version with data names...
  int nComp = a_data.nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
    {
      char labelChSt[80];
      sprintf(labelChSt, "component_%d", n);
      string label(labelChSt);
      compNames[n] = label;
    }

  WriteSigmaMappedUGHDF5(a_fileRoot,                    
                         a_grids,
                         a_data,
                         a_CoordSys,
                         compNames,
                         a_domainBox);
}




void WriteSigmaMappedUGHDF5(const string&               a_fileRoot,
                            const DisjointBoxLayout&    a_grids,
                            const LevelData<FArrayBox>& a_data,
                            const LevelSigmaCS&  a_CoordSys,
                            const Vector<string>& a_compNames,
                            const Box& a_domainBox)
{
  // implementation for LevelSigmaCS
  
  // simplest thing to do here is to just call the AMR version...
  const Vector<DisjointBoxLayout> vectGrids(1, a_grids);
  const Vector<LevelData<FArrayBox>* > vectData(1, const_cast<LevelData<FArrayBox>* >(&a_data) );
  const Vector<const LevelSigmaCS*> vectCoordSys(1, &a_CoordSys);
  Vector<int> refRatio(1,1);
  int numLevels = 1;
  
  // create names for the variables and placeholder values for dt and time,
  Real dt = 1.0;
  Real time = 1.0;
  //int nComp = a_data.nComp();
  
  WriteSigmaMappedAMRHierarchyHDF5(a_fileRoot, vectGrids, vectData,
                                   a_compNames, vectCoordSys, 
                                   a_domainBox, dt, time, refRatio,
                                   numLevels);
  
}


void
WriteSigmaMappedAMRHierarchyHDF5(const string& a_fileRoot,
                                 const Vector<DisjointBoxLayout>& a_vectGrids,
                                 const Vector<LevelData<FArrayBox>* > & a_vectData,
                                 const Vector<string>& a_vectNames,
                                 const Vector<const LevelSigmaCS* >& a_vectCoordSys,
                                 const Box& a_baseDomain,
                                 const Real& a_dt,
                                 const Real& a_time,
                                 const Vector<int>& a_vectRatio,
                                 const int& a_numLevels)
{
  // first, write out data to "regular" hdf5 file:
  char iter_str[80];

  sprintf(iter_str, "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(iter_str);

  RealVect dx = (*a_vectCoordSys[0]).dx();
#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(dataFileName,
                        a_vectGrids,
                        a_vectData,
                        a_vectNames,
                        a_baseDomain,
                        dx[0],
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);
#endif

  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(a_numLevels, NULL);
  //int dimension = (*a_vectCoordSys[0])[dit0].dimension();
  
  // make this a 3-dimensional location of the ice surface...
  int dimension = 3;
  for (int level=0; level<a_numLevels; level++)
    {
      const DisjointBoxLayout& levelGrids = a_vectGrids[level];
      // use same ghosting as cell-centered data used
      IntVect ghostVect = a_vectData[level]->ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        dimension,
                                                        ghostVect);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];
      const LevelSigmaCS& levelCS = *(a_vectCoordSys[level]);

      levelCS.getNodeRealCoordinates(levelNodeData);
          
      if (SpaceDim == 2)
        {
          // now fill in 3rd component with z_surface (for lack of a
          // better idea)
          const LevelData<FArrayBox>& zb = levelCS.getTopography();
          const LevelData<FArrayBox>& H = levelCS.getH();
          
          DataIterator dit = levelGrids.dataIterator();          
          for (dit.begin(); dit.ok(); ++dit)
            {
              NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();

              Box growBox(thisNodeFAB.box());
              growBox.grow(1);
              FArrayBox cellZsurf(growBox, 1);
              
              cellZsurf.copy(zb[dit], growBox, 0, growBox, 0, 1);
              cellZsurf.plus(H[dit], growBox, 0, 0, 1);
          
              // now average cells to nodes
              Real weight = 0.25;
              IntVect offsetX = BASISV(0);
              IntVect offsetY = BASISV(1);
              IntVect offsetXY = IntVect::Unit;

              BoxIterator bit(thisFAB.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  thisFAB(iv,2) = weight*(cellZsurf(iv,0) + cellZsurf(iv-offsetX,0)
                                          +cellZsurf(iv-offsetY,0) + cellZsurf(iv-offsetXY,0));
                }
            } // end loop over grids
        } // end if SpaceDIm == 2
            
    } // end loop over levels

  // create names
  Vector<string> locationNames(dimension);
  for (int i=0; i<dimension; i++)
    {
      if (i ==0) 
        {
          locationNames[0] = "x";
        }
      else if (i == 1)
        {
          locationNames[1] = "y";
        }
      else if (i == 2)
        {
          locationNames[2] = "z";
        }
      else if (i == 3)
        {
          locationNames[3] = "u";
        }
      else if (i == 4)
        {
          locationNames[4] = "v";
        }
      else if (i == 5)
        {
          locationNames[5] = "w";
        }
    }

  sprintf(iter_str, "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(iter_str);
#ifdef CH_USE_HDF5
  // now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(gridInfoFileName,
                        a_vectGrids,
                        vectNodeLoc,
                        locationNames,
                        a_baseDomain,
                        dx[0],
                        a_dt,
                        a_time,
                        a_vectRatio,
                        a_numLevels);
#endif
  // clean up after ourselves here
  for (int level=0; level<a_numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}

/// returns a new LevelSigmaCS which is a coarsened version of this one
/// (useful for multigrid, etc)
LevelSigmaCS* 
LevelSigmaCS::makeCoarser(int a_coarsenFactor) const                          
{
  LevelSigmaCS* coarseCS = new LevelSigmaCS;

  coarseCS->define(*this, a_coarsenFactor);

  return coarseCS;
}



void LevelSigmaCS::interpFromCoarse(const LevelSigmaCS& a_crseCoords, 
				    const int a_refinementRatio, 
				    const bool a_interpolateTopography, 
				    const bool a_interpolateThickness, 
				    const bool a_preserveMask,
				    const bool a_interpolateTopographyGhost, 
				    const bool a_interpolateThicknessGhost, 
				    const bool a_preserveMaskGhost, 
				    int a_thicknessInterpolationMethod)
{
  CH_TIME("LevelSigmaCS::interpFromCoarse");
  int ncomp = m_H.nComp();
  const IntVect& ghost = m_H.ghostVect();
  int nghost = ghost[0];
  
  FineInterp interpolator(m_grids, 
			  ncomp, 
			  a_refinementRatio , 
			  m_grids.physDomain());
  
  const LevelData<FArrayBox>& crseTopg = a_crseCoords.m_topography;
  LevelData<FArrayBox>& fineTopg = m_topography;
  const LevelData<FArrayBox>& crseThck = a_crseCoords.m_H;
  LevelData<FArrayBox>& fineThck = m_H;
  const DisjointBoxLayout& crseGrids = a_crseCoords.m_grids;
  LevelData<FArrayBox> fineSurf(m_grids, ncomp, ghost);
  LevelData<FArrayBox> fineBase(m_grids, ncomp, ghost);
  LevelData<FArrayBox> crseSurf(crseGrids, ncomp, ghost);
  LevelData<FArrayBox> crseBase(crseGrids, ncomp, ghost);
  const LevelData<BaseFab<int> >& crseMask = a_crseCoords.m_floatingMask;
  LevelData<BaseFab<int> >& fineMask = m_floatingMask;

  //need to compute surface elevation on the coarse level
  a_crseCoords.getSurfaceHeight(crseSurf);
  for (DataIterator dit(crseGrids); dit.ok(); ++dit)
    {
      crseBase[dit].copy(crseSurf[dit]);
      crseBase[dit].minus(crseThck[dit]);
    }

  //copy the mask from the coarse level
  IntFineInterp ifi(m_grids,fineMask.nComp(),
		    a_refinementRatio,
		    fineMask.ghostVect(),
		    m_grids.physDomain());
  ifi.pwcInterpToFine(fineMask, crseMask);


  //interpolate the fields, not including ghost regions
  if (a_interpolateTopography)
    {
      interpolator.interpToFine(fineTopg,crseTopg);
    }
   if (a_interpolateThickness)
     {
       CH_assert(a_thicknessInterpolationMethod < MAX_THICKNESS_INTERPOLATION_METHOD);
       
       interpolator.interpToFine(fineSurf,crseSurf);
       interpolator.interpToFine(fineBase,crseBase);

       if (a_thicknessInterpolationMethod == SMOOTH_SURFACE_THICKNESS_INTERPOLATION_METHOD)
	 {
	   // determine the thickness from the interpolated surface and the fine topography
	   // in grounded ice regions : the idea is to avoid laying a smooth thickness on top
	   // of a newly noisy bedrock and getting a noisy surface :
	   for (DataIterator dit(m_grids); dit.ok(); ++dit)
	     {
	       
	       const BaseFab<int>& mask = fineMask[dit];
	       const FArrayBox& topg = fineTopg[dit];
	       FArrayBox& lsrf = fineBase[dit];
	       const Box& b = m_grids[dit];
	       for (BoxIterator bit(b); bit.ok(); ++bit)
		 {
		   const IntVect& iv = bit();
		   if (mask(iv) == GROUNDEDMASKVAL)
		     {
		       lsrf(iv) = topg(iv);
		     }
		 }
	     };	  
	 }
     }
   else 
     {
       //assumes that thickness has been set
       getSurfaceHeight(fineSurf);
       for (DataIterator dit(m_grids); dit.ok(); ++dit)
	 {
	   fineBase[dit].copy(fineSurf[dit]);
	   fineBase[dit].minus(fineThck[dit]);
	 }
     }

  PiecewiseLinearFillPatch ghostFiller(m_grids,
				       crseGrids,
				       ncomp,
				       crseGrids.physDomain(),
				       a_refinementRatio,
				       nghost);
  // now include the ghost regions
  if (a_interpolateTopographyGhost)
    {
      ghostFiller.fillInterp(fineTopg,crseTopg,crseTopg, 0.0,0,0,ncomp);
      exchangeTopography();
    }
  if (a_interpolateThicknessGhost)
    {
      ghostFiller.fillInterp(fineSurf,crseSurf,crseSurf, 0.0,0,0,ncomp);
      ghostFiller.fillInterp(fineBase,crseBase,crseBase, 0.0,0,0,ncomp);
      fineSurf.exchange();
      fineBase.exchange();

    }
  



   bool maskEdgeInterpolate = true;
  // bool maskEdgeInterpolate = false;
  if (maskEdgeInterpolate)
    {
      //just copying the mask from the coarse to fine level leaves us with
      //coarse grid stair-stepping. However, in practical applications 
      //we only have say, 5km data but nonetheless would like smooth divisions 
      //between land ice, ice shelf, open land and open sea. Here, 
      //we attempt to achieve that, while ensuring 
      //that a margin parallel to x or y would not move 
      //(which would happen if we simply 
      //interpolate thickness and topography)
      for (DataIterator dit(m_grids); dit.ok(); ++dit)
	{
	  Box b = m_grids[dit];
	  
	  BaseFab<int>& f = fineMask[dit];
	  BaseFab<int> g(b,1);
	  g.copy(f);
	  int test[4] = {GROUNDEDMASKVAL,FLOATINGMASKVAL,OPENSEAMASKVAL,OPENLANDMASKVAL};
          int ydir = min(1,SpaceDim-1);

	  const IntVect ex = BASISV(0);
	  const IntVect ey = BASISV(ydir);
	  
	  for (BoxIterator bit(b);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      for (int i =0; i < 4; i++)
		{
		  const int a = test[i];
		  if (5 == ( ((f(iv+ex)==a)?1:0)
			    + ((f(iv-ex)==a)?1:0)
			    + ((f(iv+ex+ey)==a)?1:0)
			    + ((f(iv+ex-ey)==a)?1:0)
			    + ((f(iv-ex+ey)==a)?1:0)
			    + ((f(iv-ex-ey)==a)?1:0)
			    + ((f(iv+ey)==a)?1:0)
			    + ((f(iv-ey)==a)?1:0)))
		    {
		      g(iv) = a;
		    }
		}
	      
	    }
	  f.copy(g);
	}
      fineMask.exchange();
    }
	
  //adjust surfaces and topography such that same mask 
  //will be computed. compute the thickness accordingly
  for (DataIterator dit(m_grids); dit.ok(); ++dit)
    {
     
      Real id = m_iceDensity; Real wd = m_waterDensity; Real sl = m_seaLevel;
      if (a_preserveMask)
	{
	  Box b = m_grids[dit];
	  if (a_preserveMaskGhost)
	    b.grow(nghost);

	  
	 
	  FArrayBox oldThck(b, 1);
	  FArrayBox oldTopg(b, 1);
	  FArrayBox oldSurf(b, 1);
	  FArrayBox oldBase(b, 1);
	  oldThck.copy(fineSurf[dit]);
	  oldThck.minus(fineBase[dit]);
	  oldBase.copy(fineBase[dit]);
	  oldSurf.copy(fineSurf[dit]);
	  oldTopg.copy(fineTopg[dit]);
	  // don't make modifications outside the domain
	  b &= m_grids.physDomain().domainBox(); 

	 
	  BaseFab<int> oldMask(b, 1);
	  int anyFlt;
	  FORT_SETFLOATINGMASK(CHF_FIA1(oldMask,0),
			       CHF_CONST_FRA1(oldSurf,0),
			       CHF_CONST_FRA1(oldBase,0),
			       CHF_CONST_FRA1(oldThck,0),
			       CHF_INT(anyFlt),
			       CHF_REAL(m_iceDensity),
			       CHF_REAL(m_waterDensity),
			       CHF_REAL(m_seaLevel),
			       CHF_BOX(b));


	  FORT_PRESERVEMASK(CHF_FRA1(fineBase[dit],0),
			    CHF_FRA1(fineSurf[dit],0),
			    CHF_FRA1(fineTopg[dit],0),
			    CHF_FRA1(fineThck[dit],0),
			    CHF_CONST_FIA1(fineMask[dit],0),
			    CHF_REAL(id),
			    CHF_REAL(wd),
			    CHF_REAL(sl),
			    CHF_BOX(b));

	  const Box& lapBox =  m_grids[dit];
	  Real dx = m_dx[0]; CH_assert(SpaceDim == 1 || m_dx[0] == m_dx[1]);
	  FArrayBox lapOldSurf(lapBox,1);lapOldSurf.setVal(0.0);
	  FArrayBox lapNewSurf(lapBox,1);lapNewSurf.setVal(0.0);
	  Real alpha = 0.0; Real beta = 1.0;

	  FORT_OPERATORLAP(CHF_FRA(lapOldSurf),
                           CHF_FRA(oldSurf),
                           CHF_BOX(lapBox),
                           CHF_CONST_REAL(dx),
                           CHF_CONST_REAL(alpha),
                           CHF_CONST_REAL(beta));

	  FORT_OPERATORLAP(CHF_FRA(lapNewSurf),
                           CHF_FRA(fineSurf[dit]),
                           CHF_BOX(lapBox),
                           CHF_CONST_REAL(dx),
                           CHF_CONST_REAL(alpha),
                           CHF_CONST_REAL(beta));

	  for (BoxIterator bit(lapBox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      //avoid creating areas of open land that result in new surface minima
	      if (fineMask[dit](iv) == OPENLANDMASKVAL && oldMask(iv) == GROUNDEDMASKVAL)
		{
		  if (abs( lapOldSurf(iv)) < 0.5 * abs(lapNewSurf(iv)))
		    {
		      fineSurf[dit](iv) = oldSurf(iv);
		      fineTopg[dit](iv) = oldTopg(iv);
		      fineThck[dit](iv) = oldThck(iv);
		    } 
		}
	    }
	 
	}
      else if (!a_preserveMask && a_preserveMaskGhost)
	{
	  Box b = m_grids[dit]; 
	  b.grow(nghost);
	  // don't make modifications outside the domain
	  b &= m_grids.physDomain().domainBox(); 

	  FArrayBox tmpThck(b, 1);
	  FArrayBox tmpTopg(b, 1);
       
	  tmpThck.copy(fineSurf[dit]);
	  tmpThck.minus(fineBase[dit]);
	  tmpTopg.copy(fineTopg[dit]);
	 
	  FORT_PRESERVEMASK(CHF_FRA1(fineBase[dit],0),
			    CHF_FRA1(fineSurf[dit],0),
			    CHF_FRA1(fineTopg[dit],0),
			    CHF_FRA1(fineThck[dit],0),
			    CHF_CONST_FIA1(fineMask[dit],0),
			    CHF_REAL(id),
			    CHF_REAL(wd),
			    CHF_REAL(sl),
			    CHF_BOX(b));
	  fineThck[dit].copy(tmpThck,m_grids[dit]);
	  fineTopg[dit].copy(tmpTopg,m_grids[dit]);
	}
     
     
      //     CH_assert(fineThck[dit].min() >= 0.0);
      //   CH_assert(fineThck[dit].max() <= 1.0e+6);
    }

  // redo this, beacause we want the ghost thickness(rather than the surface or base).  
  // to be consistent with the coarse level 
  if (a_interpolateThicknessGhost)
    {
      ghostFiller.fillInterp(fineThck,crseThck,crseThck, 0.0,0,0,ncomp);
      fineThck.exchange();
    }
  //const ProblemDomain& domain = m_grids.physDomain();
  // for (DataIterator dit(m_grids); dit.ok(); ++dit)
  //   {
  //     // don't worry about ghost cells outside the domain
  //     Box validBox(fineThck[dit].box());
  //     validBox &= domain;
      
  //     CH_assert(fineThck[dit].min(validBox,0) >= 0.0);
  //     CH_assert(fineThck[dit].max(validBox,0) <= 1.0e+6);
  //   }
}
void LevelSigmaCS::exchangeTopography()
{
  m_topography.exchange();
}


void LevelSigmaCS::unitShiftExchange(LevelData<FArrayBox>& a_level)
{
  //exchange a_level, then adjust for unit shift
  a_level.exchange();
  const ProblemDomain& domain = m_Hgrids.physDomain();
  //Now correct the ghost regions
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (domain.isPeriodic(dir))
	{
	  Box lo = adjCellLo(domain.domainBox(), dir, a_level.ghostVect()[dir]);
	  Box hi = adjCellHi(domain.domainBox(), dir, a_level.ghostVect()[dir]);
	  
	  //we want to adjust corner cells too, so
	  IntVect v = a_level.ghostVect();
	  v[dir] = 0;lo.grow(v);hi.grow(v);

	  DataIterator dit = a_level.dataIterator();
	  for (dit.begin(); dit.ok(); ++dit)
	    {
	      FArrayBox& b = a_level[dit];
	      Box hiBox = b.box();
	      hiBox &= hi;
	      if (!hiBox.isEmpty())
		b.plus(m_unitShift[dir],hiBox);
	      Box loBox = b.box();
	      loBox &= lo;
	      if (!loBox.isEmpty())
		b.plus(-m_unitShift[dir],loBox);
	    }
	}

    }
  int breakpoint = 0; breakpoint++;
}




#include "NamespaceFooter.H"
