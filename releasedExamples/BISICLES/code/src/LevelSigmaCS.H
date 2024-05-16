#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _LEVELSIGMACS_H_
#define _LEVELSIGMACS_H_

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "NodeFArrayBox.H"
#include "FluxBox.H"
#include "RealVect.H"
#include "Interval.H"
#include "AMRIO.H"


#include "NamespaceHeader.H"

#define BISICLES_LAYERED 0
#define BISICLES_FULLZ 1
#if CH_SPACEDIM == 2
// for now, BISICLES_LAYERED is the same as 2D, but
// if we one day look at flow line models there could a 1D variant
#define BISICLES_Z 0
#elif CH_SPACEDIM == 3
// for now, BISICLES_FULLZ is the same as 3D, but
// if we one day look at flow line models there could a 2D variant
#define BISICLES_Z 1
#elif CH_SPACEDIM == 1
// for now, same as 2D
// for now, BISICLES_LAYERED is the same as 2D, but
// if we one day look at flow line models there could a 1D variant
#define BISICLES_Z 0
#endif

/// Basic Sigma fourth-order coordinate system on an AMR level
/**
   The SigmaCS class implements Sigma-coordinates for ice sheets 
   on a DisjointBoxLayout
*/
class LevelSigmaCS
{
public:
  /// default constructor
  LevelSigmaCS();

  /// defining constructor -- calls matching define function
  LevelSigmaCS(const DisjointBoxLayout& a_grids,
               const RealVect& a_dx, 
               const IntVect& a_ghostVect = IntVect::Unit);

  LevelSigmaCS( const DisjointBoxLayout& a_grids,
		const RealVect& a_dx,
		const LevelSigmaCS& a_fineCS, 
		int a_nRef
		);
  /**
     Destructor.
   */
  virtual ~LevelSigmaCS();

  void define(const DisjointBoxLayout& a_grids,
              const RealVect& a_dx, 
              const IntVect& a_ghostVect = IntVect::Unit);

  /// define as a coarsening of fineCS by nRef
  void define(const LevelSigmaCS& a_fineCS, 
              int a_nRef);

  /// given coordinate in mapped space, return its location in real
  /// space
  /** note that a_Xi is (sigma,x,y) while realCoord returns (x,y,z)
      (in 3d -- in 2d a_Xi is (x,y) )
   */
  RealVect realCoord(const RealVect& a_Xi, const DataIndex& a_index) const;
  
  /// given coordinate in real space, return location in mapped space  
  /** note that a_x is (x,y,z) while mappedCoord returns (sigma,x,y)
      in 3d. In 2d, mappedCoord returns (x,y)
   */  
  RealVect mappedCoord(const RealVect& a_x, const DataIndex& a_index) const;

  /// return Cartesian XYZ locations of nodes
  /** nodeCoords should have 3 coordinates
   */
  void getNodeRealCoordinates(LevelData<NodeFArrayBox>& a_nodeCoords) const;


  // return cell spacing dx
  RealVect dx() const {return m_dx;}

  /// returns modifiable cell-centered H (ice sheet thickness) for this
  LevelData<FArrayBox>& getH(); 

  /// returns const reference to cell-centered H (ice sheet thickness) 
  const LevelData<FArrayBox>& getH() const;
  
  /// returns modifiable face-centered H
  LevelData<FluxBox>& getFaceH(); 

  /// returns const reference to face-centered H
  const LevelData<FluxBox>& getFaceH() const;

  /// returns a const-reference to the cell-centered topography
  const LevelData<FArrayBox>& getTopography() const;

  /// returns a modifiable reference to the cell-centered topography
  LevelData<FArrayBox>& getTopography();

  ///sets the base height. 
  /** In practice, this will probably done by a derived class.*/
  void setTopography(const LevelData<FArrayBox>& a_topography);
  /// set the unit topography difference 
  /**  for periodic problems :  a vector with components
  (topography(x + Lx , y) - topography(x , y), 
  (topography(x , y + Ly) - topography(x , y) */
  //void setUnitShift(const RealVect& a_unitShift)
  //{
  //  m_unitShift = a_unitShift;
  //}


  /// set the background slope 
  /**  for periodic problems entirely on land we need not have periodic 
       topography. Instead we can define a vector slope such that
       
       true topography(x,y) =  stored topography(x,y) + slope_x*x + slope_y*y 
       
       can be specified. In that case, 

       surface(x,y) = stored surface(x,y) + slope_x*x + slope_y*y
       
       and

       grad(surface) = grad(stored surface(x,y)) + slope
  */

  void setBackgroundSlope(const RealVect& a_backgroundSlope)
  {
    m_backgroundSlope = a_backgroundSlope;
  }

  const RealVect& getBackgroundSlope() const
  {
    return m_backgroundSlope;
  }


  //sets the surface elevation (at cell centres)
  void setSurfaceHeight(const LevelData<FArrayBox>& a_surface);
  //sets the surface gradient (at cell centres)
  void setGradSurface(const LevelData<FArrayBox>& a_gradSurface);
  //sets the surface gradient (at cell faces)
  void setGradSurfaceFace(const LevelData<FluxBox>& a_gradSurface);
  
  
  //returns the surface elevation (at cell centres)
  const LevelData<FArrayBox>& getSurfaceHeight() const ;
  //returns the surface gradient (at cell centres)
  const LevelData<FArrayBox>& getGradSurface() const;
  //returns the surface gradient (at cell faces)
  const LevelData<FluxBox>& getGradSurfaceFace() const;
  //returns the thickness over flotation (at cell centres)
  const LevelData<FArrayBox>& getThicknessOverFlotation() const;

  /// computes ice surface height
  /** 
      if ice is grounded, then this equals z_base + H
      if ice is floating, then this equals seaLevel + H*(1 - rho_ice/rho_water)
  */
  void getSurfaceHeight(LevelData<FArrayBox>& a_zSurface) const;

  /// returns floating mask specifying whether ice is grounded or floating
  /** following mask values are defined in IceConstants.H
      floatingmaskval -- ice shelf (nonzero thickness floating on water) 
      groundedmaskval -- bottom of ice sheet is grounded.
      openseamaskval  -- zero ice thickness, base height below sea level
      openlandmaskvel -- zero ice thickness, base height above sea level
  */
  const LevelData<BaseFab<int> >& getFloatingMask() const {return m_floatingMask;}

  /// recomputes quantities which are dependent on cell-centered H. 
  /** Should be called after cell-centered H is modified to ensure 
      that things don't get out of sync. Also computes face-centered H 
      by averaging to faces */      
  void recomputeGeometry( const LevelSigmaCS* a_crseCoords, 
			  const int a_refRatio) ;
 
  /// recomputes quantities which are dependent on face-centered H. 
  /** Should be called after cell-centered H is modified to ensure 
      that things don't get out of sync. Also computes cell-centered 
      H as average of face-centered H  */      
  void recomputeGeometryFace( const LevelSigmaCS* a_crseCoords, 
			      const int a_refRatio) ;

  /// return cell-centered \f$\Delta_{\tilde{x}}\f$ and \f$\Delta_{\tilde{y}}\f$,
  /** stored as components 0 and 1, respectively. */
  const LevelData<FArrayBox>& deltaFactors() const;

  /// return face-centered \f$\Delta_{\tilde{x}}\f$ and \f$\Delta_{\tilde{y}}\f$
  /** stored as components 0 and 1, respectively. */
  const LevelData<FluxBox>& faceDeltaFactors() const;


  /// returns sea level
  Real seaLevel() const {return m_seaLevel;}

  /// set sea level
  void setSeaLevel(const Real& a_seaLevel) {m_seaLevel = a_seaLevel;}

  /// ice density -- this is probably an interim approach
  Real iceDensity() const {return m_iceDensity;}
  /// seawater density -- probably an interim approach
  Real waterDensity() const {return m_waterDensity;}
  /// acceleration due to gravity --  probably an interim approach
  Real gravity() const {return m_gravity;}

  /// interim or not, we also need to set these
  void setIceDensity(const Real& a_iceDensity) 
  {m_iceDensity = a_iceDensity;}
  void setWaterDensity(const Real& a_waterDensity) 
  {m_waterDensity = a_waterDensity;}
  void setGravity(const Real& a_gravity) 
  {m_gravity = a_gravity;}
  

  ///returns true if this object has been defined.
  bool isDefined() const {return m_isDefined;}
  
  /// returns a new LevelSigmaCS which is a coarsened version of this one
  /// (useful for multigrid, etc)
  LevelSigmaCS* makeCoarser(int a_coarsenFactor) const;
                   
  /// access function
  const DisjointBoxLayout& grids() const {return m_grids;}

  const LayoutData<bool>&  anyFloating() const { return m_anyFloating;}

  /// choice of thickness interpolation method in interpFromCoarse
  /**
     The original (default) method, STD_THICKNESS_INTERPOLATION_METHOD,
     interpolates the thickness directly, which can lead to a smooth
     thickness on a rough bedrock, and henec a rough surface

     An alternative method, SMOOTH_SURFACE_THICKNESS_INTERPOLATION_METHOD, 
     interpolates the upper and lower surfaces, but replaces the lower surface
     on grounded (on the coarse mesh) cells with the bedrock data. This
     will be conservative iff the relationship between fine and coarse
     bedrock data is conservative

   */
  enum ThicknessInterpolationMethod{
    STD_THICKNESS_INTERPOLATION_METHOD = 0,
    SMOOTH_SURFACE_THICKNESS_INTERPOLATION_METHOD = 1,
    MAX_THICKNESS_INTERPOLATION_METHOD = 2};

  /// fill the data in this by interpolation from a coarser level
  /*** 
       Interpolating thickness and topography independently tends to
       move the grounding line, and results in steep surfaces there
       which are not seen in interpolation of the surface itself. 
       interpFromCoarse interpolates topography, the upper and lower surfaces
       then adjusts base (and so thickness) to maintain the grounding line.
       
       optionally, topography and thickness can be set elsewhere,
       and the mask can be recomputed, in either the valid region, or in
       both valid region and ghost region. 
 */
  void interpFromCoarse(const LevelSigmaCS& a_crseCoords,
			const int a_refinementRatio, 
			const bool a_interpolateTopography = true, 
			const bool a_interpolateThickness = true, 
			const bool a_preserveMask = true,
			const bool a_interpolateTopographyGhost = true, 
			const bool a_interpolateThicknessGhost = true, 
			const bool a_preserveMaskGhost = true, 
			int a_thicknessInterpolationMethod = STD_THICKNESS_INTERPOLATION_METHOD);





  // exchange topography accounting for the unit shift.
  void exchangeTopography();
  // exchange a_level accounting for the unit shift.
  void unitShiftExchange(LevelData<FArrayBox>& a_level);

  /// accessor for ghost vector
  IntVect ghostVect() const {return m_ghostVect;}

protected:

  ///
  DisjointBoxLayout m_grids;

  /// "flattened" grids (in 3d) for quantities like H
  DisjointBoxLayout m_Hgrids;

  // cell spacing
  RealVect m_dx;

  /// ghost vector 
  IntVect m_ghostVect;

  // for periodic problems :  a vector with components
  // (topography(x + Lx , y) - topography(x , y), 
  // (topography(x , y + Ly) - topography(x , y)
  RealVect m_unitShift;

  RealVect m_backgroundSlope;

  /// cell-centered topography
  LevelData<FArrayBox> m_topography;

  /// cell-centered ice thickness
  LevelData<FArrayBox> m_H;

  /// face-centered ice thickness
  LevelData<FluxBox> m_faceH;

  /// cell-centered  surface elevation
  LevelData<FArrayBox> m_surface;

  /// cell-centered gradient of surface elevation
  LevelData<FArrayBox> m_gradSurface;

  /// face-centered gradient of surface elevation
  LevelData<FluxBox> m_gradSurfaceFace;

  /// cell-centered thickness over flotation
  LevelData<FArrayBox> m_thicknessOverFlotation;

  /// floating mask 
  LevelData<BaseFab<int> > m_floatingMask;

  /// is anything floating?
  LayoutData<bool> m_anyFloating; 
  ///
  LevelData<FArrayBox> m_deltaFactors;

  ///
  LevelData<FluxBox> m_faceDeltaFactors;

  /// has this object been defined?
  bool m_isDefined;


  // sea level
  Real m_seaLevel;

  /// ice density
  Real m_iceDensity;

  /// seawater density -- probably an interim approach
  Real  m_waterDensity;

  /// acceleration due to gravity -- probably an interim approach
  Real  m_gravity;

  // set a reasonable set of default values
  void setDefaultValues();
  
  /// (re)compute delta factors
  void computeDeltaFactors();

  /// (re)compute surface elevation and gradients
  void computeSurface(const LevelSigmaCS* a_crseCoords, 
		      const int a_refRatio);

  /// sets floating mask based on surface height
  void computeFloatingMask(const LevelData<FArrayBox>& a_surface);


#if BISICLES_Z == BISICLES_LAYERED

public:
  // some poor man's multidim stuff that is only useful when we
  // are treating the problem as a stack of layers with constant
  // sigma at the boundaries which never get refined

  // vector of sigma values at the layer interfaces  
  const Vector<Real>& getFaceSigma() const
  { return m_faceSigma; }

  //vector of delta sigma's (at the layer midpoints)
  const Vector<Real>& getDSigma() const
  { return m_dSigma;}

  // vector of sigma values at the layer midpoints
  const Vector<Real>& getSigma() const
  { return m_sigma; }

  void setFaceSigma(const Vector<Real>& a_faceSigma)
  {
    m_faceSigma = a_faceSigma;
    m_sigma.resize(m_faceSigma.size()-1);
    m_dSigma.resize(m_faceSigma.size()-1);
    for (unsigned int i = 0; i <  m_sigma.size(); ++i)
      {
	//layer midpoint
	m_sigma[i] = 0.5 * (a_faceSigma[i] + a_faceSigma[i+1]);
	//layer delta sigma
	m_dSigma[i] = (a_faceSigma[i+1] - a_faceSigma[i]);
      }
  }

protected:

  Vector<Real> m_faceSigma;
  Vector<Real> m_sigma;
  Vector<Real> m_dSigma;

#endif


};


// specializations of IO routines for LayoutData<SigmaCS>

void WriteSigmaMappedUGHDF5(const string&               a_fileRoot,
                            const DisjointBoxLayout&    a_grids,
                            const LevelData<FArrayBox>& a_data,
                            const LevelSigmaCS&  a_CoordSys,
                            const Box& a_domainBox);

void WriteSigmaMappedUGHDF5(const string&               a_fileRoot,
                            const DisjointBoxLayout&    a_grids,
                            const LevelData<FArrayBox>& a_data,
                            const LevelSigmaCS&  a_CoordSys,
                            const Vector<string>& a_vectNames,
                            const Box& a_domainBox);

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
                                 const int& a_numLevels);

#include "NamespaceFooter.H"

#endif
