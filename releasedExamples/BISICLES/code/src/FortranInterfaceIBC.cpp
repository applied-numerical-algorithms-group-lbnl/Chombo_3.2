#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LoadBalance.H"
#include "FortranInterfaceIBC.H"
#include "ParmParse.H"
#include "FillFromReference.H"
#include "ExtrapGhostCells.H"
#include "ExtrapBCF_F.H"
#include "ReflectGhostCells.H"
#include "FIBCF_F.H"
#include "NamespaceHeader.H"

void zeroBCValue_FIBC(Real* pos,
		     int* dir,
		     Side::LoHiSide* side,
		     Real* a_values)
{
  a_values[0]=0.0;
  a_values[1]=0.0;
}



void iceNeumannBC_FIBC(FArrayBox& a_state,
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
		  Interval NeumInterval(0,SpaceDim-1);
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Lo,
			 NeumInterval);
                }

              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
		  Interval NeumInterval(0,SpaceDim-1);
                  NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         zeroBCValue_FIBC,
                         dir,
                         Side::Hi, 
			 NeumInterval);
                }

            } // end if is not periodic in ith direction
        }
    }
}


void iceDirichletBC_FIBC(FArrayBox& a_state,
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
                         zeroBCValue_FIBC,
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
                         zeroBCValue_FIBC,
                         dir,
                         Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}

//set all components of u to zero in ghost regions
void iceDivideBC_FIBC(FArrayBox& a_state,
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
              
              if(!a_domain.domainBox().contains(ghostBoxLo))
                {
		  ghostBoxLo &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxLo, 0, a_state.nComp());
		  
                }
	      Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              if(!a_domain.domainBox().contains(ghostBoxHi))
                {
		  ghostBoxHi &= a_state.box();
		  a_state.setVal(0.0,  ghostBoxHi, 0, a_state.nComp());
                }

            } // end if is not periodic in ith direction
        }
    }

}


//ice stream parallel to x : homogeneous Neumann conditions
//along the y faces, and homogeneous Dirichlet conditions
//along the x faces. 
void iceStreamXBC_FIBC(FArrayBox& a_state,
		      const Box& a_valid,
		      const ProblemDomain& a_domain,
		      Real a_dx,
		      bool a_homogeneous)
{

  if(!a_domain.domainBox().contains(a_state.box()))
    {
      Box valid = a_valid;
      
      int xDir = 0;
      if (!a_domain.isPeriodic(xDir))
	{
	  Box ghostBoxLo = adjCellBox(valid, xDir, Side::Lo, 1);
	  
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      ghostBoxLo &= a_state.box();
	      a_state.setVal(0.0,  ghostBoxLo, 0, a_state.nComp());
	      
	    }
	  Box ghostBoxHi = adjCellBox(valid, xDir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      ghostBoxHi &= a_state.box();
	      a_state.setVal(0.0,  ghostBoxHi, 0, a_state.nComp());
	    }
	}

      if (SpaceDim > 1)
	{
	  int yDir = 1;
	  Box ghostBoxLo = adjCellBox(valid, yDir, Side::Lo, 1);
	  Box ghostBoxHi = adjCellBox(valid, yDir, Side::Hi, 1);
	  if(!a_domain.domainBox().contains(ghostBoxLo))
	    {
	      //Real bcVal = 0.0;
	      Interval NeumInterval(0,SpaceDim-1);
	      NeumBC(a_state,
		     valid,
		     a_dx,
		     a_homogeneous,
		     zeroBCValue_FIBC,
		     yDir,
		     Side::Lo,
		     NeumInterval);
	    }
	  
	  if(!a_domain.domainBox().contains(ghostBoxHi))
	    {
	      Interval NeumInterval(0,SpaceDim-1);
	      NeumBC(a_state,
		     valid,
		     a_dx,
		     a_homogeneous,
		     zeroBCValue_FIBC,
		     yDir,
		     Side::Hi, 
		     NeumInterval);
	    }
	  
	}
    }
}




// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
FortranInterfaceIBC::FortranInterfaceIBC() : m_boundaryThickness(0.0)
{
  m_isBCsetUp = false;
  m_isDefined = false;
  m_gridsSet = false;
  // set default to be true
  m_verbose = true;
  m_thicknessGhost = IntVect::Zero;
  m_topographyGhost = IntVect::Zero;
  // default is an empty vector here 
  // (don't set thickness to zero anywhere)
  m_thicknessClearRegions = Vector<Box>();
  
}

FortranInterfaceIBC::~FortranInterfaceIBC()
{
  if (m_verbose)
    {
      pout() << "in FortranInterfaceIBC destructor" << endl;
    }
}


/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
FortranInterfaceIBC::define(const ProblemDomain& a_domain,
			    const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
  m_isDefined = true;
}

// static function which defines cell centred a_fab, directly with cell centered data
// a_data_ptr if !a_nodal, or by averaging nodal data a_data_ptr.
void FortranInterfaceIBC::setFAB(Real* a_data_ptr,
				 const int* a_dimInfo,
                                 const int* a_boxlo, const int* a_boxhi, 
				 const Real* a_dew, const Real* a_dns,
				 const IntVect& a_offset,
				 const IntVect& a_nGhost,
				 FArrayBox & a_fab,
                                 FArrayBox & a_ccFab,
				 const bool a_nodal,
                                 const bool a_verbose)
{
  if (a_verbose)
    {
      pout() << "in FortranInterfaceIBC::setFAB" << endl;
      pout() << " -- doing box calculus" << endl;
    }

  IntVect loVect(D_DECL(a_boxlo[0],a_boxlo[1], a_boxlo[2]-1));
  IntVect hiVect(D_DECL(a_boxhi[0],a_boxhi[1], a_boxhi[2]-1));
  if (a_verbose)
    {
      pout () << "loVect = " << loVect << endl;
      pout () << "hiVect = " << hiVect << endl;
    }

  Box fabBox;
  if (loVect <= hiVect)
    {
      fabBox.define(loVect, hiVect);
      if (a_verbose)
        {
          pout() << "done defining box" << endl;
        }
      fabBox.shift(-a_nGhost);
      fabBox.shift(a_offset);
      if (a_verbose)
        {
          pout() << "... done! " << endl;
        }
      //    }
      
      if (a_nodal)
        {
          
          a_ccFab.define(fabBox, 1);
          Box nodeBox(fabBox);
          nodeBox.surroundingNodes();
          //nodeBox.shift(-a_nGhost);
          FArrayBox nodeFAB;
          a_fab.define(nodeBox,1,a_data_ptr);
          FORT_NODETOCELL(CHF_CONST_FRA1(a_fab,0),
                          CHF_FRA1(a_ccFab,0),
                          CHF_BOX(fabBox));
        }
      else 
        {
          // both of these are aliased to point to the input data
          a_fab.define(fabBox, 1, a_data_ptr);
          a_ccFab.define(fabBox, 1, a_data_ptr);
        }
    }
}


/// version of setFAB for (horizontal) velocity data 
void FortranInterfaceIBC::setVelFAB(Real* a_uVelPtr,
                                    Real* a_vVelPtr,
                                    const int* a_dimInfo,
                                    const int* a_boxlo, const int* a_boxhi, 
                                    const Real* a_dew, const Real* a_dns,
                                    const IntVect& a_offset,
                                    const IntVect& a_nGhost,
                                    FArrayBox & a_uFab,
                                    FArrayBox & a_vFab,
                                    FArrayBox & a_uccFab,
                                    FArrayBox & a_vccFab,
                                    const bool a_nodal,
                                    const bool a_fillData)
{
  if (m_verbose)
    {
      pout() << "in FortranInterfaceIB:setVelFAB" << endl;
      pout() << " -- doing box calculus" << endl;
    }

  IntVect loVect(D_DECL(a_boxlo[0],a_boxlo[1], a_boxlo[2]-1));
  IntVect hiVect(D_DECL(a_boxhi[0],a_boxhi[1], a_boxhi[2]-1));

  // vertical layering for velocity makes this a bit odd --
  // need to rotate everything
  int nLayers = a_dimInfo[1];

  IntVect CISMloVect(D_DECL(0,a_boxlo[0], a_boxlo[1]));
  IntVect CISMhiVect(D_DECL(nLayers-1,a_boxhi[0], a_boxhi[1]));
  int velNcomp = a_boxhi[1]-a_boxlo[1]+1;
  
  if (m_verbose)
    {
      pout () << "loVect = " << loVect << endl;
      pout () << "hiVect = " << hiVect << endl;

      pout() << "CISM Vel loVect = " << CISMloVect << endl;
      pout() << "CISM Vel hiVect = " << CISMhiVect << endl;
    }

  Box fabBox, CISMvelBox;
  if (loVect <= hiVect)
    {
      fabBox.define(loVect, hiVect);
      CISMvelBox.define(CISMloVect,CISMhiVect);
      if (m_verbose)
        {
          pout() << "done defining box" << endl;
        }
      fabBox.shift(-a_nGhost);
      fabBox.shift(a_offset);
      // first direction is "layers", which makes this a bit stranger.
      CISMvelBox.shift(1,-a_nGhost[0]);
      CISMvelBox.shift(1,a_offset[0]);
      if (m_verbose)
        {
          pout() << "... done! " << endl;
        }
    }

  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(fabBox);
      gridBox.grow(-a_nGhost);
      setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }
  if (a_nodal)
    {
      // this is going to be a bit screwy -- this FAB will have the 
      // right amount of storage
      
      a_uccFab.define(CISMvelBox, velNcomp);
      a_vccFab.define(CISMvelBox, velNcomp);
      Box nodeBox(CISMvelBox);
      nodeBox.surroundingNodes();
      // back to cell-centered in "vertical" direction
      nodeBox.enclosedCells(0);

      // 
      int nodeVelNcomp = velNcomp +1;
      //nodeBox.shift(-a_nGhost);
      FArrayBox nodeFAB;
      a_uFab.define(nodeBox,nodeVelNcomp,a_uVelPtr);
      a_vFab.define(nodeBox, nodeVelNcomp,a_vVelPtr);
      if (a_fillData)
        {
          int yGhost = a_nGhost[1];
          int yoffset = a_offset[1];
          FORT_NODETOCELLCISMVEL(CHF_CONST_FRA(a_uFab),
                                 CHF_FRA(a_uccFab),
                                 CHF_CONST_INT(yGhost),
                                 CHF_CONST_INT(yoffset),
                                 CHF_BOX(fabBox));
          
          FORT_NODETOCELLCISMVEL(CHF_CONST_FRA(a_vFab),
                                 CHF_FRA(a_vccFab),
                                 CHF_CONST_INT(yGhost),
                                 CHF_CONST_INT(yoffset),
                                 CHF_BOX(fabBox));
        }
    }
  else 
    {
      // both of these are aliased to point to the input data
      a_uFab.define(CISMvelBox, velNcomp, a_uVelPtr);
      a_vFab.define(CISMvelBox, velNcomp, a_vVelPtr);
      a_uccFab.define(CISMvelBox, velNcomp, a_uVelPtr);
      a_vccFab.define(CISMvelBox, velNcomp, a_vVelPtr);
    }
}


void
FortranInterfaceIBC::setThickness(Real* a_data_ptr,
				  const int* a_dimInfo,
                                  const int* a_boxlo, const int* a_boxhi, 
                                  const Real* a_dew, const Real* a_dns,
				  const IntVect& a_offset,
                                  const IntVect& a_nGhost,
				  const bool a_nodal)
{

  m_thicknessGhost = a_nGhost;
  m_nodalThickness = a_nodal;

  // dimInfo is (SPACEDIM, nz, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering
  //cout << "a_dimonfo" << a_dimInfo[0] << a_dimInfo[1] << endl;  

  if (m_verbose)
    {
      pout() << "In FortranInterfaceIBC::setThickness:" << endl;
      pout() << " -- entering setFAB..." << endl;
    }

  setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi,
         a_dew,a_dns,a_offset,a_nGhost,
	 m_inputThickness, m_ccInputThickness, a_nodal, m_verbose);

  if (m_verbose)
    {
      pout() << "... done" << endl;
    }

  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(m_ccInputThickness.box());
      gridBox.grow(-a_nGhost);
      setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }
  

  m_inputThicknessDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an  aliased 
  // constructor for the LevelData, but this should be fine for now....
  
  // if nodal, we'd like at least one ghost cell for the LDF
  // (since we'll eventually have to average back to nodes)
  IntVect LDFghost = m_thicknessGhost;
  if (a_nodal && (LDFghost[0] == 0))
    {
      LDFghost += IntVect::Unit;
    }
      
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, LDFghost) );

  m_inputThicknessLDF = localLDFPtr;
  // fundamental assumption that there is no more than one box/ processor 
  // don't do anything if there is no data on this processor
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box copyBox = (*m_inputThicknessLDF)[dit].box();
      copyBox &= m_inputThickness.box();
      (*m_inputThicknessLDF)[dit].copy(m_inputThickness, copyBox);
      
      // if necessary, set thickness to zero where needed
      if (m_thicknessClearRegions.size() > 0) 
        {
          for (int i=0; i<m_thicknessClearRegions.size(); i++)
            {
              Box region = m_thicknessClearRegions[i];
              // adjust for ghosting
              region.shift(-a_nGhost);
              region &= m_inputThickness.box();
              if (!region.isEmpty())
                {
                  // do this in distributed copy, rather than in the original
                  (*m_inputThicknessLDF)[dit].setVal(0.0, region, 0, 1);
                }
            } // end loop over thickness clear regions
        }
    } // end DataIterator loop

}

void
FortranInterfaceIBC::setTopography(Real* a_data_ptr,
				   const int* a_dimInfo,
                                   const int* a_boxlo, const int* a_boxhi, 
				   const Real* a_dew, const Real* a_dns,
				   const IntVect& a_offset,
                                   const IntVect& a_nGhost, 
				   const bool a_nodal)
{

  // dimInfo is (SPACEDIM, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering

  m_topographyGhost = a_nGhost;
  m_nodalTopography = a_nodal;
  setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi, a_dew,a_dns,
         a_offset, a_nGhost,m_inputTopography, m_ccInputTopography,
         a_nodal, m_verbose);

  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(m_ccInputTopography.box());
      gridBox.grow(-a_nGhost);
      setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }

  m_inputTopographyDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an aliased 
  // constructor for the LevelData, but this should be fine for now....

  // if nodal, we'd like at least one ghost cell for the LDF
  // (since we'll eventually have to average back to nodes)
  
  IntVect LDFghost = m_topographyGhost;
  if (a_nodal && (LDFghost[0] == 0) )
    {
      LDFghost += IntVect::Unit;
    }
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, LDFghost));
  m_inputTopographyLDF = localLDFPtr;
  // fundamental assumption that there is no more than one box per processor/
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box copyBox = (*m_inputTopographyLDF)[dit].box();
      copyBox &= m_inputTopography.box();
      (*m_inputTopographyLDF)[dit].copy(m_inputTopography, copyBox);
    }

}


void
FortranInterfaceIBC::setSurface(Real* a_data_ptr,
                                const int* a_dimInfo,
                                const int* a_boxlo, const int* a_boxhi, 
                                const Real* a_dew, const Real* a_dns,
                                const IntVect& a_offset,
                                const IntVect& a_nGhost, 
                                const bool a_nodal)
{

  // dimInfo is (SPACEDIM, nx, ny)

  // assumption is that data_ptr is indexed using fortran 
  // ordering from (1:dimInfo[1])1,dimInfo[2])
  // we want to use c ordering

  m_surfaceGhost = a_nGhost;
  m_nodalSurface = a_nodal;
  setFAB(a_data_ptr, a_dimInfo,a_boxlo, a_boxhi, a_dew,a_dns,
         a_offset, a_nGhost,m_inputSurface, m_ccInputSurface,
         a_nodal, m_verbose);

  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(m_ccInputSurface.box());
      gridBox.grow(-a_nGhost);
      setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }

  m_inputSurfaceDx = RealVect(D_DECL(*a_dew, *a_dns, 1));

  // now define LevelData and copy from FAB->LevelData 
  // (at some point will likely change this to be an aliased 
  // constructor for the LevelData, but this should be fine for now....

  // if nodal, we'd like at least one ghost cell for the LDF
  // (since we'll eventually have to average back to nodes)
  
  IntVect LDFghost = m_surfaceGhost;
  if (a_nodal && (LDFghost[0] == 0) )
    {
      LDFghost += IntVect::Unit;
    }
  RefCountedPtr<LevelData<FArrayBox> > localLDFPtr(new LevelData<FArrayBox>(m_grids, 1, LDFghost));
  m_inputSurfaceLDF = localLDFPtr;
  // fundamental assumption that there is no more than one box per processor/
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box copyBox = (*m_inputSurfaceLDF)[dit].box();
      copyBox &= m_inputSurface.box();
      (*m_inputSurfaceLDF)[dit].copy(m_inputSurface, copyBox);
    }

}


/// regions where we artificially set thickness to zero
/** this is done in setThickness, for lack of a better place, so 
    this needs to be set before setThickness is called. 
    a_clearRegions defines logically-rectangular regions where 
    the thickness is artificially set to zero. These regions are 
    defined relative to the original input data (i.e. before any shifting 
    due to ghost cells)
*/
void
FortranInterfaceIBC::setThicknessClearRegions(const Vector<Box>& a_clearRegions)
{
  m_thicknessClearRegions = a_clearRegions;
}


/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
IceThicknessIBC* 
FortranInterfaceIBC::new_thicknessIBC()
{
  if (m_verbose)
    {
      pout() << "in FortranInterfaceIBC::new_thicknessIBC" << endl;
    }

  FortranInterfaceIBC* retval = new FortranInterfaceIBC();

  retval->m_grids = m_grids;
  retval->m_gridsSet = m_gridsSet;
  
  // keep these as aliases, if they're actually defined
  
  if (!m_inputThickness.box().isEmpty())
    {
      retval->m_inputThickness.define(m_inputThickness.box(), 
                                      m_inputThickness.nComp(),
                                      m_inputThickness.dataPtr());
    }

  if (!m_ccInputThickness.box().isEmpty())
    {
      retval->m_ccInputThickness.define(m_ccInputThickness.box(), 
                                        m_ccInputThickness.nComp(),
                                        m_ccInputThickness.dataPtr());
    }      

  retval->m_inputThicknessLDF = m_inputThicknessLDF;
  
  retval->m_thicknessGhost = m_thicknessGhost;
  retval->m_inputThicknessDx = m_inputThicknessDx;
  
  // keep these as aliases
  if (!m_inputTopography.box().isEmpty())
    {
      retval->m_inputTopography.define(m_inputTopography.box(), m_inputTopography.nComp(), m_inputTopography.dataPtr());
    }

  if (!m_ccInputTopography.box().isEmpty())
    {
      retval->m_ccInputTopography.define(m_ccInputTopography.box(), m_ccInputTopography.nComp(), m_ccInputTopography.dataPtr());
    }

  retval->m_inputTopographyLDF = m_inputTopographyLDF;

  retval->m_inputTopographyDx = m_inputTopographyDx;
  retval->m_topographyGhost = m_topographyGhost;

  retval->m_verbose = m_verbose;

  return static_cast<IceThicknessIBC*>(retval);

}

/// Set up initial conditions
/**
 */
void
FortranInterfaceIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be here...
  MayDay::Error("FortranInterfaceIBC::initialize not implemented");
}




/// Set boundary fluxes
/**
 */
void 
FortranInterfaceIBC::primBC(FArrayBox&            a_WGdnv,
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
	    //a_WGdnv(i,0) = 0.0;
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
FortranInterfaceIBC::setBdrySlopes(FArrayBox&       a_dW,
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
FortranInterfaceIBC::artViscBC(FArrayBox&       a_F,
                             const FArrayBox& a_U,
                             const FArrayBox& a_divVel,
                             const int&       a_dir,
                             const Real&      a_time)
{
  // don't anticipate being here -- if we wind up here, need to
  // give it some thought
  MayDay::Error("FortranInterfaceIBC::artViscBC not implemented");
}


/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
FortranInterfaceIBC::velocitySolveBC()
{
  
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  
  return m_velBCs;
}


/// set non-periodic ghost cells for surface height z_s. 
void
FortranInterfaceIBC::setSurfaceHeightBCs(LevelData<FArrayBox>& a_zSurface,
                                         LevelSigmaCS& a_coords,
                                         const ProblemDomain& a_domain,
                                         const RealVect& a_dx,
                                         Real a_time, Real a_dt)
{
  if (m_extrapBoundary)
    {
      ExtrapGhostCells(a_zSurface, a_domain);
    } 
  else
    {
      //slc : i think reflection is a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      //      (and in that case, this ought to be unnecessary)
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_zSurface, a_domain, dir, Side::Hi);
	}
    }
}

/// set non-periodic ghost cells for thickness & topography
void
FortranInterfaceIBC::setGeometryBCs(LevelSigmaCS& a_coords,
				    const ProblemDomain& a_domain,
				    const RealVect& a_dx,
				    Real a_time, Real a_dt)
{
  if (m_extrapBoundary)
    {
      ExtrapGhostCells(a_coords.getH(), a_domain);
      ExtrapGhostCells(a_coords.getTopography(), a_domain);
    } 
  else
    {
      //slc : i think reflection is a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      //      (and in that case, this ought to be unnecessary)
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_coords.getH(), a_domain, dir, Side::Hi);
	  ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Lo);
	  ReflectGhostCells(a_coords.getTopography(), a_domain, dir, Side::Hi);
	}
    }

}


/// set grids using Boxes passed in from Glimmer-CISM
/** creates a DisjointBoxLayout using the grid boxes and 
    processor distribution used by CISM. */
void 
FortranInterfaceIBC::setGrids(DisjointBoxLayout& a_grids, 
                              const Box& a_gridBox,
                              const ProblemDomain& a_domain,
                              bool a_verbose)
{

  if (a_verbose)
    {
      pout() << "in setGrids, box = " << a_gridBox << endl;
    }

  int status;
  Vector<Box> filteredBoxes;
  Vector<int> filteredProcAssign;

  status = LoadBalance(filteredProcAssign,
                       filteredBoxes,
                       a_gridBox,
                       numProc());
  CH_assert(status == 0);
  // define DisjointBoxLayout
  a_grids.define(filteredBoxes, filteredProcAssign, a_domain);
  
  if (a_verbose)
    {
      pout () << "after DBL define" << endl;
    }
}

/// utility function to fill in holes in topography
/** looks for isolated values of holeVal and replaces then
    with average of neighbors */
void
FortranInterfaceIBC::fillTopographyHoles(Real a_holeVal)
{

  if (m_verbose) 
    {
      pout() << "FortranInterfaceIBC::fillTopographyHoles" << endl;
    }

  CH_assert(!m_inputTopography.box().isEmpty());
  
  Real compareThreshold = 1.0e-10;
  // neighborCutoff is number of neighbors which are not the 
  // holeValue required before we do averaging. (max number of 
  // neigbors in 2d is 8
  int neighborCutoff = 5;

  Box testBox = m_inputTopography.box();
  testBox.grow(-1);

  Box neighborBox(-IntVect::Unit, IntVect::Unit);
  IntVectSet neighborIVS(neighborBox);
  neighborIVS -= IntVect::Zero;

  IVSIterator neighborIt(neighborIVS);

  BoxIterator bit(testBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (abs(m_inputTopography(iv, 0) - a_holeVal) < compareThreshold)
        {
          Real sum = 0;
          int num = 0;
          for (neighborIt.begin(); neighborIt.ok(); ++neighborIt)
            {
              IntVect newIV = iv + neighborIt();
              if (abs(m_inputTopography(newIV, 0) - a_holeVal) > compareThreshold)
                {
                  sum += m_inputTopography(newIV, 0);
                  num++;
                }
            }

          if (num > neighborCutoff)
            {
              m_inputTopography(iv,0) = sum / num;
            }
        }
    }
          

#if 0
  // try this in c++ first
  FORT_FILLHOLES(CHF_FRA1(m_inputTopography,0),
                 CHF_REAL(a_holeVal),
                 CHF_REAL(compareThreshold),
                 CHF_INT(neighborCutoff),
                 CHF_BOX(testBox));
#endif
}


  /// set up initial ice state
  /** reads info from ParmParse and sets up ice sheet geometry
   */
void
FortranInterfaceIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                           const RealVect& a_dx,
                                           const RealVect& a_domainSize,
                                           const Real& a_time, 
					   const LevelSigmaCS* a_crseCoords,
					   const int a_refRatio)
{
  if (m_verbose) 
    {
      pout() << "FortranInterfaceIBC::initializeIceGeometry" << endl;
    }

  
  ParmParse geomPP("geometry");

  //CH_assert(m_isDefined);

  // assume that thickness and topography from glimmer are defined on the same 
  // mesh
  CH_assert(m_inputThickness.box() == m_inputTopography.box());
  CH_assert(m_inputThicknessDx == m_inputTopographyDx);

  // get scaling factor for ice thickness -- not sure if I really need this
  Real thicknessScale =1;
  geomPP.query("thickness_scale", thicknessScale);

  if (m_verbose)
    {
      pout() << " ...setting up thickness LevelData..." << endl;
    }

  
  // grab periodicity info from a_coords
  const DisjointBoxLayout& destGrids = a_coords.grids();
  const ProblemDomain& destDomain = destGrids.physDomain();
  if (geomPP.contains("basalSlope") 
      && (destDomain.isPeriodic(0) || destDomain.isPeriodic(1)) )
    {
      //slc : periodic domains need not have periodic topography b,
      //instead b = f(x,y) + Ax + By where A,B are constants and f is periodic.
      //for now, read A,B from the input file - periodic problems
      //are rare, so it is possibly not worth the effort of devising
      //something robust, e.g working out A and B from glimmer's model%geometry%topg, 
      //which is node-centred data. It would be easy to do that in serial, 
      //but in parallel the details would depend on the
      //interaction between glimmer's domain decomposition and ours.
      RealVect basalSlope(RealVect::Zero);
      Vector<Real> t(SpaceDim, 0.0);
      geomPP.getarr("basalSlope", t, 0, SpaceDim);
      D_TERM(basalSlope[0] = t[0];,
	     basalSlope[1] = t[1];,
	     basalSlope[2] = t[2];);
      if (m_verbose)
	{      
	  pout() << " ... adding constant background slope " 
		 << basalSlope ;
	}
      if (basalSlope.dotProduct(basalSlope) > 1.0e-10)
	{
	  MayDay::Error("periodic problems with non-zero basal slope not yet implemented");
	}
      a_coords.setBackgroundSlope(basalSlope);

      if (m_verbose)
	{      
	  pout() << " a_domainSize[0] = "  
		 << a_domainSize[0]
		 << " ... " << std::endl;
	}

      
    }

  if (m_verbose)
    {
      pout() << " Done with background slope" << endl;
    }

  DisjointBoxLayout ChomboGrids = a_coords.grids();
  LevelData<FArrayBox>& ChomboThickness = a_coords.getH();

  IntVect thicknessGhostVect =  ChomboThickness.ghostVect();
  IntVect sigmaGhostVect = thicknessGhostVect - IntVect::Unit;

  //LevelData<FArrayBox> ChomboTopography(ChomboGrids, 1, thicknessGhostVect);
  LevelData<FArrayBox>& ChomboTopography = a_coords.getTopography();
  Real tolerance = 1.0e-6;
  if (a_crseCoords != NULL &&
      a_refRatio * a_dx[0] <= m_inputThicknessDx[0] * (1.0 + tolerance))
    {
      // in this (common) case, interpolation from a_crseCoords is as good as it gets
      if (m_verbose)
	{
	  pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		 << a_refRatio << endl;
	}

      a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
    }
  else 
    {
      if (m_verbose)
        {
          pout() << "calling FillFromReference for thickness..." << endl;
          
        }

      FillFromReference(ChomboThickness,
			(*m_inputThicknessLDF),
			a_dx,
			m_inputThicknessDx,
			m_verbose);


      if (m_verbose)
        {
          pout() << "...and topography..." << endl;
        }

      FillFromReference(ChomboTopography,
			(*m_inputTopographyLDF),
			a_dx,
			m_inputTopographyDx,
			m_verbose);
      if (m_verbose)
        {
          pout() << "...done!" << endl;
        }

      if (a_crseCoords!= NULL)
	{
          if (m_verbose)
            {
              pout() << "filling ghost cells" << endl;
            }
	  // fill ghost cells
	  a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				    false, false, false);
          if (m_verbose) 
            {
              pout() << "done filling ghost cells" << endl;
            }
	}

    }
  
  if (m_verbose) 
    {      
      pout() << "calling exchange..." << endl;
    }
  ChomboThickness.exchange();
  //ChomboTopography.exchange();
  //LevelSigmaCS.exchangeTopography() takes care of constant topography slopes
  a_coords.exchangeTopography();
  
  // if we have ghost cells, then use extrapolation BC's (for now, at least)
  // to fill in reasonable values for topography and thickness
  m_extrapBoundary = true;
  geomPP.query("extrap_boundary", m_extrapBoundary );
  if (m_extrapBoundary)
    {
      if (m_verbose)
	{      
	  pout() << "extrapolating boundary thickness and topography" << endl;
	}
      ExtrapGhostCells(ChomboThickness, destDomain);
      ExtrapGhostCells(ChomboTopography, destDomain);
    }
  else 
    {
      //slc : I think reflection would make a better default : works for isolated
      //      islands (Antarctica, Greenland) and divides, 
      //      where grad(s) = grad(H) = grad(topo) = 0 (Pine Island).
      if (m_verbose)
	{      
	  pout() << "reflecting boundary thickness and topography" << endl;
	}
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  ReflectGhostCells(ChomboThickness, destDomain, dir, Side::Lo);
	  ReflectGhostCells(ChomboThickness, destDomain, dir, Side::Hi);
	  ReflectGhostCells(ChomboTopography, destDomain, dir,  Side::Lo);
	  ReflectGhostCells(ChomboTopography, destDomain, dir, Side::Hi);
	}
    }
  //this can wait
  //a_coords.recomputeGeometry();
  if (m_verbose)
    {      
      pout() << "leaving FortranInterfaceIBC::initializeIceGeometry" << endl;
    }
}

bool
FortranInterfaceIBC::regridIceGeometry(LevelSigmaCS& a_coords,
				       const RealVect& a_dx,
				       const RealVect& a_domainSize,
				       const Real& a_time, 
				       const LevelSigmaCS* a_crseCoords,
				       const int a_refRatio)
{

   if (m_verbose)
    {      
      pout() << "entering FortranInterfaceIBC::regridIceGeometry" << endl;
    }

   Real tolerance = 1.0e-6;

   if (a_dx[0] + tolerance < m_inputThicknessDx[0] )
    return false; // if the requested grid is finer than the stored DEM, the best approach is interpolation

   if (a_crseCoords != NULL &&
       a_refRatio * a_dx[0] <= m_inputThicknessDx[0] * (1.0 + tolerance))
     {
       // in this (common) case, interpolation from a_crseCoords is as good as it gets
       // now deprecated
       MayDay::Error("FortranInterfaceIBC::regridIceGeometry interpolation from a_crseCoords deprecated");
       if (m_verbose)
	 {
	   pout() << " ...interpolating data from coarse LevelSigmaCS with refinement ratio = " 
		  << a_refRatio << endl;
	 }
       
       a_coords.interpFromCoarse(*a_crseCoords, a_refRatio);
     }
   else
     {
       LevelData<FArrayBox>& ChomboTopography = a_coords.getTopography();
       FillFromReference(ChomboTopography,
			 m_inputTopography,
			 a_dx,
			 m_inputTopographyDx,
			 m_topographyGhost,
			 m_verbose);
       if (a_crseCoords!= NULL)
	 {
	   // fill ghost cells
	   a_coords.interpFromCoarse(*a_crseCoords, a_refRatio, 
				     false, false, false);
	 }
     }
   if (m_verbose)
     {      
       pout() << "leaving FortranInterfaceIBC::regridIceGeometry" << endl;
     }
   
   // we shouldn't ever be in here anyway on this branch
   return false;
}

/*
  flatten thickness and basal topography back to input FArrayBoxes
  re-store in m_inputThickness and m_inputTopography -- not const 
  because we modify the values in the data holders
*/
void
FortranInterfaceIBC::flattenIceGeometry(const Vector<RefCountedPtr<LevelSigmaCS> > & a_amrGeometry)
{
  
  // start with coarsest level and call FillFromReference
  int numLevels = a_amrGeometry.size();

  Vector<LevelData<FArrayBox>* > vectH(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > vectTopo(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > vectSurface(numLevels, NULL);
  Vector<RealVect> vectDx(numLevels);
  for (int lev=0; lev<numLevels; lev++)
    {
      const LevelData<FArrayBox>& levelThickness = a_amrGeometry[lev]->getH();
      const LevelData<FArrayBox>& levelTopography = a_amrGeometry[lev]->getTopography();
      const LevelData<FArrayBox>& levelSurface = a_amrGeometry[lev]->getSurfaceHeight();
      RealVect levelDx = a_amrGeometry[lev]->dx();

      // can cast away const-ness here because they will be sent in as 
      // const arguments
      vectH[lev] = const_cast<LevelData<FArrayBox>* >(&levelThickness);
      vectTopo[lev] = const_cast<LevelData<FArrayBox>* >(&levelTopography);
      vectSurface[lev] = const_cast<LevelData<FArrayBox>* >(&levelSurface);
      vectDx[lev] = levelDx;

      // these can be useful for debugging..
      bool resetInitialLevel = false;
      bool resetEachLevel = false;
      if (resetEachLevel || (resetInitialLevel && (lev == 0)))
        {
          DataIterator dit = m_inputThicknessLDF->dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              (*m_inputThicknessLDF)[dit].setVal(float(lev*111));
            }
        }

    } // end loop over levels

  if (!m_inputThicknessLDF.isNull())
    {
      flattenCellData(*m_inputThicknessLDF, m_inputThicknessDx,
                      vectH, vectDx, m_verbose);
      m_inputThicknessLDF->exchange();
    }


  if (!m_inputTopographyLDF.isNull())
    {
      flattenCellData(*m_inputTopographyLDF, m_inputTopographyDx,
                      vectTopo, vectDx, m_verbose);
      m_inputTopographyLDF->exchange();
    }
  
  if (!m_inputSurfaceLDF.isNull())
    {
      flattenCellData(*m_inputSurfaceLDF, m_inputSurfaceDx,
                      vectSurface, vectDx, m_verbose);
      m_inputSurfaceLDF->exchange();
    }

  
  if (!m_inputThicknessLDF.isNull())
    {
      // now copy back to original storage
      if (m_nodalThickness)
        {
          // have to average back to nodes
          DataIterator dit= m_inputTopographyLDF->dataIterator();
          //const DisjointBoxLayout& grids = m_inputTopographyLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              
              const Box nodeBox = m_inputThickness.box();
              Box grownCCBox = nodeBox;
              grownCCBox.enclosedCells();
              grownCCBox.grow(1);
              
              FArrayBox ccGrownFAB(grownCCBox,1);
              ccGrownFAB.copy((*m_inputTopographyLDF)[dit]);
              
              // Cell->Node looks a lot like Node->Cell...
              FORT_NODETOCELL(CHF_CONST_FRA1(ccGrownFAB,0),
                              CHF_FRA1(m_inputThickness,0),
                              CHF_BOX(nodeBox));
              
              // don't think I really need to do anything with 
              // m_ccInputThickness 
              
            } // end dataIterator "loop"
          
        } // end if nodal thickness    
      else
        {
          // simple copy will do
          DataIterator dit= m_inputThicknessLDF->dataIterator();
          const DisjointBoxLayout& grids = m_inputThicknessLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              
              Box intersectBox(grids[dit]);
              intersectBox &= m_ccInputThickness.box();
              
              m_ccInputThickness.copy((*m_inputThicknessLDF)[dit], intersectBox);
              
            } // end dataIterator "loop"
        }
    } // end if thickness was defined 
  

  // now copy back to original storage
  if (!m_inputTopographyLDF.isNull())
    {
      if (m_nodalTopography)
        {
          // have to average back to nodes
          DataIterator dit= m_inputTopographyLDF->dataIterator();
          //const DisjointBoxLayout& grids = m_inputTopographyLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box nodeBox = m_inputTopography.box();
              Box grownCCBox = nodeBox;
              grownCCBox.enclosedCells();
              grownCCBox.grow(1);
              
              FArrayBox ccGrownFAB(grownCCBox,1);
              ccGrownFAB.copy((*m_inputTopographyLDF)[dit]);
              
              // Cell->Node looks a lot like Node->Cell...
              FORT_NODETOCELL(CHF_CONST_FRA1(ccGrownFAB,0),
                              CHF_FRA1(m_inputTopography,0),
                              CHF_BOX(nodeBox));
              
              // don't think I really need to do anything with 
              // m_ccInputTopography 
              
            } // end dataIterator "loop"
          
        } // end if nodal thickness    
      else
        {
          // simple copy will do
          DataIterator dit= m_inputTopographyLDF->dataIterator();
          const DisjointBoxLayout& grids = m_inputTopographyLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              
              Box intersectBox(grids[dit]);
              intersectBox &= m_ccInputTopography.box();
              
              m_ccInputTopography.copy((*m_inputTopographyLDF)[dit], intersectBox);
              
            } // end dataIterator "loop"
        }      
    } // end if input topography is defined
  
  // now copy back to original storage

  if (!m_inputSurfaceLDF.isNull())
    {
      if (m_nodalSurface)
        {
          // have to average back to nodes
          DataIterator dit= m_inputSurfaceLDF->dataIterator();
          //const DisjointBoxLayout& grids = m_inputTopographyLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box nodeBox = m_inputSurface.box();
              Box grownCCBox = nodeBox;
              grownCCBox.enclosedCells();
              grownCCBox.grow(1);
              
              FArrayBox ccGrownFAB(grownCCBox,1);
              ccGrownFAB.copy((*m_inputSurfaceLDF)[dit]);
              
              // Cell->Node looks a lot like Node->Cell...
              FORT_NODETOCELL(CHF_CONST_FRA1(ccGrownFAB,0),
                              CHF_FRA1(m_inputSurface,0),
                              CHF_BOX(nodeBox));
              
              // don't think I really need to do anything with 
              // m_ccInputTopography 
              
            } // end dataIterator "loop"
          
        } // end if nodal surface
      else
        {
          // simple copy will do
          DataIterator dit= m_inputSurfaceLDF->dataIterator();
          const DisjointBoxLayout& grids = m_inputSurfaceLDF->getBoxes();
          for (dit.begin(); dit.ok(); ++dit)
            {
              
              Box intersectBox(grids[dit]);
              intersectBox &= m_ccInputSurface.box();
              
              m_ccInputSurface.copy((*m_inputSurfaceLDF)[dit], intersectBox);
              
            } // end dataIterator "loop"
        }
    } // end if surface has been defined

}


/// flatten an arbitrary dataset back to input FArrayBoxes
void
FortranInterfaceIBC::flattenData(Real* a_data_ptr,
                                 const int* a_dimInfo,
                                 const int* a_boxlo, const int* a_boxhi,
                                 const Real* a_dew, const Real* a_dns,
                                 const IntVect& a_offset,
                                 const Vector<LevelData<FArrayBox>* >& a_amrData,
                                 const Vector<int>& a_vectRefRatio,
                                 const Vector<Real>& a_amrDx,
                                 int a_srcComp,
                                 int a_destComp,
                                 int a_nComp,
                                 const IntVect& a_nGhost,
                                 const bool a_nodal)
{
  RealVect dx(D_DECL(*a_dew,*a_dns,*a_dew) );
              
  // first set up FArrayBox
  FArrayBox dataFab, CCdataFab;

  setFAB(a_data_ptr, a_dimInfo, a_boxlo, a_boxhi,
         a_dew, a_dns, a_offset, a_nGhost,
         dataFab, CCdataFab, a_nodal, m_verbose);
        
  // if we haven't already set the grids, do it now
  if (!gridsSet())
    {
      if (m_verbose) 
        {
          pout() << " -- entering setGrids" << endl;
        }
      Box gridBox(CCdataFab.box());
      gridBox.grow(-a_nGhost);
      setGrids(m_grids, gridBox, m_domain, m_verbose);
      m_gridsSet = true;
      if (m_verbose)
        {
          pout() << " -- out of setGrids" << endl;
        }
    }
      
  // create LevelData
  // if nodal, we'd like at least one ghost cell for the LDF
  // (since we'll eventually have to average back to nodes)
  IntVect LDFghost = m_thicknessGhost;
  if (a_nodal && (LDFghost[0] == 0))
    {
      LDFghost += IntVect::Unit;
    }
  LevelData<FArrayBox> ldf(m_grids, a_nComp, LDFghost);
              
  // start with coarsest level, calling FillFromReference
  int numLevels = a_amrData.size();
  for (int lev=0; lev<numLevels; lev++)
    {
      RealVect levelDx = a_amrDx[lev]*RealVect::Unit;
      FillFromReference(ldf, *a_amrData[lev],
                        dx, levelDx,
                        m_verbose);
    }
    
  ldf.exchange();

  // finally, either copy or average to nodes
  if (a_nodal)
    {
      // have to average back to nodes
      DataIterator dit= ldf.dataIterator();
      //const DisjointBoxLayout& grids = ldf.getBoxes();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box nodeBox = dataFab.box();
          Box grownCCBox = nodeBox;
          grownCCBox.enclosedCells();
          grownCCBox.grow(1);
          
          FArrayBox ccGrownFAB(grownCCBox,1);
          ccGrownFAB.copy(ldf[dit]);

          FORT_CELLTONODE(CHF_FRA1(dataFab,0),
                          CHF_CONST_FRA1(ccGrownFAB,0),
                          CHF_BOX(nodeBox));
          
          
        } // end dataIterator "loop"
      
    } // end if nodal thickness    
  else
    {
      // simple copy will do
      DataIterator dit= ldf.dataIterator();
      const DisjointBoxLayout& grids = ldf.getBoxes();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box intersectBox(grids[dit]);
          intersectBox &= CCdataFab.box();
                           
          CCdataFab.copy(ldf[dit], intersectBox);
          
        } // end dataIterator "loop"
    } // end if not nodal

}

/// flatten velocity field to the 3D velocity data holder
void
FortranInterfaceIBC::flattenVelocity(Real* a_uVelPtr, Real* a_vVelPtr,
                                     const int* a_dimInfo,
                                     const int* a_boxlo, const int* a_boxhi,
                                     const Real* a_dew, const Real* a_dns,
                                     const IntVect& a_offset,
                                     const Vector<LevelData<FArrayBox>* >& a_amrVel,
                                     const Vector<int>& a_vectRefRatio,
                                     const Vector<Real>& a_amrDx,
                                     const IntVect& a_nGhost,
                                     const bool a_nodal)
{

  RealVect dx(D_DECL(*a_dew,*a_dns,*a_dew) );
              
  // first set up FArrayBox
  FArrayBox uDataFab, vDataFab, uCCdataFab, vCCdataFab;
  
  setVelFAB(a_uVelPtr, a_vVelPtr,a_dimInfo, a_boxlo, a_boxhi,
            a_dew, a_dns, a_offset, a_nGhost,
            uDataFab, vDataFab, uCCdataFab, vCCdataFab, a_nodal,
            false);

  
  if (m_verbose)
    {
      pout() << "exited setVelFAB function: uDataFAB box: " 
             << uDataFab.box() << ", ncomp = " << uDataFab.nComp() << endl;
      pout() << "                           vDataFAB box: "
             << vDataFab.box() << ", ncomp = " << vDataFab.nComp() << endl;
    }

  // create LevelData
  // if nodal, we'd like one more ghost cell for the LDF than nghost
  // (since we'll need to average back to nodes)
  IntVect LDFghost = a_nGhost;
  LDFghost += IntVect::Unit;

  LevelData<FArrayBox> ldf(m_grids, SpaceDim, LDFghost);

  // need to create Vector<RealVect> form of dx
  Vector<RealVect> realVectDx(a_amrDx.size());
  for (int i=0; i<a_amrDx.size(); i++)
    {
      realVectDx[i] = a_amrDx[i]*RealVect::Unit;
    }

  MayDay::Error("FortranInterfaceIBC::flattenVelocity not fully implemented");
  //flattenCellData(ldf, dx, a_amrVel, realVectDx, m_verbose);
  
  if (m_verbose)  pout () << "before exchange" << endl;
  ldf.exchange();

  //  writeLevel(&ldf);
  if (m_verbose) pout() << "after exchange" << endl;


#if 0
  // debugging aid -- set velocity to a constant value
  {
    Real uVal = 5.0;
    Real vVal = 10.0;
    DataIterator dit = ldf.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& velFab = ldf[dit];
        Box gridBox = m_grids[dit];
        velFab.setVal(uVal, gridBox, 0);
        velFab.setVal(vVal, gridBox, 1);
      }
  }
#endif

  // finally, either copy or average to nodes
  if (a_nodal)
    {
      // have to average back to nodes
      DataIterator dit= ldf.dataIterator();
      const DisjointBoxLayout& grids = ldf.getBoxes();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box nodeBox = grids[dit];
          nodeBox.surroundingNodes();
          nodeBox.grow(a_nGhost);

          // first do u...
          int xghost = a_nGhost[0];
          int xoffset = a_offset[0];
          
          if (m_verbose) pout () << "entering cell-to-node" << endl;
          
          IntVect boxloVect(D_DECL(a_boxlo[0], a_boxlo[1], a_boxlo[2]));

          FORT_CELLTONODECISMVELNOSHEAR(CHF_FRA(uDataFab),
                                        CHF_CONST_FRA1(ldf[dit],0),
                                        CHF_CONST_INT(xghost),
                                        CHF_CONST_INT(xoffset),
                                        CHF_INTVECT(boxloVect),
                                        CHF_BOX(nodeBox));

          if (SpaceDim > 1)
            {
              // ...then do v
              int yghost = a_nGhost[1];
              int yoffset = a_offset[1];

              FORT_CELLTONODECISMVELNOSHEAR(CHF_FRA(vDataFab),
                                            CHF_CONST_FRA1(ldf[dit],1),
                                            CHF_CONST_INT(yghost),
                                            CHF_CONST_INT(yoffset),
                                            CHF_INTVECT(boxloVect),
                                            CHF_BOX(nodeBox));
            }

        } // end dataIterator "loop"
      
    } // end if nodal 
  else
    {
      // simple copy will do
      DataIterator dit= ldf.dataIterator();
      const DisjointBoxLayout& grids = ldf.getBoxes();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box intersectBox(grids[dit]);
          intersectBox &= uCCdataFab.box();
                           
          // first copy u...
          uCCdataFab.copy(ldf[dit], intersectBox, 0, 
                          intersectBox, 0, 1);

          // ... then v
          vCCdataFab.copy(ldf[dit], intersectBox, 1, 
                          intersectBox, 0, 1);
          
        } // end dataIterator "loop"
    } // end if not nodal


}


void 
FortranInterfaceIBC::setupBCs()
{
  ParmParse ppBC("bc");

  // get boundary conditions 
  Vector<int> loBCvect(SpaceDim), hiBCvect(SpaceDim);
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
          // slip-wall -- Dirichlet in normal direction, free-slip otherwise
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
          // slip-wall case -- Dirichlet in normal direction, free-slip otherwis
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
              MayDay::Error("Unknown BC type in FortranInterfaceIBC::setupBCs");
            }
        } // end loop over directions    

    } // end if using new bcs
  else
    {
      // if all BCs are the same, calll the basic BC functions
      if ((loBCvect[0] == loBCvect[1]) &&
          (hiBCvect[0] == hiBCvect[1]) &&
          (loBCvect[0] == hiBCvect[0]))
        {
          if (loBCvect[0] == FIBC_Dirichlet)
            {
              m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceDirichletBC_FIBC));
            }
          else if (loBCvect[0] == FIBC_Neumann)
            {
              m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceNeumannBC_FIBC));
            }
          else if (loBCvect[0] == FIBC_iceDivide)
            {
              m_velBCs = RefCountedPtr<BCFunction>(new BCFuncWrapper(iceDivideBC_FIBC));
            }
          else if (loBCvect[0] == FIBC_iceStreamXBC)
            {
              m_velBCs = RefCountedPtr<BCFunction>(new BCFuncWrapper(iceStreamXBC_FIBC));
            }
          else
            {
              MayDay::Error("bad BC type");
            }
        }
      else
        {
          MayDay::Error("FortranInterfaceIBC::setupBC -- bad BC type");
          
        }
    } // end if not new bc
  m_isBCsetUp = true;
}
      
void FortranInterfaceIBC::checkOK() const
{
  if (m_verbose)
    {
      pout() << "FortranInterfaceIBC::checkOK() &m_inputTopography = " 
             << &m_inputTopography << std::endl;
    }
  CH_assert(m_inputTopography.norm(0) < 1.0e+5);
}
#include "NamespaceFooter.H"
