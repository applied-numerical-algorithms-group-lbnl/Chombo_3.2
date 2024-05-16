#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
// VelocityBC.H
// ============
// Collection of  CompGridVTOBC classes

#ifndef _VELOCITYBC_H_
#define _VELOCITYBC_H_

#include "PetscCompGridVTO.H"
#include "ReflectGhostCells.H"
#include "ExtrapGhostCells.H"
#include "NamespaceHeader.H"

/// Apply no-slip conditions along one boundary
/**
   Normal velocity component = 0
   Tangent velocity component = 0
*/
class NoSlipCompGridVTOBC : public CompGridVTOBC
{

public:

  NoSlipCompGridVTOBC(int a_dir, Side::LoHiSide a_side)
  {
    m_dir = a_dir;
    m_side = a_side;
    for (int comp = 0; comp < SpaceDim; comp++)
      {
	CompGridVTOBC::setDiri(m_dir,m_side,comp,true);
      }
  }

  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    const IntVect ghostVect = IntVect::Unit;

    if (!(a_domain.isPeriodic(m_dir)))
      {
	ReflectGhostCells(a_vel, a_domain,ghostVect, m_dir, m_side);
	Box ghostBox = adjCellBox(a_valid, m_dir, m_side, 1);
	if(!a_domain.domainBox().contains(ghostBox))
	  {
	    ghostBox &= a_vel.box();
	    a_vel.mult(-1.0,ghostBox,0,SpaceDim);
	  }
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;
  NoSlipCompGridVTOBC(){;}
  NoSlipCompGridVTOBC(const NoSlipCompGridVTOBC&)
  {
    MayDay::Error(" NoSlipCompGridVTOBC(const NoSlipCompGridVTOBC&) not implemented");
  }

  NoSlipCompGridVTOBC& operator=(const NoSlipCompGridVTOBC&)
  {
    MayDay::Error(" NoSlipCompGridVTOBC& operator= (const NoSlipCompGridVTOBC&) not implemented");
    return *this;
  }
};

/// Apply ice divide conditions along one boundary
/**
   Normal velocity component = 0
   Tangent velocity component normal gradient = 0
*/
class IceDivideCompGridVTOBC : public CompGridVTOBC
{

public:

  IceDivideCompGridVTOBC(int a_dir, Side::LoHiSide a_side)
  {
    m_dir = a_dir;
    m_side = a_side;
    
    CompGridVTOBC::setDiri(m_dir,m_side,m_dir,true);
    for (int i = 1; i < SpaceDim; i++) 
      CompGridVTOBC::setDiri(m_dir,m_side,(m_dir+i)%SpaceDim,false);

  }

  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    const IntVect ghostVect = IntVect::Unit;

    if (!(a_domain.isPeriodic(m_dir)))
      {
	ReflectGhostCells(a_vel, a_domain,ghostVect, m_dir, m_side);
	Box ghostBox = adjCellBox(a_valid, m_dir, m_side, 1);
	if(!a_domain.domainBox().contains(ghostBox))
	  {
	    ghostBox &= a_vel.box();
	    a_vel.mult(-1.0,ghostBox,m_dir,1);
	  }
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;
  IceDivideCompGridVTOBC(){;}
  IceDivideCompGridVTOBC(const IceDivideCompGridVTOBC&)
  {
    MayDay::Error(" IceDivideCompGridVTOBC(const IceDivideCompGridVTOBC&) not implemented");
  }

  IceDivideCompGridVTOBC& operator=(const IceDivideCompGridVTOBC&)
  {
    MayDay::Error(" IceDivideCompGridVTOBC& operator= (const IceDivideCompGridVTOBC&) not implemented");
    return *this;
  }
};


/// Apply natural conditions along one boundary
/**
   Normal velocity component normal gradient = 0
   Tangent velocity component normal gradient = 0
*/
class NaturalCompGridVTOBC : public CompGridVTOBC
{

public:

  NaturalCompGridVTOBC(int a_dir, Side::LoHiSide a_side)
  {
    m_dir = a_dir;
    m_side = a_side;
    for (int comp = 0; comp < SpaceDim; comp++)
      {
	CompGridVTOBC::setDiri(m_dir,m_side,comp,false);
      }
  }

  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    const IntVect ghostVect = IntVect::Unit;
    
    if (!(a_domain.isPeriodic(m_dir)))
      {
	ReflectGhostCells(a_vel, a_domain,ghostVect, m_dir, m_side);
      }
  }

private:
  int m_dir;
  Side::LoHiSide m_side;
  NaturalCompGridVTOBC(){;}
  NaturalCompGridVTOBC(const NaturalCompGridVTOBC&)
  {
    MayDay::Error(" NaturalCompGridVTOBC(const NaturalCompGridVTOBC&) not implemented");
  }

  NaturalCompGridVTOBC& operator=(const NaturalCompGridVTOBC&)
  {
    MayDay::Error(" NaturalCompGridVTOBC& operator= (const NaturalCompGridVTOBC&) not implemented");
    return *this;
  }
};


/// Apply a different CompGridVTOBC along each boundary 
class PerEdgeCompGridVTOBC : public CompGridVTOBC
{
  Vector<RefCountedPtr<CompGridVTOBC> > m_loBC;
  Vector<RefCountedPtr<CompGridVTOBC> > m_hiBC;

public:

  PerEdgeCompGridVTOBC(Vector<RefCountedPtr<CompGridVTOBC> > a_loBC,
		       Vector<RefCountedPtr<CompGridVTOBC> > a_hiBC)
    :m_loBC(a_loBC), m_hiBC(a_hiBC)
  {
    CH_assert(m_loBC.size() == SpaceDim);
    CH_assert(m_hiBC.size() == SpaceDim);
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	for (int comp = 0; comp < SpaceDim; comp++)
	  {	
	    CompGridVTOBC::setDiri(dir,Side::Lo,comp,a_loBC[dir]->isDiri(dir,Side::Lo,comp));
	    CompGridVTOBC::setDiri(dir,Side::Hi,comp,a_hiBC[dir]->isDiri(dir,Side::Hi,comp));
	  }
      }
  }

  PerEdgeCompGridVTOBC()
  {
  }

  PerEdgeCompGridVTOBC(const PerEdgeCompGridVTOBC&)
  {
    MayDay::Error(" PerEdgeCompGridVTOBC(const PerEdgeCompGridVTOBC&) not implemented");
  }

  PerEdgeCompGridVTOBC& operator=(const PerEdgeCompGridVTOBC&)
  {
    MayDay::Error(" PerEdgeCompGridVTOBC& operator= (const PerEdgeCompGridVTOBC&) not implemented");
    return *this;
  }


  virtual void operator()(FArrayBox&           a_vel,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    for (int dir = 0; dir < SpaceDim; dir++)
      {
	(*m_loBC[dir])(a_vel, a_valid, a_domain, a_dx, a_homogeneous);
	(*m_hiBC[dir])(a_vel, a_valid, a_domain, a_dx, a_homogeneous);
      }
  }


};




#include "NamespaceFooter.H"

#endif
