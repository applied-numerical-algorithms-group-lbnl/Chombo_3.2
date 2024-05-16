#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRMelange.H"
#include "AmrIce.H"
#include "FineInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "DivergenceF_F.H"
#include "CellToEdge.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "LevelMappedDerivatives.H"
#include "VCAMRPoissonOp2.H"
#include "SigmaCSF_F.H"
#include "amrIceF_F.H"
#include "BisiclesF_F.H"
#include "CH_HDF5.H"
#include "IceConstants.H"
#include "BiCGStabSolver.H"
#include "NamespaceHeader.H"

MelangeIceObserver::MelangeIceObserver()
  : m_melangePtr( new AMRMelange() ),  m_next_increment_positive(true)
{
}



MelangeIceObserver::~MelangeIceObserver()
{
  if (m_melangePtr != NULL)
    {
      delete m_melangePtr;
      m_melangePtr = NULL;
    }
}

AMRMelange& MelangeIceObserver::melange() const
{
  return *m_melangePtr;
}





void MelangeIceObserver::notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce)
{

  pout() <<  "MelangeIceObserver::notify(" << a_n << ")" << std::endl;

  if (a_n == AmrIce::Observer::PreVelocitySolve)
    {
      m_melangePtr->define(a_amrIce.grids(), a_amrIce.refRatios(),  a_amrIce.finestLevel(), a_amrIce.dx(0));
    }
  else if (a_n == AmrIce::Observer::PostGeometryUpdate)
    {
      
      m_melangePtr->timestep(a_amrIce.dt(), a_amrIce); 
    }
  else if (a_n == AmrIce::Observer::PreCalving)
    {
      CH_assert(MelangeIceObserver::m_next_increment_positive);
      m_next_increment_positive = false;
      m_melangePtr->increment(a_amrIce, 1.0);
    }
  else if (a_n == AmrIce::Observer::PostCalving)
    {
      CH_assert(!MelangeIceObserver::m_next_increment_positive);
      m_next_increment_positive = true;
      m_melangePtr->increment(a_amrIce, -1.0);
    }


}
#ifdef CH_USE_HDF5
void MelangeIceObserver::addPlotVars(Vector<std::string>& a_vars)
{
  m_melangePtr->addPlotVars(a_vars);
}

void MelangeIceObserver::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  m_melangePtr->writePlotData(a_data, a_level);
}

/// fill a_var with the names of variables to add to the checkpoint file 
void MelangeIceObserver::addCheckVars(Vector<std::string>& a_vars)
{
  m_melangePtr->addCheckVars(a_vars);
}
  
/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void MelangeIceObserver::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  m_melangePtr->writeCheckData(a_handle, a_level);
}
  
/// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
void MelangeIceObserver::readCheckData(HDF5Handle& a_handle, HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  m_melangePtr->readCheckData(a_handle, a_header, a_level, a_grids);
}
#endif


AMRMelange::~AMRMelange()
{

  if (m_external_source != NULL)
    {
      delete m_external_source;
    }
  
  for (int lev = 0; lev < m_melange.size(); lev++)
    {
      if (m_melange[lev] != NULL)
	{
	  delete m_melange[lev]; m_melange[lev] = NULL;
	}
    }


 
}
AMRMelange::AMRMelange()
{
  m_time_step = 0;
  m_time = 0.0;

  m_external_source = SurfaceFlux::parse("melange_model.external_source");
  if (m_external_source == NULL)
    {
      m_external_source = new zeroFlux();
    }


  ParmParse pp( "melange_model" );

  
  m_diffusion_factor = 0.0;
  pp.query("diffusion_factor", m_diffusion_factor);
  
}

const LevelData<FArrayBox>* AMRMelange::melange(int a_level) const
{
  if (!(m_melange.size() > a_level))
    {
      std::string msg("AMRMelange::melange !(m_melange.size() > a_level)");
      pout() << msg <<endl;
      CH_assert((m_melange.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_melange[a_level];
  if (ptr == NULL)
    {
      std::string msg("AMRMelange::melange m_melange[a_level] == NULL");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }
  return ptr;
}

void AMRMelange::define
(const Vector<DisjointBoxLayout>& a_grids, 
 const Vector<int>& a_ratio,
 int a_finestLevel,
 const RealVect& a_crseDx)
{
  CH_TIME("AMRMelange::define");
  pout() <<  "AMRMelange::define" << std::endl;
  
  if (m_dx.size() > 0)
    {
      if ( !(m_dx[0] == a_crseDx) )
	{
	  std::string msg("AMRMelange::define, incompatible mesh");
	  pout() << msg << std::endl;
	  CH_assert(m_dx[0] == a_crseDx);
	  MayDay::Error(msg.c_str());
	}
    }

  //update the mesh hierarchy
  m_finestLevel = a_finestLevel;
  m_grids.resize(m_finestLevel  + 1);
  m_ratio.resize(m_finestLevel + 1, 1);
  for (int lev = 0; lev <= m_finestLevel ; lev++)
    {
      m_grids[lev] = a_grids[lev];
      if (lev < m_finestLevel )
	{
	  m_ratio[lev] = a_ratio[lev];
	}
    }
  
  m_dx.resize(m_finestLevel + 1);
  m_dx[0] = a_crseDx;
  for (int lev = 1; lev <= m_finestLevel;  lev++)
    {
      m_dx[lev] = m_dx[lev-1] / Real(m_ratio[lev-1]);
    }


  //copy any previous melange data pointers
  Vector<LevelData<FArrayBox>* > prevMelange;
  prevMelange.resize(m_melange.size());
  for (int lev = 0; lev < m_melange.size(); lev++)
    {
      prevMelange[lev] = m_melange[lev];
    }

  //initialize the melange data
  m_melange.resize(m_finestLevel + 1);
  for (int lev = 0; lev <= m_finestLevel; lev++)
    {
      m_melange[lev] = new LevelData<FArrayBox>(m_grids[lev],  MELANGE_N_COMP , MELANGE_N_GHOST * IntVect::Unit);
    }
  if (prevMelange.size() == 0)
    {
     
      // for now, set to zero, but need to set initial conditions from input data of some sort, 
      // or read checkpoints
      for (int lev = 0; lev <= m_finestLevel; lev++)
       	{
       	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       	    {
       	      FArrayBox& melange = (*m_melange[lev])[dit];
	      melange.setVal(0.0);
	      
	    }
	}
    }

  if (prevMelange.size() > 0)
    {
      // if previous melange data exists, interpolate onto the new mesh
      for (int lev = 0; lev <= m_finestLevel; lev++)
	{

	  //FIXME : only doing this to avoid stupid value outside the domain, should be a BC thing....
	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       	    {
       	      FArrayBox& melange = (*m_melange[lev])[dit];
	      melange.setVal(0.0);  
	    }


	  if (lev > 0)
	    {
	      FineInterp fi(m_grids[lev], m_melange[lev]->nComp(), 
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi.interpToFine(*m_melange[lev], *m_melange[lev-1]);

	      //need to fill the ghost cells on CF-interfaces now, because *m_melange[lev] has to be
	      //in good shape at all times if the interface with MelangeConstitutiveRelation is to work
	      
	      PiecewiseLinearFillPatch ghostFiller(m_grids[lev],m_grids[lev-1],m_melange[lev]->nComp(),
						   m_grids[lev-1].physDomain(), m_ratio[lev-1], m_melange[lev]->ghostVect()[0]);
	      ghostFiller.fillInterp(*m_melange[lev], *m_melange[lev-1],*m_melange[lev-1],1.0,0,0,m_melange[lev]->nComp());

	    }
	  Interval ival(0,m_melange[lev]->nComp()-1);

	  if (prevMelange.size() > lev && prevMelange[lev])
	    prevMelange[lev]->copyTo(ival,  *m_melange[lev], ival);
	  m_melange[lev]->exchange();

	}
      //free old data
      for (int lev =0; lev < prevMelange.size(); lev++)
	{
	  if (prevMelange[lev] != NULL)
	    {
	      delete prevMelange[lev]; prevMelange[lev] = NULL;
	    }	
	}

    }

}

void AMRMelange::increment(AmrIce& a_amrIce, Real a_scale)
{
  pout() <<  "AMRMelange::increment ice thickness * " << a_scale << std::endl;
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      const LevelData<FArrayBox>& iceThickness = a_amrIce.geometry(lev)->getH();
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	{
	  (*m_melange[lev])[dit].plus( iceThickness[dit], a_scale);
	}
    }
  
}


void AMRMelange::timestep(Real a_dt, AmrIce& a_amrIce)
{
  CH_TIME("AMRMelange::timestep");
  pout() <<  "AMRMelange::timestep" << std::endl;
  // fill ghost data
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (lev > 0)
	{
	  int nGhost = m_melange[lev]->ghostVect()[0];
	  PiecewiseLinearFillPatch pwl(m_grids[lev],  m_grids[lev-1], 1, 
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost);    
	  // since we're not subcycling, don't need to interpolate in time
	  Real time_interp_coeff = 0.0;
	  pwl.fillInterp(*m_melange[lev], *m_melange[lev-1], *m_melange[lev-1],
			 time_interp_coeff,0, 0, 1);
	}  
      m_melange[lev]->exchange();
    }
      
  Vector<LevelData<FArrayBox>* >  source;
  source.resize(m_finestLevel + 1);
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      source[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
    }
  // source terms
  computeSource(source, a_amrIce, a_dt);
  
  // if we had an advection term, might deal with it here. But for now, u_adv = 0.0

  
  // finally, update the melange thickness
  const Vector<RefCountedPtr<LevelSigmaCS> >& geometry = a_amrIce.amrGeometry();
  updateMelange(m_melange, source, geometry,  a_dt);


  //free temporary data.
  for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (source[lev] != NULL)
	{
	  delete source[lev]; source[lev] = NULL;
	}
    }

  m_time_step++;
  m_time += a_dt;

}

///Compute the source part of equation dM/dt + div( - k grad M ) = source.
/**
   The major source of melange will be calving from the ice sheet, which
   is accounted for through calls to AMRMelange::increment. The major sink
   will be user supplied, for example the UKESM ice-ocean coupler will
   provide a sink to mathc the ocean model's iceberg or freshwater sources

*/
void AMRMelange::computeSource(Vector<LevelData<FArrayBox>* >& a_source,
			       AmrIce& a_amrIce,  Real a_dt)
{
  
  for (int lev=0; lev <= m_finestLevel ; ++lev)
    {
      m_external_source->evaluate(*a_source[lev], a_amrIce, lev, a_dt);
    }
  
}



/// solve melange - div( k grad (melange) ) = melange_old + a_dt + source*a_dt;
/**
   k(x,y) is made up accoriding to 2 rules (1) no melange transport where surface > 0 (2) faster transport in deeper water

 */
void AMRMelange::updateMelange
(Vector<LevelData<FArrayBox>* >& a_melange,
 const Vector<LevelData<FArrayBox>* >& a_source,
 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
 Real a_dt)
{
  CH_TIME("AmrMelange::updateMelange");

  /// Natural BCs. 
  BCHolder bc(ConstDiriNeumBC(IntVect::Zero, RealVect::Zero,  IntVect::Zero, RealVect::Zero));

  //implicit Euler : solve (I - dt K grad ( M)  ) H = M_prev + dt * S. \todo - switch to Crank-Nicolson?
  
  Vector<RefCountedPtr<LevelData<FArrayBox> > > I(m_finestLevel + 1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > K(m_finestLevel + 1);
  //Vector<LevelData<FArrayBox>* > M(m_finestLevel + 1);
  Vector<LevelData<FArrayBox>* > rhs(m_finestLevel + 1);
  Vector<DisjointBoxLayout> grids(m_finestLevel + 1);
  
  for (int lev=0; lev <= m_finestLevel ; ++lev)
    {
      I[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids[lev], 1, IntVect::Unit));
      K[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Unit));
      rhs[lev] = new LevelData<FArrayBox>(m_grids[lev], 1, IntVect::Unit);
      //M[lev] = new LevelData<FArrayBox>(m_grids[lev], 1, IntVect::Unit);

      
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	{
	  (*I[lev])[dit].setVal(1.0);
	  ///rhs = M_old + dt * S
	  (*rhs[lev])[dit].copy((*a_melange[lev])[dit]);
	  (*rhs[lev])[dit].plus((*a_source[lev])[dit],a_dt);

	  // made up K = max(w,0) * max(1-s,0) 
	  //surface elevation
	  const FArrayBox& s = a_geometry[lev]->getSurfaceHeight()[dit];
	  // cavity thickness
	  FArrayBox w(s.box(),1); 
	  w.copy(s);
	  w -=  a_geometry[lev]->getH()[dit];
	  w -= a_geometry[lev]->getTopography()[dit];
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      Box b ( m_grids[lev][dit].surroundingNodes(dir) );
	      for (BoxIterator bit(b); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  IntVect ivp = iv + BASISV(dir);
		  Real km =  std::max( w(iv), 0.0) * std::max(1.0-s(iv) , 0.0);
		  Real kp =  std::max( w(ivp), 0.0) * std::max(1.0-s(ivp) , 0.0);
		  (*K[lev])[dit][dir](iv) = m_diffusion_factor * std::max(kp,km);
		}
	    }
	  
	}

      rhs[lev]->exchange();
      I[lev]->exchange();
     
    }

   VCAMRPoissonOp2Factory poissonOpFactory;//= new VCAMRPoissonOp2Factory;
   poissonOpFactory.define(m_grids[0].physDomain(), m_grids , m_ratio,
			   m_dx[0][0], bc, 1.0, I,  a_dt, K);

   //Plain MG
   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
   AMRMultiGrid<LevelData<FArrayBox> > mgSolver;
   mgSolver.define(m_grids[0].physDomain(), poissonOpFactory , &bottomSolver, m_finestLevel+1);
   //parse these
   mgSolver.m_eps = 1.0e-10;
   mgSolver.m_normThresh = 1.0e-10;
   
   int numMGSmooth = 4;
   mgSolver.m_pre = numMGSmooth;
   mgSolver.m_post = numMGSmooth;
   mgSolver.m_bottom = numMGSmooth;
   
   mgSolver.solve(a_melange, rhs, m_finestLevel , 0,  false);

  //dont allow melange to be negative
  for (int lev=0; lev<= m_finestLevel; lev++)
    { 
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
        {
	  FArrayBox& melange = (*m_melange[lev])[dit];
	  Real lim = 0.0;
          FORT_MAXFAB1(CHF_FRA( melange), 
                       CHF_CONST_REAL(lim), 
                       CHF_BOX(melange.box()));
        }
    }
  
}

#ifdef CH_USE_HDF5
void AMRMelange::addPlotVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("melangeThickness");
}

void AMRMelange::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  for (DataIterator dit(m_grids[a_level]);dit.ok();++dit)
    {
      a_data[dit].copy( (*m_melange[a_level])[dit],0,0,1);
    }
}



/// fill a_var with the names of variables to add to the checkpoint file 
void AMRMelange::addCheckVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("melangeThck");
}
  
/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void AMRMelange::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  write(a_handle, *m_melange[a_level], "melangeThck", m_melange[a_level]->ghostVect());
}
  
/// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
void AMRMelange::readCheckData(HDF5Handle& a_handle, HDF5HeaderData&  a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  bool containsMelangeData(false);
  map<std::string, std::string>::const_iterator i;
  for (i = a_header.m_string.begin(); i!= a_header.m_string.end(); ++i)
    {
      containsMelangeData |= (i->second == "melangeThckData");
    }

  if (containsMelangeData)
    {

      if (m_melange.size() <= a_level)
	{
	  m_melange.resize(a_level + 1, NULL);
	  if (m_melange[a_level] != NULL)
	    {
	      delete m_melange[a_level];
	    }
	  
	}
      m_melange[a_level] = new LevelData<FArrayBox>(a_grids, MELANGE_N_COMP, MELANGE_N_GHOST * IntVect::Unit); 
      int dataStatus =  read<FArrayBox>(a_handle, *m_melange[a_level], "MelangeData", a_grids);
      if (dataStatus != 0)
	{
	  MayDay::Error("failed to read melange data from checkpoint, but the header indicated its presence");
	}
    }
}
#endif


void MelangePhysIBC::define(const ProblemDomain& a_domain,
			  const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

PhysIBC* MelangePhysIBC::new_physIBC()
{
  MelangePhysIBC *ptr = new MelangePhysIBC();
  ptr->define(PhysIBC::m_domain, PhysIBC::m_dx);
  return static_cast<PhysIBC*>(ptr);
}

void MelangePhysIBC::initialize(LevelData<FArrayBox>& a_U)
{
  MayDay::Error("MelangePhysIBC::initialize not implemented");
}

void MelangePhysIBC::primBC(FArrayBox&            a_WGdnv,
			   const FArrayBox&      a_Wextrap,
			   const FArrayBox&      a_W,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const Real&           a_time)
{

  //first off, set WGdnv on the domain boundaries to a_Wextrap
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
	  for (BoxIterator bit(boundaryBox); bit.ok(); ++bit)
	    {
	      const IntVect& i = bit();
	      //All boundaries are outflows
	      a_WGdnv(i,0) = std::max(0.0,a_Wextrap(i,0));
	    }
	}
    }
 
}

void MelangePhysIBC::setBdrySlopes(FArrayBox&       a_dW,
				 const FArrayBox& a_W,
				 const int&       a_dir,
				 const Real&      a_time)
{
  //one-sided differences are fine, so do nothing
}

void MelangePhysIBC:: artViscBC(FArrayBox&       a_F,
		 const FArrayBox& a_U,
		 const FArrayBox& a_divVel,
		 const int&       a_dir,
		 const Real&      a_time)
{
  MayDay::Error("MelangePhysIBC:artViscBC not implemented");
}


#include "NamespaceFooter.H"
