 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
//===========================================================================
// UnstructuredIO.H
// Operations related to reading and writing data from/to valid regions
// of AMR Hierarchies
//===========================================================================
#ifndef _UNSTRUCTUREDIO_H_
#define _UNSTRUCTUREDIO_H_

#include "RealVect.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "netcdf.h"
#include <set>
#include "Transformation.H"
#include "NamespaceHeader.H"

class UnstructuredData
{

public:
  typedef std::set<std::pair<int, Box> >  LevelBoxSet ;
private:

  int m_nComp;
  Vector<RealVect> m_dx;
  RealVect m_x0;
  Real m_time;
  int m_nLevel;
  Vector<int> m_ratio;
  Vector<int> m_level;
  LevelBoxSet m_levelBoxSet;
  Vector<Vector<int> > m_iv; 
  Vector<Vector< Real > > m_field;
  Vector<Vector<Real > > m_x;
  Box m_crseDomain;
  bool m_isDefined;

public:

 

  UnstructuredData(int a_nComp,  const RealVect& a_crseDx,  const Box& a_crseDomain,
		   const Vector<int>& a_ratio, const Real& a_time, const RealVect& a_x0 = RealVect::Zero );
  
  void define(int a_nComp,  const RealVect& a_crseDx,  const Box& a_crseDomain,
	      const Vector<int>& a_ratio, const Real& a_time, const RealVect& a_x0 = RealVect::Zero );

  UnstructuredData():m_isDefined(false){}

  void resize(int a_size);


  void append(int a_lev, const Box& a_box, const IntVect& a_iv, 
	      const Vector<Real>& a_data);
  
  int nComp() const
  {
    return m_nComp;
  }

 size_t nCell() const
  {
    return m_level.size();
  }

  size_t nLevel() const
  {
    return m_nLevel;
  }

  Box domain(int a_lev) const
  {
    if (a_lev == 0)
      return m_crseDomain;
    else
      return domain(a_lev - 1).refine(m_ratio[a_lev-1]);
  }

  Real time() const
  {
    return m_time;
  }

  Real& time()
  {
    return m_time;
  }


  const Vector<int>& ratio() const
  {
    return m_ratio;
  }
 
  const Vector<RealVect>& dx() const
  {
    return m_dx;
  }


  const Vector<int>& level() const
  {
    return m_level;
  }

  Vector<int>& level()
  {
    return m_level;
  }
  
  const LevelBoxSet& levelBoxSet() const
  {
    return m_levelBoxSet;
  }
  
  LevelBoxSet& levelBoxSet() 
  {
    return m_levelBoxSet;
  }

  const Vector<int>& iv(int a_dir) const
  {
    return m_iv[a_dir];
  }

  Vector<int>& iv(int a_dir) 
  {
    return m_iv[a_dir];
  }

  const Vector<Real>& field(int a_comp) const
  {
    return m_field[a_comp];
  }

  Vector<Real>& field(int a_comp)
  {
    return m_field[a_comp];
  }

  const Vector<Real>& x(int a_dir) const
  {
    return m_x[a_dir];
  }

  ///fill  a_nodeCoord with the a_dir co-ordinate of node a_node. 
  /**
     \param Vector<Real>& a_nodeCoord will be resized to this.size()
     \param a_node[dir] = 0 implies lo side, 1 implies high side
  */
  void computeNodeCoord(Vector<Real>& a_nodeCoord, int a_dir, const IntVect& a_node) const;
    
};

namespace UnstructuredIO 
{

  /// populate a UnstructuredData from a  Vector<LevelData<FArrayBox>*>
  void BStoUS( UnstructuredData& a_usData , 
			  const Vector<LevelData<FArrayBox>*>& a_bsData, 
			  const Vector<int>& a_ratio, 
			  bool a_validOnly);
  
  /// build block structed data from valid data
  void USToBS ( Vector<LevelData<FArrayBox>*>& a_bsdata, 
		   const UnstructuredData& a_usData);



  int defineCFDimension(int a_ncID, size_t a_len, const std::string& a_name);

  int defineCFVar(int a_ncID, int a_nDim, int* a_dimID, nc_type a_type,
		   const std::string& a_name,
		   const std::string& a_unit,
		   const std::string& a_stdName, 
		   const std::string& a_longName);

  void readCFVar(int a_ncID, const std::string& a_name,Vector<int>& a_data);

  void readCFVar(int a_ncID, int a_varID ,Vector<int>& a_data);

  void readCFVar(int a_ncID, const std::string& a_name,Vector<Real>& a_data);

  void readCFVar(int a_ncID, int a_varID ,Vector<Real>& a_data);

  void writeCFVar(int a_ncID, const std::string& a_name, const Vector<int>& a_data);

  void writeCFVar(int a_ncID, const std::string& a_name, const Vector<Real>& a_data);

  ///write UnstructuredData to a NetCDF-CF file
  void writeCF ( const std::string&,  const UnstructuredData& , 
		 const Vector<std::string>&, const std::string&,
		 const Transformation&);

  ///read UnstructuredData from a NetCDF-CF file
  void readCF ( UnstructuredData&, Vector<std::string>&, 
		 const std::string&);



}

#include "NamespaceFooter.H"
#endif
