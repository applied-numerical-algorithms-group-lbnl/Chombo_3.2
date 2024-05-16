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
// ncio.H
// read/write FABs from netcdf files
//===========================================================================
#ifndef _FABNCIO_H_
#define _FABNCIO_H_
#include "FArrayBox.H"
#include "RealVect.H"
#include <map>
#include "CH_HDF5.H"
#include "DomainDiagnosticData.H"
#include "NamespaceHeader.H"

#ifdef HAVE_NETCDF
namespace NCIO
{


  /// write a single FAB to a (CF Compliant) netcdf file
  void writeFAB(const std::string& a_file,
		const Vector<std::string>& a_names,
		const Vector<std::string>& a_cf_standard_names,
		const Vector<std::string>& a_cf_units,
		const Vector<std::string>& a_cf_long_names,
		const FArrayBox& a_fab, 
		const Real& a_dx,
		const RealVect& a_x0, 
		int a_epsg,
		const DomainDiagnosticData& a_dd,
		const std::string& a_flattenInfo,
		const HDF5HeaderData& a_file_header);

  /// read a single FAB from a netcdf file
  void readFAB(const std::string& a_file,
	       const Vector<std::string>& a_var,
	       FArrayBox& a_fab,
	       Real& a_dx);


  
}
#endif
#include "NamespaceFooter.H"
#endif
