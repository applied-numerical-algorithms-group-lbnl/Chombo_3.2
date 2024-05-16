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
// FieldNames.H
// Provides a lookup tables between the various field names used in BISICLES
// with alternate names, units. The tables themselves (for now) are stored
// in ParmParse tables.
//===========================================================================
#ifndef _FIELDNAMES_H_
#define _FIELDNAMES_H_
#include "ParmParse.H"
#include "NamespaceHeader.H"

namespace FieldNames
{
  class CFRecord
  {
    std::string m_name; // name of field/variable in e.g netcdf files
    std::string m_stdname; // standard name (CF)
    std::string m_longname; // long name (used by CF when there is no standard name)
    std::string m_units; // 
    Real m_scale; //conversion factor 
  public:
    /// contruct a record by looking up \param a_name
    CFRecord(const std::string& a_name);
    
    const std::string& name() const {return m_name;}
    const std::string& standardName() const {return m_stdname;}
    const std::string& longName() const {return m_longname;}
    const std::string& units() const {return m_units;}
    Real scale() const {return m_scale;}
  private:
    CFRecord(){};
    
  };
  /// create CFRecords for each element of \param a_names, append to \param a_records
  void CFLookup(Vector<CFRecord>& a_records , const Vector<std::string> a_names);
  
  ///remove awkward characters from \param a_name
  void sanitize(std::string& a_name);
}

#include "NamespaceFooter.H"
#endif
