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
// FieldNames.cpp
// Provides a lookup tables between the various field names used in BISICLES
// with alternate names, units. The tables themselves (for now) are stored
// in ParmParse tables.
//===========================================================================
#include "FieldNames.H"
#include "NamespaceHeader.H"

/// create CFRecords for each element of \param a_names, append to \param a_records
void FieldNames::CFLookup(Vector<CFRecord>& a_records , const Vector<std::string> a_names)
{
  //a_records.resize(0);
  for (int i = 0; i < a_names.size(); i++)
    {
      a_records.push_back(CFRecord(a_names[i]));
    }

}

///remove awkward characters from \param a_name
void FieldNames::sanitize(std::string& a_name)
{
  const std::string awkward("/");
  const std::string repl("_");

  for (int i =0; i < a_name.size(); i++)
    {
      for (int j = 0; j < awkward.size(); j++ )
	{
	  if (a_name[i] == awkward[j])
	    {
	      a_name[i] = repl[0];
	    }
	}
    }
}

/// contruct a record by looking up \param a_name
FieldNames::CFRecord::CFRecord(const std::string& a_name)
{
  ParmParse pp(a_name.c_str());
  
  m_name = a_name;
  pp.query("name",m_name);
  sanitize(m_name);
    
  m_scale = 1.0;
  pp.query("scale",m_scale);

  ParmParse ppn(m_name.c_str());

  m_stdname = "";
  ppn.query("standard_name",m_stdname);

  m_longname = a_name;
  ppn.query("long_name",m_longname);

  m_units="";
  ppn.query("units",m_units);
}


#include "NamespaceFooter.H"
