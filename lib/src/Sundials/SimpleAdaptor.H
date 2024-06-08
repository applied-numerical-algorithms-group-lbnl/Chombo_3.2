#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SIMPLEADAPTOR_H
#define _SIMPLEADAPTOR_H

#include "ChomboSundialsAdaptor.H"

#include "LevelData.H"
#include "FArrayBox.H"

/// Class to interface a LevelData<FArrayBox> with SUNDIALS
/**
  This class implements all the required NVector methods for
  SUNDIALS to call on a LevelData<FArrayBox>, works in parallel.
  Output available using HDF5.

  Calling context must allocate memory in the LevelData, and then
  set the public Data member variable and others. This will call
  *delete* on Data in the destructor.

  Note that many of the SUNDIALS NVector operations allow the object 
  itself to be passed for an in-place update, check each function.
 */
class SimpleAdaptor : public ChomboSundialsAdaptor 
{

public:
  ///
  /**
    Constructor
   */
  SimpleAdaptor() : ChomboSundialsAdaptor() {};
  ///
  /**
    Destructor - calls delete on Data pointer.
   */
  virtual ~SimpleAdaptor()
  {
    if (m_ownData && (m_data != NULL))
      delete m_data;
  };

  /// 
  /**
     Define from an existing %LevelData<FArrayBox>.
     <ul>
       <li> lvlData: pointers to the %LevelData<FArrayBox>
       <li> ownData: (true) delete pointer in destructor
     </ul>
  */
  void define(LevelData<FArrayBox>* lvlData, bool ownData);

  /// 
  /**
     Get the %LevelData<FArrayBox> data pointer.
  */
  LevelData<FArrayBox>* getData() { return this->m_data; };

  /// 
  /**
     Get the flag for owning data deallocation
  */
  bool getOwnData() { return this->m_ownData; };

  ///
  /**
    SUNDIALS operation - calculates a*x+b*y and puts in this instance's data
   */
  virtual void linearSum(ChomboSundialsAdaptor& x, ChomboSundialsAdaptor& y,
      Real a, Real b);
  ///
  /**
    SUNDIALS operation - sets this instance's data to c
   */
  virtual void setConst(Real c);
  ///
  /**
    SUNDIALS operation - sets this instance's data to x*y
   */
  virtual void prod(ChomboSundialsAdaptor& ax, ChomboSundialsAdaptor& ay);
  ///
  /**
    SUNDIALS operation - sets this instance's data to x/y
   */
  virtual void div(ChomboSundialsAdaptor& ax, ChomboSundialsAdaptor& ay);
  ///
  /**
    SUNDIALS operation - sets this instance's data to c*x
   */
  virtual void scale(ChomboSundialsAdaptor& ax, Real c);
  ///
  /**
    SUNDIALS operation - sets this instance's data to abs(x)
   */
  virtual void abs(ChomboSundialsAdaptor& ax);
  ///
  /**
    SUNDIALS operation - sets this instance's data to 1/x
   */
  virtual void inv(ChomboSundialsAdaptor& ax);
  ///
  /**
    SUNDIALS operation - sets this instance's data to x+b
   */
  virtual void addConst(ChomboSundialsAdaptor& ax, Real b);
  ///
  /**
    SUNDIALS operation - returns MPI-sum dot product sum(data^T * input)
    across all MPI ranks 
   */
  virtual Real dotProd(ChomboSundialsAdaptor& ax);
  ///
  /**
    SUNDIALS operation - returns MPI-reduce across all MPI ranks' data max norm
   */
  virtual Real maxNorm();
  ///
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data 
    weighted RMS (root mean square) norm
   */
  virtual Real wRMSNorm(ChomboSundialsAdaptor& aw);
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual Real wRMSNormMask(ChomboSundialsAdaptor& aw, ChomboSundialsAdaptor& aid);
  ///
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data min
   */
  virtual Real min();
  ///
  /**
    SUNDIALS operation - returns MPI-reduce L1 norm (sum abs values)
   */
  virtual Real l1Norm();
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual Real wL2Norm(ChomboSundialsAdaptor& aw);
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual void compare(ChomboSundialsAdaptor& ax, Real b);
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual bool invTest(ChomboSundialsAdaptor& ax);
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual bool constrMask(ChomboSundialsAdaptor& ac, ChomboSundialsAdaptor& am);
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual Real minQuotient(ChomboSundialsAdaptor& adenom);

  ///
  /**
    SUNDIALS operation - returns the total number of interior cells
    multiplied by the number of components
   */
  virtual unsigned long getLength();
  ///
  /**
    Call exchange on the underlying LevelData<FArrayBox>
   */
  virtual void exchange()
  {
    m_data->exchange(m_data->interval());
  }
  ///
  /**
    Calls dumpLDFPar - careful, very verbose
   */
  virtual void print();
  ///
  /**
    SUNDIALS operation - not implemented, Aborts
   */
  virtual void printFile(FILE* outfile);
#ifdef CH_USE_HDF5
  ///
  /**
    Calls writeLevelname to produce HDF5 output
   */
  virtual void printFileHDF(const char* filename);
#endif

  ///
  /**
    Uses MPI to copy values from this to arg, only if DBLs are the same.
    Returns true if the copy was succesful, false otherwise.
   */
  virtual bool copyTo(ChomboSundialsAdaptor& a)
  {
    SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&a);
    if (sa == NULL) return false;
    if (!(m_data->disjointBoxLayout() == sa->m_data->disjointBoxLayout()))
      return false;
    m_data->copyTo(*sa->m_data);
    return true;
  }

  ///
  /**
   Factory, creates a clone of this and allocates new data (not copied).
   */
  virtual ChomboSundialsAdaptor* newAdaptor()
  {
    SimpleAdaptor* adaptor = new SimpleAdaptor();
    adaptor->m_dp = m_dp;
    adaptor->m_nComp = m_nComp;
    adaptor->m_ghost = m_ghost;
    adaptor->m_data = new LevelData<FArrayBox>;
    adaptor->m_data->define(m_dp, m_nComp, m_ghost);
    adaptor->m_ownData = true;
    return adaptor;
  }

private:
  // Private  Member variables
  //! Pointer to LevelData<FArrayBox> data this wraps
  LevelData<FArrayBox>* m_data = NULL;
  //! DisjointBoxLayout of the %LevelData<FArrayBox>
  DisjointBoxLayout m_dp;
  //! Number of components in the %LevelData<FArrayBox>
  int m_nComp;
  //! Number of ghost cells in the %LevelData<FArrayBox>
  IntVect m_ghost;
  //! If true, will call delete in destructor
  bool m_ownData = false;

};

#endif
