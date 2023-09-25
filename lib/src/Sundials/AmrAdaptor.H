#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRADAPTOR_H
#define _AMRADAPTOR_H

#include "ChomboSundialsAdaptor.H"

#include "LevelData.H"
#include "FArrayBox.H"

/// Data holder for %AMR hierarchy data to interact with SUNDIALS.
/**
  This class supports most of the operations needed by SUNDIALS nvector,
  mapped onto an %AMR hierarchy of LevelData<FArrayBox>. It works in parallel
  and serial, and uses masks to make sure covered cells on coarse levels don't 
  contribute to calculations needed by SUNDIALS. It is dumb, in that if any of
  the %AMR hierarchy is regridded or redefined, this must be redefined as well.

  Note that many of the SUNDIALS NVector operations allow the object itself to 
  be passed for an in-place update, check each function.
 */
class AmrAdaptor : public ChomboSundialsAdaptor 
{

public:
  /// 
  /**
     Constructor.
  */
  AmrAdaptor() : ChomboSundialsAdaptor() {};

  /// 
  /**
     Destructor.
  */
  virtual ~AmrAdaptor() {
    int nlvl = m_amrData.size();
    for (int lvl=0; lvl < nlvl; lvl++)
    {
      if (m_ownData && (m_amrData[lvl] != NULL))
        delete m_amrData[lvl];
      if (m_amrMask[lvl] != NULL)
        delete m_amrMask[lvl];
    }
  };

  /// 
  /**
     Factory.
  */
  virtual ChomboSundialsAdaptor* newAdaptor();

  /// 
  /**
     Define from an existing %AMR hierarchy.
     <ul>
       <li> amrData: vector of pointers to all the %LevelData<FArrayBox> in the hierarchy
       <li> refRatio: vector of integer refinement ratios between levels
       <li> ownData: (true) delete pointers in destructor
     </ul>
     Redefines and calculates masks. Calls delete on amrData and masks pointers in destructor.
  */
  void define(Vector<LevelData<FArrayBox>*> amrData, Vector<int> refRatio,
              bool ownData);

  /// 
  /**
     Get the %AMR data pointers.
  */
  Vector<LevelData<FArrayBox>*> getData() { return this->m_amrData; };

  /// 
  /**
     Get the refinement ratios.
  */
  Vector<int> getRefRatio() { return this->m_refRatio; };

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
    SUNDIALS operation - returns MPI-sum dot product across all ranks, 
    sum(data^T * input), masked for covered cells
  */
  virtual Real dotProd(ChomboSundialsAdaptor& ax);
  /// 
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data max norm,
    masked for covered cells
  */
  virtual Real maxNorm();
  /// 
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data 
    weighted RMS (root mean square) norm, masked for covered cells
  */
  virtual Real wRMSNorm(ChomboSundialsAdaptor& aw);
  /// 
  /**
    SUNDIALS operation - not implemented, Aborts
  */
  virtual Real wRMSNormMask(ChomboSundialsAdaptor& aw, ChomboSundialsAdaptor& aid);
  /// 
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data min,
    masked for covered cells
  */
  virtual Real min();
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
    SUNDIALS operation - returns the total number of *uncovered* %AMR cells
    multiplied by the number of components
  */
  virtual unsigned long getLength();
  /// 
  /**
    SUNDIALS operation - returns MPI-reduce across all ranks' data L1 norm
    (sum abs values), no dx or refinement ratio, masked for covered cells
  */
  virtual Real l1Norm();
  /// 
  /**
    Not implemented, Aborts
  */
  virtual void print();
  /// 
  /**
    Not implemented, Aborts
  */
  virtual void printFile(FILE* outfile);
#if CH_USE_HDF5
  /// 
  /**
    Not implemented, Aborts
  */
  virtual void printFileHDF(const char* filename);
#endif
  /// 
  /**
    Not implemented, Aborts
  */
  virtual bool copyTo(ChomboSundialsAdaptor& a);
  /// 
  /**
    Not implemented, Aborts
  */
  virtual void exchange();


  // Member variables
  //! Pointers to %AMR hierarchy data
  Vector<LevelData<FArrayBox>*> m_amrData;
  //! %AMR hierarchy refinement ratios between levels
  Vector<int> m_refRatio;
  //! %AMR mask, 1 for uncovered values, 0 for covered values
  Vector<LevelData<FArrayBox>*> m_amrMask; // To help with efficient ops
  //! Cached nvector length, sum of all unmasked values, across all ranks
  unsigned long m_length = -1;
  //! Number of components in the %LevelData<FArrayBox>'s
  int m_comp = -1;
  //! If true, will call delete in destructor
  bool m_ownData = false;

private:
};

#endif
