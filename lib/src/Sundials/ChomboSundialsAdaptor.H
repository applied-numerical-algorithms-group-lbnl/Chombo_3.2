#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHOMBOSUNDIALSADAPTOR_H
#define _CHOMBOSUNDIALSADAPTOR_H

#include <stdio.h>
#include "REAL.H"

/// Virtual base class contract to interact with SUNDIALS
/**
  Declares virtual functions needed to satisfy the SUNDIALS NVector
  functions. Derived classes need to implement these, with care taken
  to "linearize" Chombo data, e.g. avoid counting duplicate values,
  include all components, or ghost cells, etc. in NVector operations.
  Note that many operations can also take "this" object as an argument;
  in that case they should do an in-place update of the object's data.
 */
class ChomboSundialsAdaptor
{
public:
  ///
  /**
     Constructor.
  */
  ChomboSundialsAdaptor() {};

  ///
  /**
     Destructor.
  */
  virtual ~ChomboSundialsAdaptor() {};

  /// 
  /**
    SUNDIALS operation - calculates a*x+b*y and puts in this instance's data
  */
  virtual void linearSum(ChomboSundialsAdaptor& x, ChomboSundialsAdaptor& y,
      Real a, Real b) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to c
  */
  virtual void setConst(Real c) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to x*y
  */
  virtual void prod(ChomboSundialsAdaptor& ax, ChomboSundialsAdaptor& ay) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to x/y
  */
  virtual void div(ChomboSundialsAdaptor& ax, ChomboSundialsAdaptor& ay) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to c*x
  */
  virtual void scale(ChomboSundialsAdaptor& ax, Real c) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to abs(x)
  */
  virtual void abs(ChomboSundialsAdaptor& ax) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to 1/x
  */
  virtual void inv(ChomboSundialsAdaptor& ax) = 0;
  /// 
  /**
    SUNDIALS operation - sets this instance's data to x+b
  */
  virtual void addConst(ChomboSundialsAdaptor& ax, Real b) = 0;
  /// 
  /**
    SUNDIALS operation - returns dot product sum(data^T * x)
  */
  virtual Real dotProd(ChomboSundialsAdaptor& ax) = 0;
  /// 
  /**
    SUNDIALS operation - returns all data max norm
  */
  virtual Real maxNorm() = 0;
  /// 
  /**
    SUNDIALS operation - returns all data weighted RMS (root mean square) norm
  */
  virtual Real wRMSNorm(ChomboSundialsAdaptor& aw) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual Real wRMSNormMask(ChomboSundialsAdaptor& aw, ChomboSundialsAdaptor& aid) = 0;
  /// 
  /**
    SUNDIALS operation - returns all data min
  */
  virtual Real min() = 0;
  ///
  /**
    SUNDIALS operation - returns all data L1 norm
    (sum abs values)
  */
  virtual Real l1Norm() = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual Real wL2Norm(ChomboSundialsAdaptor& aw) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual void compare(ChomboSundialsAdaptor& ax, Real b) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual bool invTest(ChomboSundialsAdaptor& ax) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual bool constrMask(ChomboSundialsAdaptor& ac, ChomboSundialsAdaptor& am) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for completeness
  */
  virtual Real minQuotient(ChomboSundialsAdaptor& adenom) = 0;

  /// 
  /**
    SUNDIALS operation - returns the total number of cells in the nvector
    NOTE: careful for covered %AMR cells, ghost cells, etc. And should
    probably be multiplied by the number of components
  */
  virtual unsigned long getLength() = 0;
  /// 
  /**
    Placeholder for MPI. Does not have to be implemented, here for futures
  */
  virtual void exchange() = 0;
  /// 
  /**
    Placeholder for MPI. Does not have to be implemented, but here for futures
  */
  virtual bool copyTo(ChomboSundialsAdaptor& a) = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for futures
  */
  virtual void print() = 0;
  /// 
  /**
    SUNDIALS - does not have to be implemented, but here for futures
  */
  virtual void printFile(FILE* outfile) = 0;
  /// 
  /**
    Factory
  */
  virtual ChomboSundialsAdaptor* newAdaptor() = 0;

private:
};

#endif
