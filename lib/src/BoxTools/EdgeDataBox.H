#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EDGEDATABOX_H_
#define _EDGEDATABOX_H_

#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

/// A FArrayBox-like container for edge-centered fluxes
class EdgeDataBox
/** This is a class to contain edge-centered fluxes on a box.
 */
{

public:
  /// default constructor
  EdgeDataBox();

  /// constructs EdgeDataBox on cell-centered box with n components
  EdgeDataBox(const Box& bx, int n=1);

  /// destructor
  ~EdgeDataBox();

  /// resize EdgeDataBox similar to BaseFab::resize()
  //void resize(const Box& bx, int n=1);

  ///
  void define(const Box& bx, int n=1);

  /// this function returns the EdgeDataBox to the undefined state
  void clear();

  /// {\bf access functions}

  /// number of components
  int nComp() const;

  /// returns cell-centered box which defines EdgeDataBox
  const Box& box() const;

  /// returns edge-centered data in direction dir
  FArrayBox& getData(const int dir);

  ///
  const FArrayBox& getData(const int dir) const;

  /// returns FArrayBox in direction dir
  FArrayBox& operator[] (const int dir);

  /// constant version
  const FArrayBox& operator[] (const int dir) const;

  /// {\bf data modification functions}
  /// set all fluxes to val
  void setVal(const Real val);

  /// set fluxes in direction dir to val
  void setVal(const Real val, const int dir);

  /// more specific setVal
  void setVal(const Real val, const int dir, const int startComp,
              const int nComp);

  /// sets fluxes on edges surrounding cell-centered box bx
  void setVal(const Real val, const Box& bx);

  /*@ManMemo: most specific setVal -- sets fluxes on edges surrounding
    cell-centered box bx */
  void setVal(const Real val, const Box& bx, const int dir,
              const int startComp, const int nComp);

  /// copy from src to this EdgeDataBox -- sizes must be identical
  void copy(const EdgeDataBox& src);

  /// copy on overlap, for all directions
  void copy(const EdgeDataBox& src, const int srcComp,
            const int destComp, const int numComp);

  /// copy on overlap of EdgeDataBoxes, in direction dir
  void copy(const EdgeDataBox& src, const int dir, const int srcComp,
            const int destComp, const int numComp);

  /*@ManMemo: copies from a subsection of one box into another.
    Assumes the boxes are both in the same index space, and that
    box R is completely contained in both the src and destination
    boxes.  This function equired by BoxLayoutData */
  void copy(const Box& RegionFrom,
            const Interval& Cdest,
            const Box& RegionTo,
            const EdgeDataBox& src,
            const Interval& Csrc);

  /// Modifes this EdgeDataBox by adding src in the CELL-CENTERED sub-box
  /** Modifies this EdgeDataBox by pointwise addition of values in the argument
      EdgeDataBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
      to this EdgeDataBox's components (a_destcomp : a_destcomp+numcomp-1)
      where the domain of this EdgeDataBox intersects the a_subbox.  The actual
      directionally-dependent subbox is given by 
      a_subbox.surroundingNodes();
      a_subbox.enclosedCells(dir);
      NOTE:
      a_subbox must be contained in the cell-centered Box of this EdgeDataBox.
      Returns *this
   */
  EdgeDataBox& plus(const EdgeDataBox& a_src,
                    const Box&     a_subbox,
                    int            a_srccomp,
                    int            a_destcomp,
                    int            a_numcomp = 1);

  //  plus operator with calling sequence similar to matching copy (for use in additive exchange operations)
  void plus (const Box&        srcbox,
             const Interval&   destcomps,
             const Box&        destbox,
             const EdgeDataBox&    src,
             const Interval&   srccomps);


  
  
  /// {\bf Linearization Functions}

  /*@ManMemo: returns size, in number of bytes, of a flat linear
    representation of data in components comps in edges around
    cell-centered box R */
  size_t size(const Box& bx, const Interval& comps) const;

  /*@ManMemo:Write a linear representaion of the internal data for
    the edges surrounding cell-centered box R.  Assumes that sufficient
    memory for the buffer has already been allocated by the caller */
  void linearOut(void*buf, const Box& R, const Interval& comps) const;

  /// like linearOut, but returns current position in the buffer
  void* linearOut2(void*buf, const Box& R, const Interval& comps) const;
  
  /// 
  void linearIn(void* buf);

  /// read a linear representation of the data over the box \rm R
  /** 
      Reads in the output of LinearOut
  */
  void linearIn(void*buf, const Box& R, const Interval& comps);

  /// like linearIn, but returns current position in the buffer
  void* linearIn2(void*buf, const Box& R, const Interval& comps);

  ///
  static int preAllocatable()
  {
    return 0;
  }


protected:
  ///
  Box m_bx;

  ///
  int m_nvar;

  /// BL_SPACEDIM FArrayBoxes which hold fluxes
  Vector<FArrayBox*> m_data;


private:
  /// these are disallowed

  EdgeDataBox (const EdgeDataBox&);
  EdgeDataBox& operator = (const EdgeDataBox&);

};

#include "NamespaceFooter.H"

#endif
