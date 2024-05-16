#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _INTFINEINTERP_H_
#define _INTFINEINTERP_H_

#include "BaseFab.H"
#include "LevelData.H"
#include "NamespaceHeader.H"

/// replaces fine level BaseFAB<int> data with coarse level data
/***
    This class replaces data at a fine level of refinement with data
    given at a coarser level of refinement. In that sense, it is
    an integer analogue of FineInterp, except that it only supports
    piecewise constant interpolation.
 */
class IntFineInterp
{
public:
   ///
  /**
     Default constructor.  User must subsequently call define().
  */
  IntFineInterp();

  ///
  /**
     Destructor.
  */
  ~IntFineInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_ghost (not modified): ghost regions to fill.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

  */
  IntFineInterp(const DisjointBoxLayout& a_fine_domain,
		const int& a_numcomps,
		const int& a_ref_ratio,
		const IntVect& a_ghost,
		const Box& a_fine_problem_domain);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_ghost (not modified): ghost regions to fill.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\
  */
  IntFineInterp(const DisjointBoxLayout& a_fine_domain,
		const int& a_numcomps,
		const int& a_ref_ratio,
		const IntVect& a_ghost,
		const ProblemDomain& a_fine_problem_domain);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_ghost (not modified): ghost regions to fill.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
	 const IntVect& a_ghost,                  // the ghost region
         const Box& a_fine_problem_domain);      // problem domain

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_ghost (not modified): ghost regions to fill.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
	 const IntVect& a_ghost,
         const ProblemDomain& a_fine_problem_domain);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.

     {\bf This:}\\
     This object is not modified.
  */
  bool
  isDefined() const;

  void pwcInterpToFine(LevelData<BaseFab<int> >& a_fine_data,
		       const LevelData<BaseFab<int> >& a_coarse_data) ;

private:
  void pwcinterpGridData(BaseFab<int>& a_fine,
			 const BaseFab<int>& a_coarse,
			 const Box& a_coarsened_fine_box,
			 int a_ref_ratio) const;

protected:
  bool is_defined;
  // the refinement ratio
  int m_ref_ratio;
  // work array for the coarse level data in a domain that is the
  // outline of the fine level domain on the coarse level
  LevelData<BaseFab<int> > m_coarsened_fine_data;
  // coarse level problem domain
  ProblemDomain m_coarse_problem_domain;
  // ghost vector on the fine domain
  IntVect m_fineGhost;
};


#include "NamespaceFooter.H"
#endif
