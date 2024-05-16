#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "IntFineInterp.H"
#include "IntInterpF_F.H"

#include "NamespaceHeader.H"

IntFineInterp::IntFineInterp()
  :
  is_defined(false)
{
}

IntFineInterp::~IntFineInterp()
{
}

IntFineInterp::IntFineInterp(const DisjointBoxLayout& a_fine_domain,
			     const int&  a_numcomps,
			     const int& a_ref_ratio,
			     const IntVect& a_ghost,
			     const Box& a_fine_problem_domain)
  :
  is_defined(false)
{
  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_ghost, fineProbDomain);
}

IntFineInterp::IntFineInterp(const DisjointBoxLayout& a_fine_domain,
			     const int&  a_numcomps,
			     const int& a_ref_ratio,
			     const IntVect& a_ghost,
			     const ProblemDomain& a_fine_problem_domain)
  :
  is_defined(false)
{
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_ghost, a_fine_problem_domain);
}

void
IntFineInterp::define(const DisjointBoxLayout& a_fine_domain,
		      const int& a_numcomps,
		      const int& a_ref_ratio,
		      const IntVect& a_ghost,
		      const Box& a_fine_problem_domain)
{
  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_ghost, fineProbDomain);
}

void
IntFineInterp::define(const DisjointBoxLayout& a_fine_domain,
		      const int& a_numcomps,
		      const int& a_ref_ratio,
		      const IntVect& a_ghost,
		      const ProblemDomain& a_fine_problem_domain)
{
  CH_TIME("FineInterp::define");
  // check for consistency
  CH_assert (a_fine_domain.checkPeriodic(a_fine_problem_domain));
  m_ref_ratio = a_ref_ratio;
  m_coarse_problem_domain = coarsen(a_fine_problem_domain, m_ref_ratio);
  //
  // create the work array
  DisjointBoxLayout coarsened_fine_domain;
  coarsen ( coarsened_fine_domain,
            a_fine_domain,
            m_ref_ratio );
  m_fineGhost = a_ghost;
  IntVect crseGhost = m_fineGhost;
  crseGhost +=  m_ref_ratio * IntVect::Unit;
  crseGhost /= m_ref_ratio;
  
  m_coarsened_fine_data.define ( coarsened_fine_domain,
                                 a_numcomps,
                                 crseGhost);
  is_defined = true;
}

bool
IntFineInterp::isDefined() const
{
  return ( is_defined );
}

void IntFineInterp::pwcInterpToFine(LevelData<BaseFab<int> >& a_fine_data,
				    const LevelData<BaseFab<int> >& a_coarse_data) 
{
   CH_TIME("IntFineInterp::pwcinterpToFine");
   CH_assert(is_defined);

   // this doesnt' work on one machine but does on another.
   // I have no idea why...
   // Copier copier(a_coarse_data.disjointBoxLayout(),
   // 		 m_coarsened_fine_data.boxLayout(),
   // 		 m_coarse_problem_domain,
   // 		 m_coarsened_fine_data.ghostVect(),
   // 		 true);
	 
  
   // a_coarse_data.copyTo(a_coarse_data.interval(),
   // 			m_coarsened_fine_data,
   // 			m_coarsened_fine_data.interval(),
   // 			copier );

  a_coarse_data.copyTo(a_coarse_data.interval(),
			m_coarsened_fine_data,
			m_coarsened_fine_data.interval());


   const BoxLayout fine_domain = a_fine_data.boxLayout();
   for (DataIterator dit(fine_domain); dit.ok(); ++dit)
     {

       const Box& crseBox = m_coarsened_fine_data[dit].box();
       //Box crseBox = b;
       //crseBox.coarsen(m_ref_ratio);
       Box fineBox = crseBox;
       fineBox.refine(m_ref_ratio);
       BaseFab<int> tmp(fineBox, a_fine_data[dit].nComp());

       pwcinterpGridData(tmp,
			 m_coarsened_fine_data[dit],
			 crseBox,
			 m_ref_ratio);

       a_fine_data[dit].copy(tmp);

     }
   
}

void IntFineInterp::pwcinterpGridData(BaseFab<int>& a_fine,
				      const BaseFab<int>& a_coarse,
				      const Box& a_coarsened_fine_box,
				      int a_ref_ratio) const
{
  CH_TIME("IntFineInterp::pwcinterpGridData");
  // fill fine data with piecewise constant coarse data
  const Box& b = a_coarsened_fine_box;
  Box refbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);
  

  FORT_INTINTERPCONSTANT (CHF_FIA(a_fine),
  			  CHF_CONST_FIA(a_coarse),
  			  CHF_BOX(b),
  			  CHF_CONST_INT(a_ref_ratio),
  			  CHF_BOX(refbox));
  
}
#include "NamespaceFooter.H"
