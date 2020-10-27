#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELDATAI_H_
#define _LEVELDATAI_H_

#include <cstdlib>
#include <algorithm>
using std::sort;

#include "parstream.H"
#include "CH_Timer.H"
#include <float.h>
#include "NamespaceHeader.H"

// 
// LevelData<FluxBox> specializations

template < > void LevelData<FluxBox>::degenerateLocalOnly( LevelData<FluxBox>& a_to, const SliceSpec& a_ss ) const;

#ifdef CH_USE_PROTO

////////////////////////////////////////////////////////////////
//Utility functions for the interface functions
IntVect PointToIntVect(const ::Proto::Point& a_pt)
{
    IntVect iv;
    for(int i=0; i<CH_SPACEDIM; i++)
        iv.setVal(i,a_pt[i]);
    return iv;
}

::Proto::Point IntVectToPoint(const IntVect& a_iv)
{
    Point pt;
    for(int i=0; i<CH_SPACEDIM; i++)
        pt[i]=iv[i];
    return pt;
}

/**
 * @brief Generate a list of boxes and the corresponding processor ids from the Proto DisjointBoxLayout
 * @param a_dbl a defined disjoint box layout
 * @param a_boxes list of all boxes in a_dbl (filled by this function)
 * @param a_proc_ids list of processor ids corresponding to boxes in a_boxes. a_boxes[i] belongs to a_proc_ids[i] (filled by this function)
 */
void getBoxesAndProcessors(::Proto::DisjointBoxLayout& a_dbl, Vector<Box>& a_boxes, Vector<int>& a_proc_ids)
{
    ::Proto::Point boxSize=a_dbl.boxSize();
    const shared_ptr<::Proto::DBLInternals> dblInternals=a_dbl.dblInternals();
    const vector<pair<::Proto::Point,int>>& allBoxes=dblInternals.allBoxes();
    for(auto ptProcPair: allBoxes)
    {
        a_proc_ids.push_back(ptProcPair.second);
        const ::Proto::Point& pt=ptProcPair.first;
        const ::Proto::Point lo(pt*boxSize);
        const ::Proto::Point hi((pt+::Proto::Point::Ones())*boxSize-::Proto::Point::Ones());
        a_boxes.emplace_back(PointToIntVect(lo),PointToIntVect(hi));
    }
}

/**
 * @brief Copies src into dest. src and dest must have the same number of components and live on the same box.
 * This function is currently SLOW.
 * @param src source of copy
 * @param dest box data that lives on CPU (host)
 */
template<unsigned int C>
void copyFArrayBoxToBoxData(const FArrayBox& src, BoxData<Real,C>& dest)
{
    //This is a slow version because it calls into accessor functions instead of working with the raw buffers.
    //We can write a faster version that uses the raw buffers if we know the layout of the two data buffers
    ::Proto::BoxIterator bit(dest.box());
    for(bit.begin(); !bit.done(); ++bit)
    {
        Point& pt=*bit;
        IntVect iv=PointToIntVect(pt);
        for(int ic=0; ic<C; ic++)
            *dest.data(pt,ic)=src(iv,ic);
    }
}

/**
 * @brief Copies src into dest. src and dest must have the same number of components and live on the same box.
 * This function is currently SLOW.
 * @param src source of copy. Lives on CPU (host)
 * @param dest destination of copy (filled by function)
 */
template<unsigned int C>
void copyBoxDataToFArrayBox(const BoxData<Real,C>& src, FArrayBox& dest)
{
    ::Proto::BoxIterator bit(src.box());
    for(bit.begin(); !bit.done(); ++bit)
    {
        Point& pt=*bit;
        IntVect iv=PointToIntVect(pt);
        for(int ic=0; ic<C; ic++)
            dest(pt,ic)=*src.data(iv,ic);
    }
}


///////////////////////////////////////////////////////////////
//Interface functions
/**
 * @brief Define a LevelData ld with the same box layout as the LevelBoxData lbd.
 * In particular, after this function is called, ld and lbd will be defined on the same set of boxes and
 * the same distribution of boxes across processors.
 * @param ld undefined LevelData. After the function call, ld will be defined and have the following properties
 *   1) same box layout as lbd
 *   2) same number of components as lbd
 *   3) same number of ghost cells (in each direction) as lbd
 * @param lbd defined LevelBoxData . Note: lbd lives on the host and has cell-centered data.
 */
template<unsigned int C>
void LevelData<FArrayBox>::define(const ::Proto::LevelBoxData<Real,C>& a_lbd)
{
    getBoxesAndProcessors(lbd->getDBL(),boxes,proc_ids);
    ld.define(DisjointBoxLayout(boxes,proc_ids),C,lbd->getGhostVector());
    ld.setProtoDBL(lbd->getDBL());
}

/**
 * @brief Copy the contents (valid+ghost data) of a_lbd into this object.
 * This object must be defined using the same Proto DisjointBoxLayout object as a_lbd (i.e. getProtoDBL()==a_lbd.getDBL()) 
 * It must also have C components and the same number of ghost cells as a_lbd
 * @param a_lbd LevelBoxData that lives on the host and contains cell-centered data
 */
template<unsigned int C>
void LevelData<FArrayBox>::copyFrom(const ::Proto::LevelBoxData<Real,C>& a_lbd)
{/**FILL IN**/}

/**
 * @brief Copy the contents (valid+ghost data) of this object into a_lbd.
 * This object must be defined using the same Proto DisjointBoxLayout object as a_lbd (i.e. getProtoDBL()==a_lbd.getDBL()) 
 * It must also have C components and the same number of ghost cells as a_lbd
 * @param a_lbd LevelBoxData that lives on the host and contains cell-centered data
 */
template<unsigned int C>
void LevelData<FArrayBox>::copyTo(const ::Proto::LevelBoxData<Real,C>& a_lbd)
{
    CH_assert(src.getProtoDBL()==dest.getDBL());
    //TODO: Put checks that src and dest have the same number of components and the same number of ghost cells
    DataIterator dit_ld(src.disjointBoxLayout());
    ::Proto::DataIterator dit_lbd(dest.getDBL());
    for(dit_lbd.begin(); *dit_lbd!=dit_lbd.end(); ++dit_lbd)
    {
        DataIndex di_ld=dit_ld[dit_lbd.getInt()]; //TODO: Check this: it may not be right because there are indices related to the local processor (the index that matters for looping through box data) and the indices related to the global box layout (all the boxes over all the processors)
        ::Proto::BoxData<Real,C>& dest_bd=dest[*dit_lbd];
        FArrayBox& src_bf=src[di_ld];
        copyFArrayBoxToBoxData(src_bf,dest_bd);
    }
}
#endif

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::LevelData()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::~LevelData()
{
  CH_TIME("~LevelData");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::LevelData(const DisjointBoxLayout& dp, int comps, const IntVect& ghost,
                        const DataFactory<T>& a_factory)
  : m_disjointBoxLayout(dp), m_ghost(ghost)
{
#ifdef CH_MPI
  this->numSends = 0;
  this->numReceives = 0;
#endif
  this->m_boxLayout = dp;
  this->m_comps = comps;
  this->m_isdefined = true;
  //no need for exchange copier if no ghost
  if(m_ghost != IntVect::Zero)
    {
      m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
    }

  if (!dp.isClosed())
    {
      MayDay::Error("non-disjoint DisjointBoxLayout: LevelData<T>::LevelData(const DisjointBoxLayout& dp, int comps)");
    }

  this->m_threadSafe = a_factory.threadSafe();
  //this->m_threadSafe = false;

  Interval interval(0, comps-1);
  this->allocateGhostVector(a_factory, ghost);
  this->setVector(*this, interval, interval); // Does nothing.
}
//-----------------------------------------------------------------------

// Since I need to thwart the user from invoking the
// 'define' methods that use a general BoxLayout, I cannot
// get at said functions myself now. Ha!  So, I have to recode
// them here.

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const DisjointBoxLayout& dp, int comps, const IntVect& ghost,
                          const DataFactory<T> & a_factory)
{
  CH_TIME("LevelData<T>::define(dbl,comps,ghost,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }

  this->m_isdefined = true;
  if (!dp.isClosed())
    {
      MayDay::Error("non-disjoint DisjointBoxLayout: LevelData<T>::define(const DisjointBoxLayout& dp,....)");
    }
  if (comps<=0)
    {
      MayDay::Error("LevelData::LevelData(const BoxLayout& dp, int comps)  comps<=0");
    }
  this->m_comps = comps;
  this->m_boxLayout = dp;

  this->m_threadSafe = a_factory.threadSafe();
  //this->m_threadSafe = false;

  m_disjointBoxLayout = dp;
  m_ghost = ghost;
  //no need for exchange copier if no ghost
  //do not do this if it is periodic because bad things can happen
  //with large number of ghost cells.
  bool periodic = m_disjointBoxLayout.physDomain().isPeriodic();
  if((m_ghost != IntVect::Zero) && !periodic)
    {
      m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
    }

  // Interval interval(0, comps-1);
  this->allocateGhostVector(a_factory, ghost);
  //  this->setVector(*this, interval, interval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const LevelData<T>& da,  const DataFactory<T> & a_factory)
{
  CH_TIME("LevelData<T>::define(LevelData<T>,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }
  this->m_isdefined = true;
  if (this == &da) return;
  m_disjointBoxLayout = da.m_disjointBoxLayout;
  this->m_boxLayout  = da.m_disjointBoxLayout;
  this->m_comps     = da.m_comps;
  this->m_threadSafe = a_factory.threadSafe();
  //this->m_threadSafe = false;

  m_ghost     = da.m_ghost;

  //no need for exchange copier if no ghost
  if(m_ghost != IntVect::Zero)
    {
      m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
    }
  this->allocateGhostVector(a_factory, m_ghost);
  da.copyTo(*this);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const LevelData<T>& da, const Interval& comps,
                          const DataFactory<T>& a_factory)
{
  CH_TIME("LevelData<T>::define(LevelData<T>,comps,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }
  this->m_isdefined = true;
  this->m_threadSafe = a_factory.threadSafe();
  //this->m_threadSafe = false;
  if (this == &da)
  {
    MayDay::Error(" LevelData<T>::define(const LevelData<T>& da, const Interval& comps) called with 'this'");
  }
  CH_assert(comps.size()>0);
  // this line doesn't make any sense!
  //CH_assert(comps.end()<=this->m_comps);
  CH_assert(comps.begin()>=0);

  m_disjointBoxLayout = da.m_disjointBoxLayout;
  this->m_boxLayout  = da.m_disjointBoxLayout;

  this->m_comps = comps.size();

  m_ghost = da.m_ghost;
  //no need for exchange copier if no ghost
  if(m_ghost != IntVect::Zero)
    {
      m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
    }

  Interval dest(0, this->m_comps-1);

  this->allocateGhostVector(a_factory, m_ghost);

  this->setVector(da, comps, dest);

}
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          BoxLayoutData<T>& dest,
                          const Interval& destComps) const
{
  CH_TIME("copyTo");
  if ((BoxLayoutData<T>*)this == &dest) return;
  
  if (this->boxLayout() == dest.boxLayout())
    {
      CH_TIME("local copy");
      // parallel direct copy here, no communication issues
      DataIterator it = this->dataIterator();
      int nbox = it.size();
#pragma omp parallel for if(this->m_threadSafe)
      for (int ibox = 0; ibox < nbox; ibox++)
        {
          dest[it[ibox]].copy(this->box(it[ibox]),
                          destComps,
                          this->box(it[ibox]),
                          this->operator[](it[ibox]),
                          srcComps);
        }
      return;
    }

  Copier copier(m_disjointBoxLayout, dest.boxLayout());
  copyTo(srcComps, dest, destComps, copier);
}

template<class T>
void LevelData<T>::localCopyTo(const Interval& srcComps,
                               LevelData<T>& dest,
                               const Interval& destComps) const
{ CH_TIME("localCopyTo");
  if ((BoxLayoutData<T>*)this == &dest) return;

  if(this->disjointBoxLayout() == dest.disjointBoxLayout())
    {
      DataIterator it = this->dataIterator();
      int nbox = it.size();
#pragma omp parallel for if(this->m_threadSafe)
      for (int ibox = 0; ibox < nbox; ibox++)
        {
          Box srcBox = this->box(it[ibox]);
          Box dstBox = this->box(it[ibox]);
          srcBox.grow(this->ghostVect());
          dstBox.grow(dest.ghostVect());
          Box minBox = srcBox;
          minBox &= dstBox;
          dest[it[ibox]].copy(minBox,
                          destComps,
                          minBox,
                          this->operator[](it[ibox]),
                          srcComps);
        }
    }
  else
    {
      MayDay::Error("localCopyTo only works with identical DBLs");
    }

  return;
}

template<class T>
void LevelData<T>::localCopyTo(LevelData<T>& dest) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->localCopyTo(this->interval(), dest, dest.interval());
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(BoxLayoutData<T>& dest) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          LevelData<T>& dest,
                          const Interval& destComps) const
{
  if (this == &dest)
  {
    MayDay::Error("src == dest in copyTo function. Perhaps you want exchange ?");
  }

  if (this->boxLayout() == dest.boxLayout()  && dest.ghostVect() == IntVect::Zero)
    {
      // parallel direct copy here, no communication issues
      DataIterator it = this->dataIterator();
      int nbox = it.size();
#pragma omp parallel for if(this->m_threadSafe)
      for (int ibox = 0; ibox < nbox; ibox++)
        {
          dest[it[ibox]].copy(this->box(it[ibox]),
                          destComps,
                          this->box(it[ibox]),
                          this->operator[](it[ibox]),
                          srcComps);
        }
      return;
    }

  Copier copier(m_disjointBoxLayout, dest.getBoxes(), dest.m_ghost);
  copyTo(srcComps, dest, destComps, copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(LevelData<T>& dest) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          BoxLayoutData<T>& dest,
                          const Interval& destComps,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_TIME("copyTo");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier copyTo");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif

  this->makeItSo(srcComps, *this, dest, destComps, copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(BoxLayoutData<T>& dest,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval(), copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          LevelData<T>& dest,
                          const Interval& destComps,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_TIME("copyTo");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier copyTo");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  this->makeItSo(srcComps, *this, dest, destComps, copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(LevelData<T>& dest,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval(), copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Interval& comps)
{
  CH_TIME("exchange+copier");
  // later on we can code this part as a direct algorithm
  // by copying and pasting the code from the Copier::define code
  // for now, just do the easy to debug approach.
  if(m_ghost != IntVect::Zero)//no need for exchange copier if no ghost
    {
      if(!m_exchangeCopier.isDefined())
        {
          CH_TIME("defining_copier");
          m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
        }
        {
          CH_TIME("actual_exchange");
          exchange(comps, m_exchangeCopier);
        }
    }

  // if there aren't any ghost cells, there isn't really anything
  // to do here (also, if m_ghost == Zero, m_exchangeCopier
  // wasn't defined!
  //if (m_ghost != IntVect::Zero)
  //this->makeItSo(comps, *this, *this, comps, m_exchangeCopier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(void)
{
  exchange(this->interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Interval& comps,
                            const Copier& copier,
                            const LDOperator<T>& a_op)
{
  CH_TIME("exchange");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier exchange");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  this->makeItSo(comps, *this, *this, comps, copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Copier& copier)
{
  exchange(this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeNoOverlap(const Copier& copier)
{
  CH_TIME("exchangeNoOverlap");
  this->makeItSoBegin(this->interval(), *this, *this, this->interval(), copier);
  this->makeItSoEnd(*this, this->interval());
  this->makeItSoLocalCopy(this->interval(), *this, *this, this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeBegin(const Copier& copier)
{
  CH_TIME("exchangeBegin");
  this->makeItSoBegin(this->interval(), *this, *this, this->interval(), copier);
  this->makeItSoLocalCopy(this->interval(), *this, *this, this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeEnd()
{
  CH_TIME("exchangeEnd");
  this->makeItSoEnd(*this, this->interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayout& dp, int comps,  const DataFactory<T>& a_factory)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayout& dp)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayoutData<T>& da, const DataFactory<T>& a_factory )
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayoutData<T>& da, const Interval& comps,
                          const DataFactory<T>& a_factory)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::apply( void (*a_Function)(const Box& box, int comps, T& t) )
{
  DataIterator it = this->dataIterator();
  int nbox = it.size();
#pragma omp parallel for
  for (int ibox = 0; ibox < nbox; ibox++)
  {

    a_Function(m_disjointBoxLayout.get(it[ibox]), this->m_comps, *(this->m_vector[it[ibox].datInd()]));
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::apply( const ApplyFunctor & f )
{
  DataIterator it = this->dataIterator();
  int nbox = it.size();
#pragma omp parallel for
  for (int ibox = 0; ibox < nbox; ibox++)
  {
      //     unsigned int index = this->m_boxLayout.lindex(it());
    f(m_disjointBoxLayout.get(it[ibox]), this->m_comps, *(this->m_vector[it[ibox].datInd()]));
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T> void
LevelData<T>::degenerate( LevelData<T>& a_to,
                          const SliceSpec& a_sliceSpec ) const
{
    DisjointBoxLayout toDBL;
    m_disjointBoxLayout.degenerate( toDBL, a_sliceSpec );
    IntVect toGhost;
    for ( int i=0;i<CH_SPACEDIM;++i )
    {
        if ( i != a_sliceSpec.direction )
        {
            toGhost[i] = m_ghost[i];
        } else
        {
            toGhost[i] = 0;
        }
    }
    a_to.define( toDBL, this->nComp(), toGhost );
    copyTo( a_to );
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T> void
LevelData<T>::degenerateLocalOnly( LevelData<T>& a_to,
                                   const SliceSpec& a_sliceSpec ) const
{
    DisjointBoxLayout toDBL;
    m_disjointBoxLayout.degenerate( toDBL, a_sliceSpec, true );
    IntVect toGhost;
    for ( int i=0;i<CH_SPACEDIM;++i )
    {
        if ( i != a_sliceSpec.direction )
        {
            toGhost[i] = m_ghost[i];
        } else
        {
            toGhost[i] = 0;
        }   
    }
    a_to.define( toDBL, this->nComp(), toGhost );

    // manage copyTo ourselves to maintain locality
    DataIterator fromDit = this->dataIterator();
    DataIterator toDit = a_to.dataIterator();
    for (toDit.begin(); toDit.ok(); ++toDit)
      {
        fromDit.begin();
        bool done = false;
        while (!done)
          {
            if (m_disjointBoxLayout[fromDit].contains(toDBL[toDit]))
              {
                // boxes intersect, do copy
                FArrayBox& toFAB = a_to[toDit];
                const FArrayBox& fromFAB = this->operator[](fromDit);
                // note that we're including ghost cells here
                Box intersectBox = toFAB.box();
                intersectBox &= fromFAB.box();
                // now do copy
                toFAB.copy(fromFAB, intersectBox);
                done = true;
              } // end if we found an intersection
            else
              {
                ++fromDit;
                // probably want to check for lexigraphical sorting done-ness

                // are we done looping through fromBoxes?
                if (!fromDit.ok()) done = true;
              } // end if we didn't find a match
          } // end while loop over fromBoxes
          
      } // end loop over toBoxes
    
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
#endif
