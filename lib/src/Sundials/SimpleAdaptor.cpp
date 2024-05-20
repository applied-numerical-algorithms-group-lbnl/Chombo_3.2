#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "SimpleAdaptor.H"
#include <cassert>

#include "DebugOut.H"
#include "AMRIO.H"

unsigned long SimpleAdaptor::getLength()
{
  return m_dp.numCells()*m_nComp;
}

// Define from a pointer, 
void SimpleAdaptor::define(LevelData<FArrayBox>* lvlData, bool ownData)
{
  m_data = lvlData;
  m_dp = lvlData->disjointBoxLayout();
  m_nComp = lvlData->nComp();
  m_ghost = lvlData->ghostVect();
  m_ownData = ownData;
}

// Puts a*x+b*y in this instance's data
void SimpleAdaptor::linearSum(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay, Real a, Real b)
{
  SimpleAdaptor* sax = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sax != NULL);
  SimpleAdaptor* say = dynamic_cast<SimpleAdaptor*>(&ay);
  assert(say != NULL);
  assert(m_dp == sax->m_dp);
  assert(m_dp == say->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& xfab = (*sax->m_data)[dit()];
    const FArrayBox& yfab = (*say->m_data)[dit()];
    fab.axby(xfab, yfab, a, b); // all comps, ghost zones too
  }
}

// Sets this instance's data to constant c
void SimpleAdaptor::setConst(Real c)
{
  DataIterator dit = m_dp.dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
    (*m_data)[dit()].setVal(c);
  }
}

// Sets this instance's data to x*y on valid regions
void SimpleAdaptor::prod(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay)
{
  SimpleAdaptor* sax = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sax != NULL);
  SimpleAdaptor* say = dynamic_cast<SimpleAdaptor*>(&ay);
  assert(say != NULL);
  bool xinplace=(sax == this); // check if ax pointed to this 
  assert(m_dp == sax->m_dp);
  assert(m_dp == say->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& xfab = (*sax->m_data)[dit()];
    const FArrayBox& yfab = (*say->m_data)[dit()];
    if (!xinplace)
      fab.copy(xfab); 
    fab.mult(yfab);
  }
}

// Sets this instance's data to x/y on valid regions
void SimpleAdaptor::div(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay)
{
  SimpleAdaptor* sax = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sax != NULL);
  SimpleAdaptor* say = dynamic_cast<SimpleAdaptor*>(&ay);
  assert(say != NULL);
  bool xinplace=(sax == this); // check if ax pointed to this 
  assert(m_dp == sax->m_dp);
  assert(m_dp == say->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& xfab = (*sax->m_data)[dit()];
    const FArrayBox& yfab = (*say->m_data)[dit()];
    if (!xinplace)
      fab.copy(xfab); 
    fab.divide(yfab); 
  }
} 

// Puts c*x in this instance's data
void SimpleAdaptor::scale(ChomboSundialsAdaptor& ax, Real c)
{
  SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sa != NULL);
  bool inplace=(sa == this); // check if ax pointed to this 
  assert(m_dp == sa->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& afab = (*sa->m_data)[dit()];
    if (!inplace)
      fab.copy(afab); 
    fab.mult(c); 
  }
}

// Puts abs(x) in this instance's data
void SimpleAdaptor::abs(ChomboSundialsAdaptor& ax)
{
  SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sa != NULL);
  bool inplace=(sa == this); // check if ax pointed to this 
  assert(m_dp == sa->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& afab = (*sa->m_data)[dit()];
    if (!inplace)
      fab.copy(afab); 
    fab.abs(); 
  }
}

// Puts 1/x in this instance's data
void SimpleAdaptor::inv(ChomboSundialsAdaptor& ax)
{
  SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sa != NULL);
  bool inplace=(sa == this); // check if ax pointed to this 
  assert(m_dp == sa->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& afab = (*sa->m_data)[dit()];
    if (!inplace)
      fab.copy(afab); 
    fab.invert(1.0); 
  }
  // Note - no check!
}

// Puts x+b in this instance's data
void SimpleAdaptor::addConst(ChomboSundialsAdaptor& ax, Real b)
{
  SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sa != NULL);
  assert(m_dp == sa->m_dp);
  bool inplace=(sa == this); // check if ax pointed to this
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& afab = (*sa->m_data)[dit()];
    if (!inplace)
      fab.copy(afab);
    fab.plus(b);
  }
}

// Returns MPI-sum dot product across all ranks, sum(data^T * input)
Real SimpleAdaptor::dotProd(ChomboSundialsAdaptor& ax)
{
  Real dotProd=0;
  SimpleAdaptor* sa = dynamic_cast<SimpleAdaptor*>(&ax);
  assert(sa != NULL);
  assert(m_dp == sa->m_dp);
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& box = m_dp[dit()];
    FArrayBox& fab = (*m_data)[dit()];
    const FArrayBox& xfab = (*sa->m_data)[dit()];
    dotProd += fab.dotProduct(xfab, box); // all components?
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &dotProd, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return dotProd;
}

// Do an MPI-reduce to return all ranks' data max norm
Real SimpleAdaptor::maxNorm()
{
  Real maxNorm = -BASEFAB_REAL_SETVAL;
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& b = m_dp[dit()];
    const FArrayBox& fab = (*m_data)[dit()];
    maxNorm = max(fab.norm(b, 0, 0, m_nComp), maxNorm);
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
#endif
  return maxNorm;
}

// Returns MPI-sum weighted RMS norm
Real SimpleAdaptor::wRMSNorm(ChomboSundialsAdaptor& aw)
{
  Real sumwsq = 0;
  SimpleAdaptor* saw = dynamic_cast<SimpleAdaptor*>(&aw);
  assert(saw != NULL);
  assert(m_dp == saw->m_dp);
  DataIterator dit = m_dp.dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
    // ignore ghost cells
    const Box& box = m_dp[dit()];
    FArrayBox tmp(box, m_nComp); // TODO - get rid of this allocated temp?
    const FArrayBox& wfab = (*saw->m_data)[dit()];
    const FArrayBox& fab = (*m_data)[dit()];
    Real localsumwsq = 0;
    tmp.copy(wfab, box);
    tmp.mult(fab, box, 0, 0, m_nComp); // w*fab
    tmp.mult(tmp, box, 0, 0, m_nComp); // (w*fab)^2
    localsumwsq += tmp.sum(box, 0, m_nComp); // sum (w*fab)^2
    sumwsq += localsumwsq;
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &sumwsq, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return std::sqrt(sumwsq / getLength());
}

// Returns MPI-sum signed weighted RMS norm
Real SimpleAdaptor::wRMSNormMask(ChomboSundialsAdaptor& aw,
    ChomboSundialsAdaptor& id)
{
  MayDay::Abort("SimpleAdaptor::wRMSNormMask not implemented yet");
  return 0;
}

// Do an MPI-reduce to return all ranks' data min
Real SimpleAdaptor::min()
{
  Real minVal = BASEFAB_REAL_SETVAL;
  DataIterator dit = m_dp.dataIterator();
  dit.reset();
  for (dit.reset(); dit.ok(); ++dit) {
    const Box& b = m_dp[dit()];
    const FArrayBox& fab = (*m_data)[dit()];
    for (int comp=0; comp<m_nComp; comp++)
      minVal = std::min(fab.min(b, comp), minVal);
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &minVal, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
#endif
  return minVal;
}

// Do an MPI-sum to return all-ranks weighted L2 norm
Real SimpleAdaptor::wL2Norm(ChomboSundialsAdaptor& aw)
{
  MayDay::Abort("SimpleAdaptor::weightedL2Norm not implemented yet");
  return 0;
}

// Do an MPI-sum to return all ranks' p-norm
Real SimpleAdaptor::l1Norm()
{
  Real norm = 0;
  DataIterator dit = m_dp.dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
    // ignore ghost cells
    const Box& region = m_dp[dit()];
    norm += (*m_data)[dit()].norm(region, 1, 0, m_nComp);
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return norm;
}

// Return 1 if x>=b, 0 otherwise, in this instance's data
void SimpleAdaptor::compare(ChomboSundialsAdaptor& ax, Real b)
{
  MayDay::Abort("SimpleAdaptor::compare not implemented yet");
}

// Return 1 if x>=b, 0 otherwise, in this instance's data
bool SimpleAdaptor::invTest(ChomboSundialsAdaptor& ax)
{
  MayDay::Abort("SimpleAdaptor::invTest not implemented yet");
  return false;
}

// Not sure what this does, but requires an MPI-reduce min
bool SimpleAdaptor::constrMask(ChomboSundialsAdaptor& ac,
    ChomboSundialsAdaptor& am)
{
  MayDay::Abort("SimpleAdaptor::constrMask not implemented yet");
  return false;
}

// Not sure what this does, but requires an MPI-reduce min
Real SimpleAdaptor::minQuotient(ChomboSundialsAdaptor& adenom)
{
  MayDay::Abort("SimpleAdaptor::minQuotient not implemented yet");
  return 0;
}

// Print this to stdout
void SimpleAdaptor::print()
{
  dumpLDFPar(m_data);
}

void SimpleAdaptor::printFile(FILE* outfile)
{
  MayDay::Error("SimpleAdaptor::printFile not implemented");
}

#if CH_USE_HDF5
void SimpleAdaptor::printFileHDF(const char* filename)
{
  writeLevelname(m_data, filename);
}
#endif
