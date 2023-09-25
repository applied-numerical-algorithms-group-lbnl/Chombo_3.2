#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AmrAdaptor.H"
#include "N_VectorOps_F.H"

#ifdef CH_MPI
#include <sundials/sundials_mpi_types.h>
#endif

// Define AMR data and masks for doing multi-level operations
void AmrAdaptor::define(Vector<LevelData<FArrayBox>*> amrData, 
    Vector<int> refRatio, bool ownData)
{
  CH_TIME("AmrAdaptor::define");
  this->m_ownData = ownData;
  size_t nlvl = amrData.size();
  CH_assert(nlvl > 0);
  // TODO more checks if paranoid
  m_amrData = amrData;
  CH_assert(refRatio.size() >= nlvl-1);
  m_refRatio = refRatio;
  m_comp = m_amrData[0]->nComp();

  // Build the mask
  m_amrMask = Vector<LevelData<FArrayBox>*>(nlvl,NULL);
  unsigned long masksum = 0;
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    m_amrMask[lvl] = new LevelData<FArrayBox>(dbl, 1, IntVect::Zero);
    LevelData<FArrayBox>& mask = *m_amrMask[lvl];

    // Loop over coarse dbl to set masks
    for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit) {
      // initialize this level mask to 1
      const Box& box = dbl[dit()];
      FArrayBox& fab = mask[dit()];
      fab.setVal(1.0);
      size_t boxsum = box.numPts();

      // For each box in this dbl, mask=0 on each coarsened fine dbl box
      if (lvl < nlvl-1) { // if there is a finer level
        // loop over ALL finer level boxes, processes -> layout iterator
        DisjointBoxLayout fdbl = m_amrData[lvl+1]->disjointBoxLayout();
        for (LayoutIterator flit= fdbl.layoutIterator(); flit.ok(); ++flit) {
          // set mask to 0 in intersection of any of coarse fine boxes
          Box maskb(fdbl[flit()]);
          maskb.coarsen(m_refRatio[lvl]);
          maskb &= box;
          if (!maskb.isEmpty())
          {
            fab.setVal(0.0, maskb, 0);
            boxsum -= maskb.numPts();
          }
        }
      }
      // Check with mask sum vs. length
      CH_assert(boxsum == fab.sum(0));
      pout() << "  Level " << lvl << ", box " << dbl[dit()] <<
        ", sum of mask: " << boxsum << endl;
      masksum += boxsum;
    }
  }

#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &masksum, 1, MPI_UNSIGNED_LONG, MPI_SUM, Chombo_MPI::comm);
#endif
  pout() << "Sum of masked level values: " << masksum << endl;
  m_length = masksum;
}

bool AmrAdaptor::copyTo(ChomboSundialsAdaptor& a)
{
  MayDay::Abort("AmrAdaptor::copyTo not implemented yet");
  return false;
}

void AmrAdaptor::exchange()
{
  MayDay::Abort("AmrAdaptor::exchange not implemented yet");
}

ChomboSundialsAdaptor* AmrAdaptor::newAdaptor()
{
  CH_TIME("AmrAdaptor::newAdaptor");
  AmrAdaptor* retval = new AmrAdaptor();
  size_t nlvl = m_amrData.size();
  CH_assert(nlvl > 0);
  retval->m_amrData = Vector<LevelData<FArrayBox>*>(nlvl,NULL);
  retval->m_amrMask = Vector<LevelData<FArrayBox>*>(nlvl,NULL);
  retval->m_refRatio = m_refRatio;
  retval->m_length = m_length;
  retval->m_comp = m_comp;
  retval->m_ownData = true; // will be responsible for deleting memory
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    IntVect ghostVect = ldf.ghostVect();
    retval->m_amrData[lvl] = new LevelData<FArrayBox>(dbl, m_comp, ghostVect);
    // Create and copy the mask data
    retval->m_amrMask[lvl] = new LevelData<FArrayBox>(dbl, 1);
    m_amrMask[lvl]->copyTo(*(retval->m_amrMask[lvl]));
  }

  return retval; 
}

// Return the total number of *uncovered* AMR cells X comps
unsigned long AmrAdaptor::getLength()
{
  size_t nlvl = m_amrData.size();
  CH_assert(nlvl > 0);
  return m_length*m_comp;
}

// Puts a*x+b*y in this instance's data
void AmrAdaptor::linearSum(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay, Real a, Real b)
{
  CH_TIME("AmrAdaptor::linearSum");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  AmrAdaptor* y = dynamic_cast<AmrAdaptor*>(&ay);
  CH_assert(y != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(y->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  CH_assert(m_length == y->m_length);
  // Add scale input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    LevelData<FArrayBox>& yldf = *(y->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      const FArrayBox& yfab = yldf[dit()];
      fab.axby(xfab, yfab, a, b); // all comps, ghost zones too
    }
  }
}

// Sets this instance's data to constant c
void AmrAdaptor::setConst(Real c)
{
  CH_TIME("AmrAdaptor::setConst");
  // Just set all the un/covered/ghost cells
  size_t maxlvl = m_amrData.size();
  for (size_t lvl=0; lvl < maxlvl; lvl++)
  {
    DataIterator dit = (*m_amrData[lvl]).dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      (*m_amrData[lvl])[dit()].setVal(c);
    }
  }
}

// Sets this instance's data to x*y on valid regions
void AmrAdaptor::prod(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay)
{
  CH_TIME("AmrAdaptor::prod");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  AmrAdaptor* y = dynamic_cast<AmrAdaptor*>(&ay);
  CH_assert(y != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(y->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  CH_assert(m_length == y->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add scale input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    LevelData<FArrayBox>& yldf = *(y->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      const FArrayBox& yfab = yldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.mult(yfab); 
    }
  }
}

// Sets this instance's data to x/y on valid regions
void AmrAdaptor::div(ChomboSundialsAdaptor& ax,
  ChomboSundialsAdaptor& ay)
{
  CH_TIME("AmrAdaptor::div");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  AmrAdaptor* y = dynamic_cast<AmrAdaptor*>(&ay);
  CH_assert(y != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(y->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  CH_assert(m_length == y->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add scale input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    LevelData<FArrayBox>& yldf = *(y->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      const FArrayBox& yfab = yldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.divide(yfab); 
    }
  }
} 

// Puts c*x in this instance's data
void AmrAdaptor::scale(ChomboSundialsAdaptor& ax, Real c)
{
  CH_TIME("AmrAdaptor::scale");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add scale input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.mult(c); 
    }
  }
}

// Puts abs(x) in this instance's data
void AmrAdaptor::abs(ChomboSundialsAdaptor& ax)
{
  CH_TIME("AmrAdaptor::abs");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add abs input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.abs(); 
    }
  }
}

// Puts 1/x in this instance's data
void AmrAdaptor::inv(ChomboSundialsAdaptor& ax)
{
  CH_TIME("AmrAdaptor::inv");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add 1 / input in all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.invert(1.0); 
    }
  }
}

// Puts x+b in this instance's data
void AmrAdaptor::addConst(ChomboSundialsAdaptor& ax, Real b)
{
  CH_TIME("AmrAdaptor::addConst");
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  bool inplace=(x == this); // check if ax pointed to this 
  // Add input + constant to all valid boxes, including covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      if (!inplace)
        fab.copy(xfab); 
      fab.plus(b); 
    }
  }
}

// Returns MPI-sum dot product across all ranks, sum(data^T * input)
// Masked for covered cells
Real AmrAdaptor::dotProd(ChomboSundialsAdaptor& ax)
{
  CH_TIME("AmrAdaptor::dotProd");
  Real dotProd=0;
  AmrAdaptor* x = dynamic_cast<AmrAdaptor*>(&ax);
  CH_assert(x != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(x->m_amrData.size() == nlvl);
  CH_assert(m_length == x->m_length);
  // Add dot this and input in all valid boxes, mask covered cells
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& xldf = *(x->m_amrData[lvl]);
    LevelData<FArrayBox>& mldf = *m_amrMask[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      const Box& box = dbl[dit()];
      const FArrayBox& fab = ldf[dit()];
      const FArrayBox& xfab = xldf[dit()];
      const FArrayBox& mask = mldf[dit()];
      FORT_MASKDOTPROD(CHF_CONST_FRA(fab), CHF_CONST_FRA(xfab), 
          CHF_CONST_FRA(mask), CHF_BOX(box), CHF_REAL(dotProd));
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &dotProd, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return dotProd;
}

// Do an MPI-reduce to return all ranks' data max norm
Real AmrAdaptor::maxNorm()
{
  CH_TIME("AmrAdaptor::maxNorm");
  Real maxNorm = 0;
  size_t maxlvl = m_amrData.size();
  for (size_t lvl=0; lvl < maxlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& mldf = *m_amrMask[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      const Box& b = dbl[dit()];
      const FArrayBox& fab = ldf[dit()];
      const FArrayBox& mask = mldf[dit()];
      FORT_MASKMAXNORM(CHF_CONST_FRA(fab), CHF_CONST_FRA(mask), CHF_BOX(b), 
          CHF_REAL(maxNorm));
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
#endif
  return maxNorm;
}

// Returns MPI-sum weighted RMS norm
Real AmrAdaptor::wRMSNorm(ChomboSundialsAdaptor& aw)
{
  CH_TIME("AmrAdaptor::wRMSNorm");
  Real sumwsq = 0;
  AmrAdaptor* w = dynamic_cast<AmrAdaptor*>(&aw);
  CH_assert(w != NULL);
  // quick check on length, number of levels
  size_t nlvl = m_amrData.size();
  CH_assert(w->m_amrData.size() == nlvl);
  CH_assert(m_length == w->m_length);
   for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& wldf = *(w->m_amrData[lvl]);
    LevelData<FArrayBox>& mldf = *m_amrMask[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      const Box& b = dbl[dit()];
      const FArrayBox& fab = ldf[dit()];
      const FArrayBox& wfab = wldf[dit()];
      const FArrayBox& mask = mldf[dit()];
      FORT_MASKWTDSQ(CHF_CONST_FRA(fab), CHF_CONST_FRA(wfab), 
          CHF_CONST_FRA(mask), CHF_BOX(b), CHF_REAL(sumwsq));
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &sumwsq, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return std::sqrt(sumwsq / getLength());
}

// Returns MPI-sum signed weighted RMS norm
Real AmrAdaptor::wRMSNormMask(ChomboSundialsAdaptor& aw,
    ChomboSundialsAdaptor& id)
{
  MayDay::Abort("AmrAdaptor::wRMSNormMask not implemented yet");
  return 0;
}

// Do an MPI-reduce to return all ranks' data min
Real AmrAdaptor::min()
{
  CH_TIME("AmrAdaptor::min");
  Real minVal = BASEFAB_REAL_SETVAL;
  size_t nlvl = m_amrData.size();
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& mldf = *m_amrMask[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      const Box& b = dbl[dit()];
      const FArrayBox& fab = ldf[dit()];
      const FArrayBox& mask = mldf[dit()];
      FORT_MASKMIN(CHF_CONST_FRA(fab), CHF_CONST_FRA(mask), CHF_BOX(b), 
          CHF_REAL(minVal));
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &minVal, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
#endif
  return minVal;
}


// Do an MPI-sum to return all-ranks weighted L2 norm
Real AmrAdaptor::wL2Norm(ChomboSundialsAdaptor& aw)
{
  MayDay::Abort("AmrAdaptor::wL2Norm not implemented yet");
  return 0;
}

// Do an MPI-sum to return all ranks' p-norm sums across components,
// masking AMR-covered cells and without any dx / refinement ratio weighting
Real AmrAdaptor::l1Norm()
{
  CH_TIME("AmrAdaptor::l1Norm");
  Real sumxp = 0;
  size_t nlvl = m_amrData.size();
  for (size_t lvl=0; lvl < nlvl; lvl++)
  {
    LevelData<FArrayBox>& ldf = *m_amrData[lvl];
    LevelData<FArrayBox>& mldf = *m_amrMask[lvl];
    DisjointBoxLayout dbl = ldf.disjointBoxLayout();
    DataIterator dit = ldf.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
      const Box& b = dbl[dit()];
      const FArrayBox& x = ldf[dit()];
      const FArrayBox& mask = mldf[dit()];
      Real fabsumxp = 0;
      int p = 1;
      FORT_MASKSUMXP(CHF_CONST_FRA(x), CHF_CONST_FRA(mask), 
          CHF_BOX(b), CHF_REAL(fabsumxp), CHF_INT(p));
      sumxp += fabsumxp;
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &sumxp, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  return sumxp;
}

// Return 1 if x>=b, 0 otherwise, in this instance's data
void AmrAdaptor::compare(ChomboSundialsAdaptor& ax, Real b)
{
  MayDay::Abort("AmrAdaptor::compare not implemented yet");
}

// Return 1 if x>=b, 0 otherwise, in this instance's data
bool AmrAdaptor::invTest(ChomboSundialsAdaptor& ax)
{
  MayDay::Abort("AmrAdaptor::invTest not implemented yet");
  return false;
}

bool AmrAdaptor::constrMask(ChomboSundialsAdaptor& ac,
    ChomboSundialsAdaptor& am)
{
  MayDay::Abort("AmrAdaptor::constrMask not implemented yet");
  return false;
}

Real AmrAdaptor::minQuotient(ChomboSundialsAdaptor& adenom)
{
  MayDay::Abort("AmrAdaptor::minQuotient not implemented yet");
  return 0;
}

// Print this to stdout
void AmrAdaptor::print()
{
  MayDay::Error("AmrAdaptor::print not implemented");
}

void AmrAdaptor::printFile(FILE* outfile)
{
  MayDay::Error("AmrAdaptor::printFile not implemented");
}

#if CH_USE_HDF5
void AmrAdaptor::printFileHDF(const char* filename)
{
  MayDay::Error("AmrAdaptor::printFileHDF not implemented");
}
#endif

