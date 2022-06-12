#pragma once
#include "Lapack.H"

//-----------------------------------------------------------------------------
// There are lots of flavors of batch 1D solvers that we will need to consider
// the whole big foundation of non-linear, PDI, compact, EB generalities
// - Batch vs. not
// - 3- or 5- diag? 1-, 2-, or 3D?
// - Linear (factor once / reuse)
// - All dirs vs. 1 or 2 sub-dirs
// - Back ends (LAPACK, MKL, tri-gpu, etc.) and formats ("compact")
// - Mix-and-match all above?
// - Arch / platform, affects libs
// - Multi-component or single or interval
// - Callbacks to expose RHS / diag setting? Or just inheritance
// - Brick extensions?
//-----------------------------------------------------------------------------

class BatchSolver
{
public:
  virtual void copyFABRhsDir(const int l, const FArrayBox& rhs) {};
  virtual void setFABBatchDiagsDir(const int l) {};
  virtual void solveFABBatchDir(const int l) {};
  virtual void copyFABSolnDir(const int l, FArrayBox& soln) {};

// private:
public:
  int m_verbosity; // 0 - no output, 1 - high-level, 2 - details, 3 - lots
 };


class VC2lineGSRB : public BatchSolver
{
public:
  // Note, if a_dx is not provided unit dx is used
  VC2lineGSRB(Box a_b, RealVect a_dx=RealVect::Unit) : m_b(a_b), m_dx(a_dx)
  {
    assert(SpaceDim == 3); // Only works for 3D

    // TODO - generalize this, dim, # comp?
    m_diags.define(m_b,3); // diags, sub/comp 0, diag/comp 1, sup/comp 2
    m_tmp.define(m_b,1); // tmp for the LAPACK rhs's --> soln's

    // Get pointers for diags, rhs
    // NOTE: these are just chunks of memory, with l-major indexing
    // --> Don't index by [iv]
    m_sub = m_diags.dataPtr(0);
    m_diag = m_diags.dataPtr(1); 
    m_sup = m_diags.dataPtr(2); 
    m_rhs = m_tmp.dataPtr(0); 
  }

  // Does a line relaxation for this op in direction d:
  // L phi = alpha * acoef_i * phi_i - beta * 
  //         \sum_d(  b_d_i+1/2 * (phi_i+1 - phi_i)/dx_d^2
  //                - b_d_i-1/2 * (phi_i - phi_i-1)/dx_d^2 )
  void lineRelaxRB(const int a_d, FArrayBox& a_phi, const FArrayBox& a_rhs,
      const RealVect& a_dx_vect, const Real& a_alpha, 
      const FArrayBox& a_aCoef, const Real& a_beta, 
      const FluxBox& bCoef, const int whichPass)
  {
    CH_TIME("VC2lineGSRB::lineRelaxRB");
    assert(a_rhs.box().contains(m_b));

    int N = m_b.size(a_d); // size in dir of solve
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(a_d)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction = " << a_d << ", size = " << N << endl;

    // Holds which direction is 0,1,2
    IntVect id(a_d % 3, (a_d+1) % 3, (a_d+2) % 3);
    RealVect dx(a_dx_vect[id[0]], a_dx_vect[id[1]], a_dx_vect[id[2]]);
    IntVect lsz(bsize[id[0]], bsize[id[1]], bsize[id[2]]);
    if (m_verbosity >= 2)
    {
      pout() << "  Looping over " << bloop << endl;
      pout() << "  Shifted sizes " << lsz << endl;
    }
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      if ((iv.sum() % 2) != whichPass) // Do all R then B
        continue; 
      IntVect liv(iv[id[0]], iv[id[1]], iv[id[2]]);
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      Real dx2inv = 1./(dx[0]*dx[0]);
      for (int i=0; i < N; ++i)
      {
        IntVect ix = iv + i*BASISV(a_d);
        // 1. Copy rhs into direction rhs
        m_rhs[offset + i] = a_rhs(ix,0); // TODO - more than 1 comp?

        // L phi = alpha * acoef_i * phi_i - beta * 
        //         \sum_d(  b_d_i+1/2 * (phi_i+1 - phi_i)/dx_d^2
        //                - b_d_i-1/2 * (phi_i - phi_i-1)/dx_d^2 )
        // 2. Copy system into diagonals
        Real sup = -a_beta*bCoef[id[0]](ix+BASISV(a_d),0)*dx2inv; // i+1 coef
        Real sub = -a_beta*bCoef[id[0]](ix,0)*dx2inv; // i-1 coef
        Real diag = a_alpha * a_aCoef(ix,0) - sup - sub; // i coef
        m_sub[offset + i] = sub;
        m_diag[offset + i] = diag;
        m_sup[offset + i] = sup;
      }
      // Solve

      /*
      // Copy phi back
      for (int i=0; i < N; ++i)
      {
        IntVect ix = iv + i*BASISV(a_d);
        // 1. Copy rhs into direction rhs
        m_rhs[offset + i] = a_rhs(ix,0); // TODO - more than 1 comp?
      }
      */
    }
  }

#if 0
  virtual void copyFABRhsDir(const int l, const FArrayBox& rhs)
  {
    if (m_verbosity >= 2)
      pout() << "  Copying FAB into RHS for direction l = " << l << endl;
    assert(rhs.box().contains(m_b));

  }

  virtual void setFABBatchRhsDir(const int l)
  {
    if (m_verbosity >= 2)
      pout() << "  Setting RHS direction l = " << l << endl;

    // Model problem to solve: diags, rhs
    int N = m_b.size(l); // size in l dir
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(l)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction l = " << l << ", size = " << N << endl;

    IntVect lsz(bsize[l % 3], bsize[(l+1) % 3], bsize[(l+2) % 3]);
    if (m_verbosity >= 2)
    {
      pout() << "  Looping over " << bloop << endl;
      pout() << "  Shifted sizes " << lsz << endl;
    }
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect liv(iv[l % 3], iv[(l+1) % 3], iv[(l+2) % 3]);
      Real val = 1+liv[1] + lsz[1]*liv[2]; 
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      if (m_verbosity >= 3)
        pout() << "    liv = " << liv << ", offset = " << offset 
          << ", val = " << val << endl;
      for (int i=0; i < N; ++i)
      {
        m_rhs[offset + i] = -val;
      }
    }
  }

  virtual void setFABBatchDiagsDir(const int l)
  {
    if (m_verbosity >= 1)
      pout() << "  Setting Diags direction l = " << l << endl;

    // Model problem to solve: diags, rhs
    int N = m_b.size(l); // size in l dir
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(l)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction l = " << l << ", size = " << N << endl;

    // Get pointers for diags, rhs
    // NOTE: these are just chunks of memory, with l-major indexing
    // --> Don't index by [iv]

    IntVect lsz(bsize[l % 3], bsize[(l+1) % 3], bsize[(l+2) % 3]);
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect liv(iv[l % 3], iv[(l+1) % 3], iv[(l+2) % 3]);
      Real val = 1+liv[1] + lsz[1]*liv[2]; 
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      if (m_verbosity >= 3)
        pout() << "    liv = " << liv << ", offset = " << offset 
          << ", val = " << val << endl;
      Real coef = val/(m_dx[l]*m_dx[l]); // diag = laplacian * pencil id
      for (int i=0; i < N; ++i)
      {
        m_sub[offset + i] = coef;
        m_diag[offset + i] = -2.0*coef;
        m_sup[offset + i] = coef;
      }
    }
  }

  virtual void solveFABBatchDir(const int l)
  {
    if (m_verbosity >= 1)
      pout() << "  Solving direction l = " << l << endl;

    int N = m_b.size(l); // size in l dir
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(l)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    IntVect lsz(bsize[l % 3], bsize[(l+1) % 3], bsize[(l+2) % 3]);

    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect liv(iv[l % 3], iv[(l+1) % 3], iv[(l+2) % 3]);
      if (m_verbosity >= 3)
        pout() << "  Solving on line " << liv << " to "
          << (liv + (N-1)*BASISV(l)) << endl;
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      double* sub = m_sub+offset;
      double* diag = m_diag+offset; 
      double* sup = m_sup+offset; 
      double* rhs = m_rhs+offset; 

      // Call lapack for Ax = b, A = diag[sub,diag,sup] 
      int nrhs = 1, ldb = N;
      int info = 0;
      if (m_verbosity >= 3)
      {
        pout() << "    Inputs to dgtsv ..." << endl;
        for (int i=0; i < N; ++i)
        { 
          pout() << "    [" << i << "], ";
          pout() << "    sub = " << sub[i] << ", ";
          pout() << "    diag = " << diag[i] << ", ";
          pout() << "    sup = " << sup[i] << ", ";
          pout() << "    rhs = " << rhs[i] << ", ";
          pout() << endl;
        }
      }
      LAPACK(GTSV,gtsv)(&N, &nrhs, sub, diag, sup, rhs, &ldb, &info);
      if (m_verbosity >= 3)
      {
        pout() << "    info = " << info << endl;
        pout() << "    soln = ";
        for (int i=0; i < N; ++i)
        {
          pout() << rhs[i] << ", ";
        }
        pout() << endl;
      }
      // TODO: copy results into original array
    }
  }

  virtual void copyFABSolnDir(const int l, FArrayBox& soln)
  {
    if (m_verbosity >= 1)
      pout() << "  Copying solution into FAB for direction l = " << l << endl;
    assert(soln.box().contains(m_b));

    // Model problem to solve: diags, rhs
    int N = m_b.size(l); // size in l dir
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(l)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction l = " << l << ", size = " << N << endl;

    IntVect lsz(bsize[l % 3], bsize[(l+1) % 3], bsize[(l+2) % 3]);
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect liv(iv[l % 3], iv[(l+1) % 3], iv[(l+2) % 3]);
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      for (int i=0; i < N; ++i)
      {
        IntVect ix = iv + i*BASISV(l);
        soln(ix,0) = m_rhs[offset + i]; // TODO - more than 1 comp?
      }
    }
  }

  virtual void doFABBatchSolve(FArrayBox& soln, FArrayBox& rhs)
  {
    CH_TIME("doFABBatchSolve");
    assert(SpaceDim == 3); // Only works for 3D

    if (m_verbosity >= 1)
      pout() << "Starting Batch solve on " << m_b << endl;

    // Grid spacing
    if (m_verbosity >= 2)
      pout() << "Grid spacing dx=" << m_dx << endl;

    // For each component of the input vector
    // for (int d=0; d < nComp; d++) // loop over soln components
    int d=0;
    {
      if (m_verbosity >= 2)
        pout() << "Component d=" << d << endl;
      // For each direction, do an ADI pass
      for (int l=0; l < SpaceDim; l++) // loop over line solve dir
      {
        // Setup rhs and diags
        setFABBatchRhsDir(l);
        setFABBatchDiagsDir(l);
       
        // Solve and copy to soln
        solveFABBatchDir(l);
      }
    }
    if (m_verbosity >= 1)
      pout() << "Batch solve complete!" << endl;
  }
#endif

private: 

// TODO change these to private eventually
public:
  Box m_b;
  RealVect m_dx;
  FArrayBox m_tmp; // for both rhs and soln, since LAPACK overwrites
  FArrayBox m_diags; // diags, a = sub/comp 0, b = diag/comp 1, c = sup/comp 2
  double *m_sub, *m_diag, *m_sup, *m_rhs;
};

#if 0
class Test2BatchSolver : public TestBatchSolver
{
public:
  Test2BatchSolver(Box a_b, RealVect a_dx) : TestBatchSolver(a_b, a_dx)
  {};
 
  virtual void setFABBatchDiagsDir(const int l)
  {
    if (m_verbosity >= 1)
      pout() << "  Setting Test2 Diags direction l = " << l << endl;

    int N = m_b.size(l); // size in l dir
    Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(l)); // 2D box
    IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction l = " << l << ", size = " << N << endl;

    IntVect lsz(bsize[l % 3], bsize[(l+1) % 3], bsize[(l+2) % 3]);
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect liv(iv[l % 3], iv[(l+1) % 3], iv[(l+2) % 3]);
      Real val = 1+liv[1] + lsz[1]*liv[2]; 
      size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      if (m_verbosity >= 3)
        pout() << "    liv = " << liv << ", offset = " << offset 
          << ", val = " << val << endl;
      Real coef = 1/(m_dx[l]*m_dx[l]); // for laplacian
      for (int i=0; i < N; ++i)
      {
        // 1D laplacian
        m_sub[offset + i] = coef;
        m_diag[offset + i] = -2.0*coef;
        m_sup[offset + i] = coef;
      }
      // Dirichlet bc's
      m_diag[offset] = -3.0*coef;
      m_diag[offset+N-1] = -3.0*coef;
    }
  }

  void setTestRHS(const int d, FArrayBox& frhs, RealVect dx)
  {
    // TODO - update for more than 1 comp?
    pout() << "Setting test RHS in a FAB" << endl;
    Box b = frhs.box();
    frhs.setVal(0);
    // rhs values for exact soln = 1 + i + 2*j + 3*k w/ Dirichlet (-3) bc's
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect i=bit();
      if (i[d] == b.smallEnd(d))
      {
        Real val = (1 + i[0] + 2*i[1] + 3*i[2]);
        frhs(i,0) += val*(-2.0)/(dx[d]*dx[d]);
      }
      if (i[d] == b.bigEnd(d))
      {
        IntVect id = i + BASISV(d);
        Real val = (1 + id[0] + 2*id[1] + 3*id[2]);
        frhs(i,0) += val*(-2.0)/(dx[d]*dx[d]);
      }
    }
  }

  void checkTestSoln(const int d, FArrayBox& fsoln, Box b, RealVect dx)
  {
    // TODO - update for more than 1 comp?
    pout() << "Checking test Soln in a FAB on " << b << endl;
    // exact soln = 1 + i + 2*j + 3*k w/ Dirichlet (-3) bc's
    // evaluated on cell centers
    Real tol = 1e-15;
    bool passed = true;
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect i=bit();
      RealVect c = RealVect(i) + BASISREALV(d)*.5; // cell center locations
      Real val = 1 + c[0] + 2*c[1] + 3*c[2];
      if (abs(fsoln(i,0) - val) > tol*val)
      {
        pout() << "  result mismatch at iv= " << i << ", exact = " << val 
          << ", soln = " << fsoln(i,0) << endl;
        passed = false;
      }
    }
    if (passed)
      pout() << "  Test passed!" << endl;
    else
      pout() << "  Test failed!" << endl;
  }

};
#endif
