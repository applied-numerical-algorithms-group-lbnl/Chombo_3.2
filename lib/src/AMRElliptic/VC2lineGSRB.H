#pragma once
#include "Lapack.H"

// Special line solver for variable coef AMRPoissonOp, "VC2"
class VC2lineGSRB
{
public:
  // Note, if a_dx is not provided unit dx is used
  VC2lineGSRB(Box a_b, RealVect a_dx) : m_b(a_b), m_dx(a_dx)
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
  //         \sum_d (  b_d_i+1/2 * (phi_i+1 - phi_i)/dx_d^2
  //                 - b_d_i-1/2 * (phi_i - phi_i-1)/dx_d^2 )
  // Line relaxation in direction d to solve Lphi=b is:
  // f = \phi - \lambda (\sum_d' offdiag(L_d') \phi - b)
  // M = L_d + \sum_d' diag(L_d')
  // phi := (I + \lambda M) \ f 
  void lineRelaxRB(const int d, FArrayBox& a_phi, const FArrayBox& a_rhs,
      const Real& a_alpha, const FArrayBox& a_aCoef, 
      const Real& a_beta, const FluxBox& a_bCoef, const int whichPass)
  {
    CH_TIME("VC2lineGSRB::lineRelaxRB");
    assert(a_rhs.box().contains(m_b));

    const int N = m_b.size(d); // size in dir of solve
    const Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(d)); // 2D box
    const IntVect bsize = m_b.size(); // 3D box size for offsets
    if (m_verbosity >= 2)
      pout() << "  line direction = " << d << ", size = " << N << endl;

    // Holds which direction is 0,1,2
    const IntVect id(d % 3, (d+1) % 3, (d+2) % 3);
    const RealVect dx(m_dx[id[0]], m_dx[id[1]], m_dx[id[2]]);
    const IntVect lsz(bsize[id[0]], bsize[id[1]], bsize[id[2]]);
    if (m_verbosity >= 2)
    {
      pout() << "  Looping over " << bloop << endl;
      pout() << "  Shifted sizes " << lsz << endl;
    }
    // Basis index vectors in line directions
    const IntVect e0 = BASISV(id[0]);
    const IntVect e1 = BASISV(id[1]);
    const IntVect e2 = BASISV(id[2]);
    for (BoxIterator bit(bloop); bit.ok(); bit.next())
    {
      const IntVect iv = bit();
      if ((iv.sum() % 2) != whichPass) // Do all R then B
        continue; 
      const IntVect liv(iv[id[0]], iv[id[1]], iv[id[2]]);
      const size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]); 
      double* sub = m_sub+offset;
      double* diag = m_diag+offset; 
      double* sup = m_sup+offset; 
      double* rhs = m_rhs+offset; 
      const double tol = 1e-14;
      if (m_verbosity >= 2)
        pout() << "  line solve on " << iv << "->" << (iv+(N-1)*e0) << endl;
      if (m_verbosity >= 3)
        pout() << "Solving this system [sub, diag, sup] = rhs:" << endl;
      for (int i=0; i < N; ++i)
      {
        const IntVect ix = iv + i*e0;
        // TODO - more than 1 comp?        
        // 1. Calculate line relax rhs, copy into direction rhs
        // f = \phi - \lambda (\sum_d' offdiag(L_d') \phi - b)
        Real lres = -a_rhs(ix,0);
        lres += -a_beta*a_bCoef[id[1]](ix+e1,0)*a_phi(ix+e1,0)/(dx[1]*dx[1]);
        lres += -a_beta*a_bCoef[id[1]](ix,0)*a_phi(ix-e1,0)/(dx[1]*dx[1]);
        lres += -a_beta*a_bCoef[id[2]](ix+e2,0)*a_phi(ix+e2,0)/(dx[2]*dx[2]);
        lres += -a_beta*a_bCoef[id[2]](ix,0)*a_phi(ix-e2,0)/(dx[2]*dx[2]);
        // 1.1 What is the line relaxation coef - just the d' diag terms?
        Real Ldiagnotd = 0;
        Ldiagnotd += a_beta*a_bCoef[id[1]](ix+e1,0)/(dx[1]*dx[1]);
        Ldiagnotd += a_beta*a_bCoef[id[1]](ix,0)/(dx[1]*dx[1]);
        Ldiagnotd += a_beta*a_bCoef[id[2]](ix+e2,0)/(dx[2]*dx[2]);
        Ldiagnotd += a_beta*a_bCoef[id[2]](ix,0)/(dx[2]*dx[2]);
        Real lambda = (abs(Ldiagnotd) > tol) ? (1./4.)/Ldiagnotd : 0; // line RB weight
        // Form the line rhs
        rhs[i] = a_phi(ix,0) - lambda*lres;

        // 2. Copy system into diagonals
        // M = L_d + \sum_d' diag(L_d')
        // For M
        // diag = \sum_d ( alpha * acoef_i + beta * 
        //                 (b_d_i+1/2 + b_d_i-1/2)/dx_d^2 )
        // sub =  -beta * b_d_i-1/2 / dx_d^2
        // sup =  -beta * b_d_i+1/2 / dx_d^2
        //
        Real Msup = -a_beta*a_bCoef[id[0]](ix+e0,0)/(dx[0]*dx[0]); // i+1 coef
        Real Msub = -a_beta*a_bCoef[id[0]](ix,0)/(dx[0]*dx[0]); // i-1 coef
        Real Mdiag = a_alpha * a_aCoef(ix,0) - Msup - Msub + Ldiagnotd; // diag
        // solve A \ f = (I + \lambda M) \ f 
        sub[i] = lambda*Msub;
        diag[i] = 1 + lambda*Mdiag;
        sup[i] = lambda*Msup;
        if (m_verbosity == 3)
          pout() << "[" << sub[i] << ", " << diag[i] << ", " << sup[i] << "] = "
            <<  rhs[i] << endl;
      }

      // Solve
      int nrhs = 1, ldb = N;
      int info = 0;
      LAPACK(GTSV,gtsv)(&ldb, &nrhs, sub, diag, sup, rhs, &ldb, &info);

      // Copy updated soln back into phi
      for (int i=0; i < N; ++i)
      {
        IntVect ix = iv + i*e0;
        a_phi(ix,0) = rhs[i]; // TODO - more than 1 comp?
      }
    }
  }

private: 

// TODO change these to private eventually
public:
  Box m_b;
  RealVect m_dx;
  FArrayBox m_tmp; // for both rhs and soln, since LAPACK overwrites
  FArrayBox m_diags; // diags, a = sub/comp 0, b = diag/comp 1, c = sup/comp 2
  double *m_sub, *m_diag, *m_sup, *m_rhs;
  int m_verbosity = 0; // 0 - no output, 1 - high-level, 2 - details, 3 - lots
};

