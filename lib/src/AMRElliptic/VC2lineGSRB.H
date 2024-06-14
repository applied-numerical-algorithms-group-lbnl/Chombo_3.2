#pragma once
#include "Lapack.H"
#include <cassert>
#include "CH_Timer.H"
#include "NamespaceHeader.H"
// Special line solver for variable coef AMRPoissonOp, "VC2"
class VC2lineGSRB
{
public:
  // Note, if a_dx is not provided unit dx is used
  VC2lineGSRB(Box a_b, RealVect a_dx) : m_b(a_b), m_dx(a_dx)
  {
    CH_TIME("VC2lineGSRB constructor");
    assert(SpaceDim == 3);   // Only works for 3D

    /*
    // TODO - generalize this, dim, # comp?
    m_diags.define(m_b,3);   // diags, sub/comp 0, diag/comp 1, sup/comp 2
    m_tmp.define(m_b,1);     // tmp for the LAPACK rhs's --> soln's

    // Get pointers for diags, rhs
    // NOTE: these are just chunks of memory, with l-major indexing
    // --> Don't index by [iv]
    m_sub  = m_diags.dataPtr(0);
    m_diag = m_diags.dataPtr(1);
    m_sup  = m_diags.dataPtr(2);
    m_rhs  = m_tmp.dataPtr(0);
    */
    const int Nmax = a_b.longside();
    m_sub = new double[Nmax];
    m_diag = new double[Nmax];
    m_sup = new double[Nmax];
    m_rhs = new double[Nmax];
  }

  ~VC2lineGSRB()
  {
    if (m_sub != NULL)
      delete m_sub;
    if (m_diag != NULL)
      delete m_diag;
    if (m_sup != NULL)
      delete m_sup;
    if (m_rhs != NULL)
      delete m_rhs;
  }

  // Does a line relaxation for this op in direction d:
  // L phi = alpha * acoef_i * phi_i - beta *
  //         \sum_d (  b_d_i+1/2 * (phi_i+1 - phi_i)/dx_d^2
  //                 - b_d_i-1/2 * (phi_i - phi_i-1)/dx_d^2 )
  // Line relaxation in direction d to solve Lphi=b is:
  // f = b + beta
  //     \sum_d' (b_d'_i+1/2 \phi^k_i+1 + b_d'_i-1/2 \phi^k_i-1)/dx_d'^2
  // M = L_d + \sum_d' diag(L_d')
  // phi := M \ f
  void lineRelaxRB(const int d, FArrayBox& a_phi, const FArrayBox& a_rhs,
                   FArrayBox& a_res, const Real& a_alpha, const FArrayBox& a_aCoef,
                   const Real& a_beta, const FluxBox& a_bCoef, const int whichPass)
  {
#if CH_SPACEDIM == 3
    CH_TIMERS("VC2lineGSRB::lineRelaxRB");
    CH_TIMER("VC2lineGSRB::lineRelaxRB::setDiagsRHS", t2);
    CH_TIMER("VC2lineGSRB::lineRelaxRB::LAPACK", t3);
    CH_TIMER("VC2lineGSRB::lineRelaxRB::copySoln", t4);
    CH_TIMER("VC2lineGSRB::lineRelaxRB::resid", t5);
    assert(a_rhs.box().contains(m_b));

    const int N = m_b.size(d); // size in dir of solve
    const Box bloop(m_b.smallEnd(), m_b.bigEnd()-(N-1)*BASISV(d)); // 2D box
    const IntVect bsize = m_b.size(); // 3D box size for offsets

    if (m_verbosity >= 2) {
        pout() << "  line direction = " << d << ", size = " << N << endl;
    }

    // Holds which direction is 0,1,2
    const IntVect id(d % 3, (d+1) % 3, (d+2) % 3);
    const RealVect dx(m_dx[id[0]], m_dx[id[1]], m_dx[id[2]]);
    const IntVect lsz(bsize[id[0]], bsize[id[1]], bsize[id[2]]);

    if (m_verbosity >= 2) {
        pout() << "  Looping over " << bloop << endl;
        pout() << "  Shifted sizes " << lsz << endl;
    }

    // Basis index vectors in line directions
    const IntVect e0 = BASISV(id[0]);
    const IntVect e1 = BASISV(id[1]);
    const IntVect e2 = BASISV(id[2]);
    for (BoxIterator bit(bloop); (whichPass < 2) && bit.ok(); bit.next()) {
        const IntVect iv = bit();
        if ((iv.sum() % 2) != whichPass) { // Do all R then B
            continue;
        }
        const IntVect liv(iv[id[0]], iv[id[1]], iv[id[2]]);
        // const size_t offset = liv[0] + lsz[0]*(liv[1] + lsz[1]*liv[2]);
        const size_t offset = 0;
        double* sub  = m_sub + offset;
        double* diag = m_diag + offset;
        double* sup  = m_sup + offset;
        double* rhs  = m_rhs + offset;

        if (m_verbosity >= 2) {
            pout() << "  line solve on " << iv << "->" << (iv+(N-1)*e0) << endl;
        }
        if (m_verbosity >= 3) {
            pout() << "Solving this system [sub, diag, sup] = rhs:" << endl;
        }

        CH_START(t2);
        for (int i=0; i < N; ++i) {
            const IntVect ix = iv + i*e0;
            // TODO - more than 1 comp?
            // 1. Calculate line relax rhs, copy into direction rhs
            // f = b + beta
            //     \sum_d' (b_d'_i+1/2 \phi^k_i+1 + b_d'_i-1/2 \phi^k_i-1)/dx_d'^2
            Real lres = a_rhs(ix,0);
            lres += a_beta*a_bCoef[id[1]](ix+e1,0) * a_phi(ix+e1,0)/(dx[1]*dx[1]);
            lres += a_beta*a_bCoef[id[1]](ix,0)    * a_phi(ix-e1,0)/(dx[1]*dx[1]);
            lres += a_beta*a_bCoef[id[2]](ix+e2,0) * a_phi(ix+e2,0)/(dx[2]*dx[2]);
            lres += a_beta*a_bCoef[id[2]](ix,0)    * a_phi(ix-e2,0)/(dx[2]*dx[2]);
            // 1.1 Correct for ghost cells in this direction
            if (i==0) {
                lres += a_beta*a_bCoef[id[0]](ix,0)*a_phi(ix-e0,0)/(dx[0]*dx[0]);
            }
            if (i==N-1) {
                lres += a_beta*a_bCoef[id[0]](ix+e0,0)*a_phi(ix+e0,0)/(dx[0]*dx[0]);
            }
            rhs[i] = lres;

            // 2. Copy system into diagonals
            // M = L_d + \sum_d' diag(L_d')
            // For M
            // diag = alpha * acoef_i
            //        + \sum_d&d' ( beta * (b_d_i+1/2 + b_d_i-1/2)/dx_d^2 )
            // sub =  -beta * b_d_i-1/2 / dx_d^2
            // sup =  -beta * b_d_i+1/2 / dx_d^2
            //
            Real Ldiagnotd = 0;
            Ldiagnotd += a_beta*a_bCoef[id[1]](ix+e1,0)/(dx[1]*dx[1]);
            Ldiagnotd += a_beta*a_bCoef[id[1]](ix,0)/(dx[1]*dx[1]);
            Ldiagnotd += a_beta*a_bCoef[id[2]](ix+e2,0)/(dx[2]*dx[2]);
            Ldiagnotd += a_beta*a_bCoef[id[2]](ix,0)/(dx[2]*dx[2]);
            Real Msup = -a_beta*a_bCoef[id[0]](ix+e0,0)/(dx[0]*dx[0]); // i+1 coef
            Real Msub = -a_beta*a_bCoef[id[0]](ix,0)/(dx[0]*dx[0]); // i-1 coef
            Real Mdiag = a_alpha * a_aCoef(ix,0) - Msup - Msub + Ldiagnotd; // diag
            if (i>0) { // sub index is shifted to start at 0
                sub[i-1] = Msub;
            }
            diag[i] = Mdiag;
            sup[i] = Msup;

            if (m_verbosity == 3) {
                if (i>0) { // sub index is shifted to start at 0
                    pout() << "[" << sub[i-1] << ", ";
                } else {
                    pout() << "[ - , ";
                }
                pout()  << diag[i] << ", ";
                if (i<N-1) {
                    pout() << sup[i] << "]";
                } else {
                    pout() << " - ]";
                }
                pout() << " = " << rhs[i] << endl;
            }

        }
        CH_STOP(t2);

        CH_START(t3);
        // Solve
        int nrhs = 1, ldb = N;
        int info = 0;
        LAPACK(GTSV,gtsv)(&ldb, &nrhs, sub, diag, sup, rhs, &ldb, &info);
        CH_STOP(t3);

        CH_START(t4);
        // Copy updated soln back into phi
        for (int i=0; i < N; ++i) {
            IntVect ix = iv + i*e0;
            a_phi(ix,0) = rhs[i]; // TODO - more than 1 comp?
            if (m_verbosity == 3) {
                pout() << "soln[" << i << "] = " << rhs[i] << endl;
            }
        }
        CH_STOP(t4);
    }

    /*
    CH_START(t5);
    IntVect b[SpaceDim];
    b[0] = BASISV(0);
    b[1] = BASISV(1);
    b[2] = BASISV(2);
    for (BoxIterator bit(m_b); bit.ok(); bit.next()) {
        const IntVect iv = bit();
        a_res(iv,0) = a_alpha*a_aCoef(iv,0)*a_phi(iv,0) - a_rhs(iv,0);
        for (int l=0; l < SpaceDim; ++l) {
            a_res(iv,0) -= a_beta*a_bCoef[l](iv+b[l],0)
                           *(a_phi(iv+b[l],0) - a_phi(iv,0)) / (dx[l]*dx[l]);
            a_res(iv,0) += a_beta*a_bCoef[l](iv,0)
                           *(a_phi(iv,0) - a_phi(iv-b[l],0)) / (dx[l]*dx[l]);
        }
    }
    CH_STOP(t5);
    */
#endif
  }

  /* Is this used ? */
  void relaxGSRB(FArrayBox& a_phi, const FArrayBox& a_rhs,
                 FArrayBox& a_res, const Real& a_alpha, const FArrayBox& a_aCoef,
                 const Real& a_beta, const FluxBox& a_bCoef, const int whichPass)
  {
#if CH_SPACEDIM == 3
      CH_TIME("VC2lineGSRB::relaxGSRB");
      assert(a_rhs.box().contains(m_b));

      // Calculate the residual
      IntVect b[SpaceDim];
      b[0] = BASISV(0);
      b[1] = BASISV(1);
      b[2] = BASISV(2);
      for (BoxIterator bit(m_b); bit.ok(); bit.next()) {
          const IntVect iv = bit();
          if ((iv.sum() % 2) != whichPass) { // Do all R then B
              continue;
          }
          a_res(iv,0) = a_alpha*a_aCoef(iv,0)*a_phi(iv,0) - a_rhs(iv,0);
          Real rcoef = a_alpha*a_aCoef(iv,0);
          for (int d=0; d < SpaceDim; ++d) {
              a_res(iv,0) -= a_beta*a_bCoef[d](iv+b[d],0)
                             *(a_phi(iv+b[d],0) - a_phi(iv,0)) / (m_dx[d]*m_dx[d]);
              a_res(iv,0) += a_beta*a_bCoef[d](iv,0)
                             *(a_phi(iv,0) - a_phi(iv-b[d],0)) / (m_dx[d]*m_dx[d]);
              rcoef += a_beta*(a_bCoef[d](iv,0) + a_bCoef[d](iv+b[d],0))
                / (m_dx[d]*m_dx[d]);
          }
          // Update
          a_phi(iv,0) -= a_res(iv,0) / rcoef;
      }
#endif
  }

private:

// TODO change these to private eventually
public:
  Box m_b;
  RealVect m_dx;
  FArrayBox m_tmp;   // for both rhs and soln, since LAPACK overwrites
  FArrayBox m_diags; // diags, a = sub/comp 0, b = diag/comp 1, c = sup/comp 2
  double *m_sub, *m_diag, *m_sup, *m_rhs;
  int m_verbosity = 0; // 0 - no output, 1 - high-level, 2 - details, 3 - lots
};

#include "NamespaceFooter.H"
