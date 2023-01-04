

#include "PR_CH_AMRPoissonOp.H"

/****/
/****/
PROTO_KERNEL_START 
unsigned int
gsrbResidF(int                       a_pt[DIM],
           Proto::Var<double, 1>     a_phi,
           Proto::Var<double, 1>     a_res,
           double     a_lambda
           int     a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    double diagvalu = -2.*DIM/(a_dx*a_dx);
    double realdiag = a_alpha + a_beta*diagvalu;

    double lambda = 1./realdiag;

    double phival = a_phi(0);
    double resval = a_res(0);

    a_phi(0) = phival + lambda*resval;
  }
  return 0;
}
PROTO_KERNEL_END(gsrbResidF, gsrbResid)

namespace Chombo
{
  ///
  void
  PR_CH_AMRPoissonOp::
  relax(pr_lbd& a_phi, const pr_lbd& a_rhs, int a_maxiter)
  {

    double lambda = m_lambda;
    for(int iter = 0; iter < a_maxiter; iter++)
    {
      for(int iredblack = 0; iredblack < 2; iredblack++)
      {
        m_resid.setVal(0.);
        residual(m_resid, a_phi, a_rhs, true);
        
        for(auto dit = a_phi.begin(); dit != a_phi.end(); ++dit)
        {
          auto& resfab = m_resid[dit];
          auto& phifab = a_phi[dit];
          Proto::forall_i(gsrbResid, resfab.box(), phi, res, lambda, iredblack);
        } //end loop over boxes

      } //end loop over red and black

    }// end loop over iteratioons
  }
  ///
  void
  PR_CH_AMRPoissonOp::
  restrictResidual(pr_lbd& a_resCoar, const pr_lbd& a_phiFine, const pr_lbd & a_rhsFine)
  {

    auto &   residFine = m_resid;
    residual(residFine, a_phiFine, a_rhsFine, true);
    auto dbl = a_phiFine.layout();
    //this one is all on proc
    for(auto dit = a_phi.begin(); dit != a_phi.end(); ++dit)
    {
      auto valid = dbl[dit];
      auto& resfabFine = residFine[dit];
      auto& resfabCoar = a_resCoar[dit];
      resfabCoar |= m_restrict(resfabFine, valid);
    } //end loop over boxes
  }
  
  ///
  void
  PR_CH_AMRPoissonOp::
  prolongIncrement(pr_lbd& a_phiFine, const pr_lbd& a_corCoar)
  {
    auto dbl = a_phiFine.layout();
    //this one is all on proc
    for(auto dit = a_phi.begin(); dit != a_phi.end(); ++dit)
    {
      auto valid = dbl[dit];
      auto& phifabFine = a_phiFine[dit];
      auto& corfabCoar = a_corCoar[dit];
      phifabFine |= m_prolong[icolor](corfabCoar, valid);
    } //end loop over boxes
  }
    
}

