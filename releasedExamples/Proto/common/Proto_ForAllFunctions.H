#ifndef FUNCTIONS_FOR_FORALL_BITS
#define FUNCTIONS_FOR_FORALL_BITS


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

#endif
