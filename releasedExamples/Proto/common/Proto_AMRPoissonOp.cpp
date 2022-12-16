

#include "Proto_AMRPoissonOp.H"
#include "ProtoForAllFunctions.H"
  ///
   void
   Proto_AMRPoissonOp::
   relax(pr_lbd& a_phi, const pr_lbd& a_rhs, int a_maxiter);
   {

     double diagvalu = -2.*DIM/(a_dx*a_dx);
     double realdiag = a_alpha + a_beta*diagvalu;
     double lambda = 1./realdiag;
     for(int iter = 0; iter < a_maxiter; iter++)
     {
       for(int iredblack = 0; iredblack < 2; iredblack++)
       {
         auto & resid = const_cast<pr_lbd &>(m_resid);
         resid.setVal(0.);
         residual(resid, a_phi, a_rhs, true);
        
         for(auto dit = a_phi.begin(); dit != a_phi.end(); ++dit)
         {
           auto& resfab = resid[dit];
           auto& phifab = a_phi[dit];
           Proto::forall_i(gsrbResid, resfab.box(), phi, res, lambda, iredblack);
         } //end loop over boxes

       } //end loop over red and black

     }// end loop over iteratioons
   }


#endif
