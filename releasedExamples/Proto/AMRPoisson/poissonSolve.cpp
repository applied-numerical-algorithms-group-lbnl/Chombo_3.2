#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"

#include "CONSTANTS.H"
#include "memusage.H"
#include "Proto.H"
#include "ProtoTestCommon.H"


///this space is mostly so I can use the typedefs without putting the include in too open a scope.
namespace Introspection
{
#include "PR_CH_TypeDefs.H" 

/**
 **/
  shared_ptr<pr_amr_poisson>  >
  getAMRPoissonSolver(const shared_ptr<pr_amrgrid>& pr_hierarchy)
  {
    const auto& hierarchy = *pr_hierarchy;
    unsigned int domsize = hierarchy[0].problemDomain.sizes()[0];
    Real dx  = 1./Real(domsize);
    std::array<double, DIM> dxarray;
    dxarray.fill(dx);
    shared_ptr<pr_amr_poisson>  retval(new pr_amr_poisson(hierarchy, dxarray));
  
    return retval;
  }
/**
 **/
  void fillChomboData(shared_ptr<ch_data> a_phi,
                      shared_ptr<ch_data> a_rhs,
                      shared_ptr<ch_grid> a_grid)
  {
    int nlev = a_grid->m_grids.size();
    for(int ilev = 0; ilev < nlev; ilev++)
    {
      auto  dbl = a_grid->m_grids[ilev];
      auto  dit = dbl.dataIterator();
      auto& phi = (*(a_phi->m_data))[ilev];
      auto& rhs = (*(a_rhs->m_data))[ilev];
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        phi[dit[ibox]].setVal(0.);
        rhs[dit[ibox]].setVal(1.);
      }
    }
  }
/**
 **/
  static void runSolver()
  {
    //Proto has to drive the grid generation
    //get the grid hierarchy
    unsigned int nghost = 4;  unsigned int ncomp = 1;
    shared_ptr<pr_grid>   pr_hier =  PrChCommon::getTelescopingProtoAMRGrid();
    
    //get the chombo equivalent of the grid hierarchy
    shared_ptr<ch_grid>   ch_hier =  PrChCommon::getChomboAMRGrid(pr_grid);
    
    Point pghost = Proto::Point::Ones(nghost);
    //define the data on both sides of the divide
    shared_ptr<pr_data>   pr_rhs(new pr_data(*pr_hier, pghost,      ));
    shared_ptr<pr_data>   pr_phi(new pr_data(*pr_hier, pghost,      ));
    shared_ptr<ch_data>   ch_rhs(new ch_data(*ch_hier, nghost, ncomp));
    shared_ptr<ch_data>   ch_phi(new ch_data(*ch_hier, nghost, ncomp));
    Introspection::fillChomboData(ch_rhs, ch_phi, ch_grid);
    
    //set rhs and initial guess to phi
    PrChCommon::copyChomboDataToProto(pr_rhs, ch_rhs);
    PrChCommon::copyChomboDataToProto(pr_phi, ch_phi);

    unsigned int maxiter = 100;  Real tolerance = 1.0e-10;
    
    //get the poisson solver
    shared_ptr<pr_amr_solver> pr_pois = Introspection::getAMRPoissonSolver(pr_hier);
    pr_pois->solve(pr_phi, pr_rhs, maxiter, tolerance);
  
    TestCommon::copyProtoDataToChombo(ch_rhs, pr_rhs);
    TestCommon::copyProtoDataToChombo(ch_phi, pr_phi);

    Chombo::writeDataToFile(ch_rhs,  string("rhs.hdf5"));
    Chombo::writeDataToFile(ch_phi,  string("phi.hdf5"));
  }
}
/**
 **/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    Introspection::runSolver();
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return 0;
}
