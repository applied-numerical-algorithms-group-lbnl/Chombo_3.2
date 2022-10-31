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
#include "BCFunc.H"
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "Proto.H"

namespace Chombo
{
  //just group stuff the same way that Proto does
  class AMRGrid
  {
  public:
    Vector<DisjointBoxLayout> m_grids;
    Vector<int>               m_refRat;
  }
}


/**
  Begin bits pirated from Proto/examples and Proto/tests
 **/

AMRGrid telescopingGrid(
        Point crseDomainSize,
        unsigned int numLevels,
        std::vector<Point> refRatios,
        std::vector<Point>   boxSizes,
        std::array<bool, DIM> periodicity)
{
    std::vector<DisjointBoxLayout> layouts;
    layouts.resize(numLevels);
    Box domainBox(crseDomainSize);
    ProblemDomain domain(crseDomainSize, periodicity);
    layouts[0].define(domain, domainBox, boxSizes[0]);
    for (int lvl = 1; lvl < numLevels; lvl++)
    {
        domain = domain.refine(refRatios[lvl-1]);
        domainBox = domainBox.grow(-domainBox.sizes()/4).refine(refRatios[lvl-1]);
        layouts[lvl].define(domain, domainBox, boxSizes[lvl]); 
    }
    return AMRGrid(layouts, refRatios, numLevels);
}

AMRGrid telescopingGrid(
        int domainSize,
        unsigned int numLevels,
        Point refRatio,
        Point boxSize)
{
    std::vector<Point> refRatios(numLevels-1, refRatio);
    std::vector<Point> boxSizes(numLevels, boxSize);
    Point crseDomainSize = Point::Ones(domainSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    return telescopingGrid(crseDomainSize, numLevels, refRatios, boxSizes, periodicity);
}

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Laplace : public BoxOp<T, 1, 0, MEM>
{
    public:

    // These functions are not optional
    inline static Point ghost() { return Point::Ones(2);}
    inline static Point auxGhost() { return Point::Zeros();}
    inline static constexpr int order() { return 4; }

    inline BoxOp_Laplace() : BoxOp<T, 1, 0, MEM>() {}
    inline BoxOp_Laplace(double a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}
    inline BoxOp_Laplace(std::array<double, DIM> a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}

    inline double spectralRadius() const {return (2.0*DIM) / pow(this->dxMin(), 2); }
    inline void flux(
            BoxData<T, 1, MEM>& a_flux,
            const BoxData<T, 1, MEM>& a_state,
            int a_dir) const
    {
        auto GRAD = Stencil<T>::DiffCellToFace(a_dir, Side::Lo);
        //auto GRAD = 1.0*Shift::Zeros() - 1.0*Shift::Basis(a_dir, -1);
        a_flux |= GRAD(a_state, 1.0 / this->dx()[a_dir]); 
    }
};

typedef Chombo::FArrayBox                                     ch_fab;
typedef Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef Chombo::DisjointBoxLayout                             ch_dbl;
typedef Chombo::Vector<ch_dbl>                                ch_hierarchy;
typedef Chombo::ProblemDomain                                 ch_probdom;
                                                              
typedef  Proto::AMRGrid                                       pr_amrgrid;
typedef  Proto::ProblemDomain                                 pr_probdom;
typedef  Proto::BoxOp_Laplace<double>                         pr_poissonop;
typedef  Proto::AMRSolver_FASMultigrid<pr_poissonop, double>  pr_amrsolver;


/**
End pirated bits.   Back to ChomboLand.
**/
shared_ptr<pr_amrgrid>  > 
getProtoAMRGrid()
{
  unsigned int domsize, boxsize, refrat, numLevels, pt;
  ParmParse pp(amrgrid);
  pp.get("domain_size"     , domsize);
  pp.get("fixed_box_size"  , boxsize);
  pp.get("refinement_ratio". refirat);

  ///Proto defines the grid hierarchy
  Point refiRatPt = Point::Ones(refirat);
  Point boxSizePt = Point::Ones(boxsize);
  shared_ptr<pr_amrgrid>
    retval(new  telescopingGrid(domsize, numLevels, refRatPt, boxSizePt));
  return retval;
}
  
shared_ptr<pr_amr_poisson>  >
getAMRPoissonSolver(const shared_ptr<pr_amrgrid>& a_amrgrid)
{
  vector<int> refrat;
  int errcode = getAMRGrids(hierarchy, refrat);
  CH_assert(errcode == 0);
  Real dx  = 1./Real(domsize);
  std::array<double, DIM> dxarray;
  dxarray.fill(dx);
  typedef BoxOp_Laplace<double> POISSON_OP;
  
  shared_ptr<AMRSolver_FASMultigrid<POISSON_OP, double> >  retval =
    new AMRSolver_FASMultigrid<POISSON_OP, double>(hierarchy, dxarray);
  
  return retval;
}

int runFASSolver()
{
  //Proto has to drive the grid generation
  shared_ptr<pr_amrgrid> = getProtoAMRGrid();
  shared_ptr<ch_amrgrid> 
  
}
/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runFASSolver();
    status += solverStatus;
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
