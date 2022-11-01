#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef PROTO_TEST_COMMON_H
#define PROTO_TEST_COMMON_H

#include "LevelData.H"
#include "FArrayBox.H"

#include "Vector.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"

#include "Proto.H"

/**
   Begin bits pirated from Proto/examples and Proto/tests
**/
namespace Proto
{
  inline AMRGrid
  telescopingGrid(
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


  inline AMRGrid
  telescopingGrid(
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
} //end namespace proto
/**
End pirated bits.   Back to ChomboLand.
**/


//get grids useful for testing
shared_ptr<Proto::AMRGrid>  > 
getProtoAMRGrid()
{
  #include "PR_CH_typedefs.H"
  unsigned int domsize, boxsize, refrat, maxilev, pt;
  ParmParse pp(amrgrid);
  pp.get("domain_size"     , domsize);
  pp.get("fixed_box_size"  , boxsize);
  pp.get("refinement_ratio", refirat);
  pp.get("max_level"       , maxilev);
  unsigned int numLevels = maxilev+1;
  
  ///Proto defines the grid hierarchy
  Point refiRatPt = Point::Ones(refirat);
  Point boxSizePt = Point::Ones(boxsize);
  shared_ptr<pr_amrgrid>
    retval(new  telescopingGrid(domsize, numLevels, refRatPt, boxSizePt));
  return retval;
}
namespace Chombo
{
  ///
  /**
     Meta data for AMR calculation in Chombo Land
   */
  class AMRGrid
  {
  public:
    Vector<DisjointBoxLayout> m_grids;
    Vector<int>               m_ratio;
  };

  ///
  /**
     Meta data and data all in one place.
     Why we never made this awfully convenient container is beyond me.
   */
  class AMRHierarchy
  {
  public:

  private:
    //not a good candidate for weak construction
    AMRHierarchy();
    void operator=(const AMRHierarchy & a_amr);


    Vector<LevelData<FArrayBox>* > m_data;
    AMRGrid                        m_hierarchy;

    
  };
}
/**
 **/
inline  Chombo::DisjointBoxLayout
  getChomboDBL(const Proto::DisjointBoxLayout & a_pr_grids)
  
{
#include "PR_CH_TYPEDEFS.H"
  Vector<
  DisjointBoxLayout retval(boxes, loads);
  return retval;
}
/**
 **/
shared_ptr<ch_amr_grid>
  getChomboAMRGrid(const shared_ptr<Proto::AMRGrid>& a_pr_hierarchy)
{

  shared_ptr<ch_amr_grid> retval(newChombo::AMRGrid());
  unsigned int numlev = a_pr_hierarchy-numLevels();
  retval->m_grids.resize(numlev);
  retval->m_ratio.resize(numlev, 2);
  for(int ilev = 0; ilev < numlev; ilev++)
  {
    auto ptrefrat = a_pr_hierarchy->refRatio(ilev)
      
    (*retval).m_grids= getChomboDBL((*a_pr_hierarchy)[ilev]);
    (*retval).m_ratio= ptrefrat[0]
    
  }
  return retval;
}
#endif
