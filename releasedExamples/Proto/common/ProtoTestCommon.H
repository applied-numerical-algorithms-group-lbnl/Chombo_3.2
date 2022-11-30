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
#include "amr/Proto_AMRGrid.H"
/**
   Begin bits pirated from Proto/examples and Proto/tests
**/
namespace Proto
{
  Proto::AMRGrid
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
    periodicity.fill(false);
    return telescopingGrid(crseDomainSize, numLevels, refRatios, boxSizes, periodicity);
  }

  class BCHolder
  {
  public:
    BCHolder(const Box                 & a_validBox,
             const Box                 & a_domain,
             const Point               & a_srcGhost,
             string                      a_domainBCName[2*DIM])
    {
      m_validBox = a_validBox;
      m_domain   = a_domain  ;
      m_srcGhost = a_srcGhost;
      for(int ivec = 0; ivec < 2*DIM; ivec++)
      {
        m_domainBCName[ivec] = a_domainBC[ivec];
      }
    
      Point lobox = a_validBox.low();
      Point hibox = a_validBox.high();
      Point lodom = a_domain.low();
      Point hidom = a_domain.high();
      Box srcBox = a_validBox.grow(a_srcGhost);

      for(int idir = 0; idir < DIM; idir++)
      {
        int numghost = a_srcGhost[idir]; //need a bc stencil for each ghost cell
        Point loapp  = srcBox.low();
        Point hiapp  = srcBox.high();
        loapp[idir]  = lodom[idir]-1;
        hiapp[idir]  = lodom[idir]-1;

        if(lobox[idir] == lodom[idir])
        {
          int loindex = ebp_index(idir, Side::Lo);
          m_regStencilBC[loindex].resize(numghost);
          m_BCApplyBoxes[loindex].resize(numghost);
          string bcname = a_domainBCName[loindex];
          if((bcname == string("Dirichlet")) || (bcname == string("Neumann")))
          {
            double  stenSign = getStenSign(bcname);
            Box appBox0(loapp, hiapp);

            for(int igh = 0; igh < numghost; igh++)
            {
              m_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, -igh);
              m_regStencilBC[loindex][igh] = (stenSign)*Shift(Point::Basis(idir, 2*igh+1));
            }
          }
          else if(bcname == string("LinearExtrapolation"))
          {
            //just going to do one ghost cell as extrapolation.  the rest are piecewise constant
            Box appBox0(loapp, hiapp);
            m_BCApplyBoxes[loindex][0] = appBox0;
            m_regStencilBC[loindex][0] = (2.)*Shift(Point::Basis(idir, 1)) + (-1.)*Shift(Point::Basis(idir, 2));

            for(int igh = 1; igh < numghost; igh++)
            {
              m_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, -igh);
              m_regStencilBC[loindex][igh] = (1.)*Shift(Point::Basis(idir, 1));
            }
          }
          else if(bcname == string("no_bcs"))
          {
            //do nothing
          }
          else
          {
            std::cout << "bogus bc name = " << bcname << std::endl;
            PROTO_ASSERT(false, "getsecondorderbcs: bogus bc string lo");
          }
        }
        if(hibox[idir] == hidom[idir])
        {
          int hiindex = ebp_index(idir, Side::Hi);
          m_regStencilBC[hiindex].resize(numghost);
          m_BCApplyBoxes[hiindex].resize(numghost);
          string bcname = a_domainBCName[hiindex];
          Point loapp = srcBox.low();
          Point hiapp = srcBox.high();
          loapp[idir] = hidom[idir]+1;
          hiapp[idir] = hidom[idir]+1;

          if((bcname == string("Dirichlet")) || (bcname == string("Neumann")))
          {
            double  stenSign = getStenSign(bcname);
            Box appBox0(loapp, hiapp);
        
            for(int igh = 0; igh < numghost; igh++)
            {
              m_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir,  igh);
              m_regStencilBC[hiindex][igh] = (stenSign)*Shift(Point::Basis(idir, -2*igh-1));
            }
          }
          else if(bcname == string("LinearExtrapolation"))
          {
            //just going to do one ghost cell as extrapolation.  the rest are piecewise constant
            Box appBox0(loapp, hiapp);
            m_BCApplyBoxes[hiindex][0] = appBox0;
            m_regStencilBC[hiindex][0] = (2.)*Shift(Point::Basis(idir, -1)) + (-1.)*Shift(Point::Basis(idir, -2));

            for(int igh = 1; igh < numghost; igh++)
            {
              m_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir, igh);
              m_regStencilBC[hiindex][igh] = (1.)*Shift(Point::Basis(idir, -1));
            }
          }
          else if(bcname == string("no_bcs"))
          {
            //do nothing
          }
          else
          {
            std::cout << "bogus bc name = " << bcname << std::endl;
            PROTO_ASSERT(false, "getsecondorderbcs: bogus bc string hi");
          }
        }
      }
    }
      
    
    vector<Stencil<double> >    m_regStencilBC[2*DIM];
    vector<Box>                 m_BCApplyBoxes[2*DIM];       
    string                      m_domainBCName[2*DIM];               
    Box                         m_validBox;
    Box                         m_domain  ;
    Point                       m_srcGhost;

  private:
    BCHolder();
  };
  

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
      //monkey with the ghost cells
      BoxData<T, 1, MEM>& caststate = (BoxData<T, 1, MEM>&)(a_state);
      applyDomainBoundaryConditions(caststate, a_dir);
      auto GRAD = Stencil<T>::DiffCellToFace(a_dir, Side::Lo);
      //auto GRAD = 1.0*Shift::Zeros() - 1.0*Shift::Basis(a_dir, -1);
      a_flux |= GRAD(a_state, 1.0 / this->dx()[a_dir]); 
    }


    ///
    inline void applyDomainBoundaryConditions(BoxData<Real, 1>& a_src)
    {

        BoxData<data_t, 1>& castSrc = const_cast<BoxData<data_t, 1> &>(a_src.getRegData());
        for(int ibc = 0; ibc < 2*DIM; ibc++)
        {
          for(int ivec = 0; ivec < m_regStencilBCS[ibc].size(); ivec++)
          {
            Box bcsrcbox = m_BCApplyBoxes[ibc][ivec];
            m_regStencilBCS[ibc][ivec].apply(castSrc, castSrc,
                                             bcsrcbox, initToZero, 1.0);
          }
        }
    }

  private:
    
   
  };
} //end namespace proto


///
namespace Chombo
{

  ///these are  the primary tools one would use to adapt a Chombo3 code into using proto  
  /**
     Meta data for AMR calculation in Chombo Land
   */
  class AMRGrid
  {
  public:
    Vector<DisjointBoxLayout> m_grids;
    Vector<int>               m_ratio;
    unsigned int size() const
    {
      return m_grids.size();
    }
  };

  ///
  /**
     The reason we have not made this awfully convenient container before now is beyond me.
   */
  class AMRData
  {
  public:
    typedef LevelData<FArrayBox> ch_ldfab;
    
    AMRData(shared_ptr<AMRGrid> a_amrgrid, unsigned int a_nghost, unsigned int a_ncomp)
    {
      IntVect      ghostVect = a_nghost*IntVect::Unit;
      unsigned int numLevels = m_amrgrid.m_grids.size();
      m_amrgrid = a_amrgrid;
      m_data.resize(numLevels);
      m_ncomp = a_ncomp;
      for(int ilev = 0; ilev < numLevels; ilev++)
      {
        auto dbl = m_amrgrids.m_grids[ilev];
        m_data[ilev] = shared_ptr<ch_ldfab>(new ch_ldfab(dbl, a_ncomp, ghostVect));
      }
      
    }
    //beats writing access functions
    vector<shared_ptr<ch_ldfab> > m_data;
    AMRGrid                       m_grid;
    unsigned int                  m_ncomp; //strictly speaking redundant but still useful
  private:
    //not a good candidate for weak construction
    AMRData();
    void operator=(const AMRData& a_amr);
  };

  inline void
  writeDataToFile(shared_ptr<AMRData> a_data, shared_ptr<AMRGrid> a_grid, string a_filename)
  {

    Vector<string> varnames(a_data->m_ncomp);
    for(int icomp = 0; icomp < a_data->m_ncomp; icomp++)
    {
      varnames[icomp] = string("variable_") + to_string(icomp);
    }
    Real scal = 1; //cannot be zero, used for dx, time and dt
    Box domain = a_data->m_grid[0].physDomain().domainBox();
    
    writeAMRHierarchyHDF5(a_filename,
                          a_data->m_data,
                          a_data->m_grid->m_grids,
                          varnames, domain,
                          scal, scal, scal,  //dx dt time 
                          a_data->m_grid->m_ratio,
                          numlev);
  }
}

namespace PrChCommon
{
#include "PR_CH_TypeDefs.H"
  
  shared_ptr<pr_grid>  > 
  getTelescopingProtoAMRGrid()
  {
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
      retval(new  Proto::telescopingGrid(domsize, numLevels, refRatPt, boxSizePt));
    return retval;
  }

  inline void
  copyProtoDataToChombo(ch_data       & a_ch_data,
                        const pr_data & a_pr_data)
  {
    for(int ilev = 0; ilev < a_ch_data.m_grid.size(); ilev++)
    {
      auto dbl = a_ch_data->m_grid[ilev];
      auto ch_dit = dbl.dataIterator();
      auto&       ch_levdat = (*a_ch_data).m_data[ilev];
      const auto& pr_levdat = (*a_pr_data)[ilev];
      int ibox = 0;
      for(auto pr_dit: pr_levdat)
      {
        const auto& pr_fab = pr_levdat[pr_dit];
        auto      & ch_fab = ch_levdat[ch_dit[ibox]];

        double*       dataPtrHost = ch_fab.dataPtr(0);
        const double* dataPtrDevi = pr_fab.dataPtr(0);
#ifdef PROTO_CUDA
        protoMemcpy(HOST, dataPtrDevi, dataPtrHost, datSize,  protoMemcpyDeviceToHost);
        protoError err = protoGetLastError();
        if (err != protoSuccess)
        {
          fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, protoGetErrorString(err));
          std::abort();
        }
#else
        memcpy(dataPtrDevi, dataPtrHost, datSize);
#endif
        ibox++;
      }
    }
  }

  inline void
  copyChomboDataToProto(ch_data       & a_ch_data,
                        const pr_data & a_pr_data)
  {
    for(int ilev = 0; ilev < a_ch_data.m_grid.size(); ilev++)
    {
      auto dbl = a_ch_data->m_grid[ilev];
      auto ch_dit = dbl.dataIterator();
      const auto&       ch_levdat = (*a_ch_data).m_data[ilev];
      auto      & pr_levdat = (*a_pr_data)[ilev];
      int ibox = 0;
      for(auto pr_dit: pr_levdat)
      {
        auto       & pr_fab = pr_levdat[pr_dit];
        const auto & ch_fab = ch_levdat[ch_dit[ibox]];

        double*       dataPtrHost = ch_fab.dataPtr(0);
        const double* dataPtrDevi = pr_fab.dataPtr(0);
        size_t datSize = sizeof(Real)*ch_fab.box().numPts() * ch_fab.nComp();

        
#ifdef PROTO_CUDA
        protoMemcpy(DEVICE, dataPtrHost, dataPtrDevi, datSize,  protoMemcpyHostToDevice);
        protoError err = protoGetLastError();
        if (err != protoSuccess)
        {
          fprintf(stderr, "protoGetLastError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, protoGetErrorString(err));
          std::abort();
        }
#else
        memcpy(dataPtrDevi, dataPtrHost, datSize);
#endif
        ibox++;
      }
    }

  }

/**
 **/
  inline  Chombo::DisjointBoxLayout
  getChomboDBL(const Proto::DisjointBoxLayout & a_pr_grids)
  
  {
#include "PR_CH_TypeDefs.H"
    const std::vector<pair<Point, unsigned int>>  box_pairs= a_pr_grids.boxes();
    Vector<int> procs(box_pairs.size());
    Vector<int> boxes(box_pairs.size());
    for(int ibox = 0; ibox < box_pairs.size(); ibox++)
    {
      boxes[ibox] = box_pairs[ibox].first;
      procs[ibox] = box_pairs[ibox].second;
    }
    
    ch_dbl retval(boxes, procs);
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
      auto ptrefrat = a_pr_hierarchy->refRatio(ilev);
      
      (*retval).m_grids[ilev] = getChomboDBL((*a_pr_hierarchy)[ilev]);
      (*retval).m_ratio[ilev] = ptrefrat[0];
      
    }
    return retval;
  }
}
#endif
