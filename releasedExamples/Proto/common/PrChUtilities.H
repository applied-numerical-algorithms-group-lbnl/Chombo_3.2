#ifndef __PR_CH_UTILIITIES_
#define __PR_CH_UTILIITIES_
#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "Proto.H"
#include "Proto_Helmholtz_Op.H"
#include "Proto_Conductivity_Op.H"
#include "Proto_Viscous_Tensor_Op.H"
#include "Proto_FastVTO.H"
#include "Proto_Resistivity_Op.H"
#include "memusage.H"
#include "Vector.H"

namespace Chombo
{
  ///
  /**
     This may look like yet another class just to keep the syntax clean and it may be that.
     This class is also a pretty wrapper around proto/include/transition/Proto_Transition_Utilities (PTU).   
     PTU on this side of the fence is a plain-old-data interface.   
     This class prettifies that somewhat and provides some code reuse for the test codes.
  **/
  template <int ncomp>
  class PrChUtilities
  {
  public:
    ///no data. all static functions
    PrChUtilities(){;}
    ~PrChUtilities(){;}

    
    typedef         LevelData<FArrayBox>           ch_ldf_cell;
    typedef         LevelData<FluxBox>             ch_ldf_flux;
    typedef         ProblemDomain                  ch_dom;
    typedef  Proto::ProblemDomain                  pr_dom;
    typedef  Proto::Point                          pr_pt;
    typedef         IntVect                        ch_iv;
    typedef  Proto::Box                            pr_box;
    typedef         Box                            ch_box;
    typedef         DisjointBoxLayout              ch_dbl;
    typedef  Proto::DisjointBoxLayout              pr_dbl;
    typedef  Proto::BoxData<double, 1    >         pr_box_data;
    typedef  Proto::LevelBoxData< double, ncomp, Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   >      pr_lbd_cell;
    typedef  Proto::LevelBoxData< double, 1,     Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_0 >      pr_lbd_xfac; 
    typedef  Proto::LevelBoxData< double, 1,     Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_1 >      pr_lbd_yfac; 
    typedef  Proto::LevelBoxData< double, 1,     Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_2 >      pr_lbd_zfac;

    typedef         AMRLevelOp<pr_lbd_cell>        ch_amrlevelop;
    typedef         AMRLevelOpFactory<pr_lbd_cell> ch_factory;
    typedef         MayDay                         ch_mayday;
    typedef         RefCountedPtr< ch_ldf_cell  >   ch_ldf_cell_ptr;
    typedef         RefCountedPtr< ch_ldf_flux  >   ch_ldf_flux_ptr;
    typedef         shared_ptr<    pr_lbd_cell  >   pr_lbd_cell_ptr;
    typedef         shared_ptr<    pr_lbd_xfac  >   pr_lbd_xfac_ptr;
    typedef         shared_ptr<    pr_lbd_yfac  >   pr_lbd_yfac_ptr;    
    typedef         shared_ptr<    pr_lbd_zfac  >   pr_lbd_zfac_ptr;    
    
    ///AMRMultigrid's interface requires raw pointers
    static Vector<ch_ldf_cell*>
    getRawVectorCh(const Vector<RefCountedPtr<ch_ldf_cell> > a_rhs_shared)
    {
      CH_TIME("PrChUtilities::getRawVectorCh");
      Vector<ch_ldf_cell *>  raw_rhs(a_rhs_shared.size(), NULL);
      for(int ilev = 0; ilev <  a_rhs_shared.size() ; ilev++)
      {
        raw_rhs[ilev] = (ch_ldf_cell*)(&(*a_rhs_shared[ilev]));
      }
      return raw_rhs;
    }
    
    ///AMRMultigrid's interface requires raw pointers
    static Vector<pr_lbd_cell*>
    getRawVectorPr(const vector<shared_ptr<pr_lbd_cell> > a_rhs_shared)
    {
      CH_TIME("PrChUtilities::getRawVectorPr");
      Vector<pr_lbd_cell *>  raw_rhs(a_rhs_shared.size(), NULL);
      for(int ilev = 0; ilev <  a_rhs_shared.size() ; ilev++)
      {
        raw_rhs[ilev] = (ch_ldf_cell*)(&(*a_rhs_shared[ilev]));
      }
      return raw_rhs;
    }
    ///
    /**
       For performance reasons, proto requires fixed-size boxes.   
       This checks that and return the stuff needed for a layout.
    **/
    static void
    getRawLayout(vector<pair<pr_pt, unsigned int> > & a_raw_layout,
                 pr_pt                              & a_box_size,
                 pr_dom                             & a_domain,
                 const ch_dbl                       & a_grids)
    {
      CH_TIME("PrChUtilities::getRawLayout");
      a_raw_layout.resize(0);
      a_domain = ProtoCh::getProtoDomain(a_grids.physDomain());
      a_box_size = pr_pt::Ones(4586);
      Vector<ch_box> boxes = a_grids.boxArray();
      Vector<int>    procs = a_grids.procIDs();
      if(boxes.size() > 0)
      {
        //proto requires fixed size grids
        pr_box  first_box_pr = ProtoCh::getProtoBox(boxes[0]);
        pr_pt box_size_first = first_box_pr.sizes();
        ///just checking the box size is constant
        for(int ibox = 1; ibox < boxes.size(); ibox++)
        {
          pr_box box_pr = ProtoCh::getProtoBox(boxes[ibox]);
          pr_pt box_size_loc  = box_pr.sizes();
          if(box_size_loc != box_size_first)
          {
            ch_mayday::Error("PrChUtilities::getRawLayout requires that boxes be fixed size");
          }
        }
        
        a_box_size = box_size_first;
        a_raw_layout.resize(boxes.size());
        for(int ibox = 0; ibox < boxes.size(); ibox++)
        {
          pr_box box_pr = ProtoCh::getProtoBox(boxes[ibox]);
          pr_box box_pt = box_pr.coarsen(a_box_size);
          if(box_pt.size() != 1)
          {
            ch_mayday::Error("PrChUtilities::getRawLayout logic error");
          }
          a_raw_layout[ibox].first  = box_pt.low();
          a_raw_layout[ibox].second = procs[ibox];
        }
      }        
    }

    ///
    static shared_ptr<pr_dbl>
    getProtoLayout(const ch_dbl& a_grids)
    {
      CH_TIME("PrChUtilities::getProtoLayout");
      vector<pair<pr_pt, unsigned int> >  raw_layout;
      pr_pt                               box_size;
      pr_dom                              domain;
      getRawLayout(raw_layout, box_size, domain, a_grids);
      shared_ptr<pr_dbl> retval =
        Proto_Transition::MetaUtilities::getLayout(raw_layout, box_size, domain);
      return retval;
      
    }
    
    ///just always tag IntVect::Zero;
    static void
    setupLLCornerGrids(Vector<ch_dbl>& a_amrGrids,
                       Vector<ch_dom>& a_amrDomains,
                       Vector<int>& a_refRatios,
                       Vector<Real>& a_amrDx,
                       int& a_max_level)
    {
      CH_TIME("PrChUtilities::setupGridsWithTagLocation");

      ParmParse ppGrids("setupLLCornerGrids");

      // get grid generation parameters
      int  maxBoxSize, blockFactor;
      Real fillRatio;
      

      ppGrids.get("max_level", a_max_level);

      ppGrids.get("max_box_size",maxBoxSize);

      ppGrids.get("block_factor", blockFactor);

      int refined_denom = 8; //1/8th of the domain gets refined
      //query because I got sick of  updating input files.
      ppGrids.query("refined_denom", refined_denom);
      ppGrids.get("fillRatio", fillRatio);
      int maxNumLevels = a_max_level +1;

      a_refRatios.resize(maxNumLevels);
      ppGrids.getarr("ref_ratio", a_refRatios, 0, maxNumLevels);

      Vector<int>  is_periodic_int;
      bool is_periodic[SpaceDim];
      ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        is_periodic[dir] = (is_periodic_int[dir] == 1);
      }

      IntVect numCells;
      Vector<int> incells(SpaceDim);
      ppGrids.getarr("num_cells", incells, 0, SpaceDim);
      numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                                 incells[3],incells[4],incells[5]) );

      RealVect domainSize = RealVect::Unit;
      if (ppGrids.contains("domain_size"))
      {
        Vector<Real> insize(SpaceDim);
        ppGrids.getarr("domain_size", insize, 0, SpaceDim);
        domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                                      insize[3],insize[4],insize[5]) );
      }

      // resize dataholders
      a_amrGrids.resize(maxNumLevels);
      a_amrDomains.resize(maxNumLevels);
      a_amrDx.resize(maxNumLevels,-1);

      // assumes dx=dy=dz
      a_amrDx[0] = domainSize[0]/numCells[0];

      IntVect domLo = IntVect::Zero;
      IntVect domHi  = numCells - IntVect::Unit;

      ch_dom baseDomain(domLo, domHi, is_periodic);
      a_amrDomains[0] = baseDomain;

      // set up refined domains, etc
      for (int lev=1; lev<= a_max_level; lev++)
      {
        a_amrDomains[lev] = a_amrDomains[lev-1];
        a_amrDomains[lev].refine(a_refRatios[lev-1]);
        a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
      }

  
      Vector<Vector<Box> > vectBoxes(a_max_level+1);
      {
        CH_TIME("setupGrids: BaseGridCreation");
        // generate base level grids

        domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

        Vector<int> procAssign(vectBoxes[0].size(), 0);

        LoadBalance(procAssign, vectBoxes[0]);

        ch_dbl baseGrids(vectBoxes[0], procAssign, baseDomain);

        a_amrGrids[0] = baseGrids;
      }


      if (a_max_level > 0)
      {
        CH_TIME("setupGrids: meshRefine and all that");
        int bufferSize = 0;
        BRMeshRefine meshGen(a_amrDomains[0],
                             a_refRatios,
                             fillRatio,
                             blockFactor,
                             bufferSize,
                             maxBoxSize);
      
        // to be used by MeshRefine...
        Vector<Vector<Box> > old_meshes(a_max_level+1);
        old_meshes[0] = vectBoxes[0];
        for (int lev=1; lev< old_meshes.size(); lev++)
        {
          Box whole_dom = a_amrDomains[lev].domainBox();
          old_meshes[lev] = Vector<Box>(1, whole_dom);
        }

        Vector<IntVectSet>  tags(maxNumLevels);
        for(int ilev = 0; ilev < maxNumLevels; ilev++)
        {
          Box tag_box(IntVect::Zero, 7*IntVect::Unit);
          tags[ilev] = IntVectSet(tag_box);
        }
        Vector<Vector<Box> > new_meshes;

        meshGen.regrid(new_meshes, tags,  0, a_max_level,  old_meshes);

        ch_dom dom_lev = baseDomain;
        for (int ilev=0; ilev <= a_max_level; ilev++)
        {
          const auto& boxes = new_meshes[ilev];
          Vector<int> procs;
          LoadBalance(procs, boxes);
          a_amrGrids[ilev] = ch_dbl(boxes, procs, dom_lev);
          if(ilev < a_max_level)
          {
            dom_lev.refine(a_refRatios[ilev]);
          }
        }
      }
    } //end setupGridsWithTagLocation
    
    ///Gets the data pointers from fabs.
    static vector<double*>
    getRawHostData(const ch_ldf_cell     & a_host_data) 
    {
      CH_TIME("PrChUtilities::getRaHostData");
      auto grids =  a_host_data.disjointBoxLayout();
      auto dit   =  grids.dataIterator();
      vector<double*> rawHostData(dit.size());
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const auto  & fab  =  a_host_data[dit[ibox]];
        rawHostData[ibox] =  (double*)(fab.dataPtr(0)); 
      }
      return rawHostData;
    }
    
    /// Copy the contents of a proto data holder into a Chombo data holder.
    static void
    copyToHost(ch_ldf_cell       & a_host_data,
               const pr_lbd_cell & a_devi_data)
    {
      CH_TIME(  "PrChUtilities::copyToHost");
      std::vector<double*> rawHostData  = getRawHostData(a_host_data);
      Proto_Transition::DataUtilities<ncomp>::copyToHost(rawHostData,
                                                         a_devi_data);
    }

    /// Copy the contents of a Chombo data holder into a proto data holder.
    static void
    copyToDevice(pr_lbd_cell             & a_devi_data,
                 const ch_ldf_cell       & a_host_data)
    {
      CH_TIME(  "PrChUtilities::copyToDevice-Cell");
      std::vector<double*> rawHostData    = getRawHostData(a_host_data);
      Proto_Transition::DataUtilities<ncomp>::copyToDevice(a_devi_data,
                                                           rawHostData);
    }

    /// Copy the contents of a Chombo flux data holder into a DIM face-based proto data holders .
    static void
    copyToDevice(pr_lbd_xfac_ptr             & a_devi_xfac_ptr,
                 pr_lbd_yfac_ptr             & a_devi_yfac_ptr,
                 pr_lbd_zfac_ptr             & a_devi_zfac_ptr,
                 const ch_ldf_flux_ptr       & a_host_flux_ptr)
    {
      CH_TIME(  "PrChUtilities::copyToDevice-Flux");

      auto grids    =  a_host_flux_ptr->disjointBoxLayout();
      auto dit_pr   =  a_devi_xfac_ptr->dataIterator();
      auto dit_ch   =  grids.dataIterator();
      
      for(int ibox  = 0; ibox < dit_ch.size(); ibox++)
      {
        auto ind_ch = dit_ch[ibox];
        auto ind_pr = dit_pr[ibox];
        {
          const auto  &   xflux_ch  = (*a_host_flux_ptr)[ind_ch][0];
          auto        &   xflux_pr  = (*a_devi_xfac_ptr)[ind_pr];
          
          const auto  &   flux_box  = xflux_ch.box();
          const double*   host_data = xflux_ch.dataPtr(0);
          double      *   devi_data = xflux_pr.data();
          
          size_t data_size = ncomp*sizeof(double)*(flux_box.numPts());
          Proto_Transition::DataUtilities<ncomp>::copyToDevice(devi_data, host_data, data_size);
        }
        {
          const auto  &   yflux_ch  = (*a_host_flux_ptr)[ind_ch][1];
          auto        &   yflux_pr  = (*a_devi_yfac_ptr)[ind_pr];
          
          const auto  &   flux_box  = yflux_ch.box();
          const double*   host_data = yflux_ch.dataPtr(0);
          double      *   devi_data = yflux_pr.data();
          
          size_t data_size = ncomp*sizeof(double)*(flux_box.numPts());
          Proto_Transition::DataUtilities<ncomp>::copyToDevice(devi_data, host_data, data_size);
        }
#if DIM==3        
        {
          const auto  &   zflux_ch  = (*a_host_flux_ptr)[ind_ch][2];
          auto        &   zflux_pr  = (*a_devi_zfac_ptr)[ind_pr];
          
          const auto  &   flux_box  = zflux_ch.box();
          const double*   host_data = zflux_ch.dataPtr(0);
          double      *   devi_data = zflux_pr.data();

          size_t data_size = ncomp*sizeof(double)*(flux_box.numPts());
          Proto_Transition::DataUtilities<ncomp>::copyToDevice(devi_data, host_data, data_size);
        }
#endif        
      }
    }
    
    ///
    static void
    setupSolver(shared_ptr< AMRMultiGrid< pr_lbd_cell > >   a_amrSolver,
                shared_ptr< LinearSolver< pr_lbd_cell > >   a_bottomSolver,
                const Vector<ch_dbl>& a_amrGrids,
                const Vector<ch_dom>& a_amrDomains,
                const Vector<int>& a_refRatios,
                const Vector<Real>& a_amrDx,
                shared_ptr<ch_factory> a_solver_factory_ptr)
    {
      CH_TIME("PrChUtilities::setupSolver");

      ParmParse ppSolver("PrChUtil_setupSolver");

      a_amrSolver->define(a_amrDomains[0], *a_solver_factory_ptr,
                          &(*a_bottomSolver), a_amrGrids.size());

      // multigrid solver parameters
      int numSmooth, numMG, maxIter;
      Real eps, hang;
      ppSolver.get("num_smooth", numSmooth);
      ppSolver.get("num_mg",     numMG);
      ppSolver.get("max_iterations", maxIter);
      ppSolver.get("tolerance", eps);
      ppSolver.get("hang",      hang);

      Real normThresh = 1.0e-30;
      a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                       numMG, maxIter, eps, hang, normThresh);
      a_amrSolver->m_verbosity = 3;

      // optional parameters
      ppSolver.query("num_pre", a_amrSolver->m_pre);
      ppSolver.query("num_post", a_amrSolver->m_post);
      ppSolver.query("num_bottom", a_amrSolver->m_bottom);
    }//end function setupSolver
  }; //end class PrChUtilities

  /// stuff for Helmholtz and conductivity
  class HelmholtzUtilities
  {
  public:

    typedef         LevelData<FArrayBox>           ch_ldf_cell;
    typedef         LevelData<FluxBox>             ch_ldf_flux;
    typedef         ProblemDomain                  ch_dom;
    typedef  Proto::ProblemDomain                  pr_dom;
    typedef  Proto::Point                          pr_pt;
    typedef         IntVect                        ch_iv;
    typedef  Proto::Box                            pr_box;
    typedef         Box                            ch_box;
    typedef         DisjointBoxLayout              ch_dbl;
    typedef  Proto::DisjointBoxLayout              pr_dbl;
    //typedef  Proto::LevelBoxData< double, ncomp >  pr_lbd_cell; 

    typedef  Proto::LevelBoxData< double, 1, Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   >  pr_lbd_cell;
    ///sometimes you need to specify
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_0 >      pr_lbd_xfac; 
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_1 >      pr_lbd_yfac; 
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_2 >      pr_lbd_zfac;

    typedef         AMRLevelOp<pr_lbd_cell>        ch_amrlevelop;
    typedef         AMRLevelOpFactory<pr_lbd_cell> ch_factory;
    typedef         MayDay                         ch_mayday;
    typedef         RefCountedPtr< ch_ldf_cell  >   ch_ldf_cell_ptr;
    typedef         RefCountedPtr< ch_ldf_flux  >   ch_ldf_flux_ptr;
    typedef         shared_ptr<    pr_lbd_cell  >   pr_lbd_cell_ptr;
    typedef         shared_ptr<    pr_lbd_xfac  >   pr_lbd_xfac_ptr;
    typedef         shared_ptr<    pr_lbd_yfac  >   pr_lbd_yfac_ptr;    
    typedef         shared_ptr<    pr_lbd_zfac  >   pr_lbd_zfac_ptr;    

    ///
    static shared_ptr<ch_factory>
    getProtoHelmholtzOpFactory(const ch_dom         &        a_coarsest_dom,
                               const Vector<int>    &        a_ref_ratios,
                               const Vector<ch_dbl> &        a_grids,
                               const double         &        a_coarsest_dx,
                               const string         &        a_domBC,
                               const double         &        a_alpha,
                               const double         &        a_beta)               
    {

      typedef PrCh_AMR_Elliptic::Proto_Helmholtz_Op_Factory prch_helm_fact;
      
      prch_helm_fact* derived_ptr =
        new prch_helm_fact(a_coarsest_dom,
                           a_grids,
                           a_ref_ratios,
                           a_coarsest_dx,        
                           a_domBC,              
                           a_alpha,              
                           a_beta);
      
      ch_factory* base_ptr =  static_cast<ch_factory*>(derived_ptr);
      return shared_ptr<ch_factory>(base_ptr);
    } //end function getProtoHelmholtzFactory

    ///
    static void
    protoizeConductivityArguments(vector<pr_dbl>                        & a_grids_pr,
                                  vector<pr_lbd_cell_ptr>               & a_acoef_pr,
                                  vector<pr_lbd_xfac_ptr>               & a_bco_x_pr,
                                  vector<pr_lbd_yfac_ptr>               & a_bco_y_pr,
                                  vector<pr_lbd_zfac_ptr>               & a_bco_z_pr,
                                  vector<int>                           & a_refrat_pr,
                                  const Vector<ch_dbl>                  & a_grids_ch,
                                  const Vector< ch_ldf_cell_ptr >       & a_acoef_ch,
                                  const Vector< ch_ldf_flux_ptr >       & a_bcoef_ch,
                                  const Vector<int>                     & a_refrat_bigV)
    {
      CH_TIME("PrChUtilities::protoizeConductivityArguments");
      a_refrat_pr.resize(      a_refrat_bigV.size());
      for(int ilev = 0; ilev < a_refrat_bigV.size(); ilev++)
      {
        a_refrat_pr[ilev] = a_refrat_bigV[ilev];
      }

      a_grids_pr.resize(        a_grids_ch.size());
      for( int ilev = 0; ilev < a_grids_ch.size(); ilev++)
      {
        shared_ptr<pr_dbl>  lev_ptr = PrChUtilities<1>::getProtoLayout(a_grids_ch[ilev]);
        a_grids_pr[ilev] = *lev_ptr;
      }

      a_acoef_pr.resize(       a_grids_pr.size());
      a_bco_x_pr.resize(       a_grids_pr.size());
      a_bco_y_pr.resize(       a_grids_pr.size());
      a_bco_z_pr.resize(       a_grids_pr.size());
      for(int ilev = 0; ilev < a_grids_pr.size(); ilev++)
      {
        {
          ch_iv ghost_iv = a_acoef_ch[ilev]->ghostVect();
          pr_pt ghost_pt = ProtoCh::getPoint(ghost_iv);
          a_acoef_pr[ilev] = pr_lbd_cell_ptr( new pr_lbd_cell(a_grids_pr[ilev], ghost_pt));
          PrChUtilities<1>::copyToDevice(*a_acoef_pr[ilev], *a_acoef_ch[ilev]);
        }
        {
          ch_iv ghost_iv = a_bcoef_ch[ilev]->ghostVect();
          pr_pt ghost_pt = ProtoCh::getPoint(ghost_iv);
          a_bco_x_pr[ilev] = pr_lbd_xfac_ptr( new pr_lbd_xfac(a_grids_pr[ilev], ghost_pt)); // 
          a_bco_y_pr[ilev] = pr_lbd_yfac_ptr( new pr_lbd_yfac(a_grids_pr[ilev], ghost_pt)); //
#if DIM==3          
          a_bco_z_pr[ilev] = pr_lbd_zfac_ptr( new pr_lbd_zfac(a_grids_pr[ilev], ghost_pt)); //
#endif          
          PrChUtilities<1>::copyToDevice(a_bco_x_pr[ilev],
                                         a_bco_y_pr[ilev],
                                         a_bco_z_pr[ilev],
                                         a_bcoef_ch[ilev]);
        }
      }
    }//end function protoizeConductivityArguments;
    ///
    static shared_ptr<ch_factory>
    getProtoConductivityFactory(const ch_dom                             & a_coarsest_dom,
                                const Vector<int>                        & a_ref_ratios,
                                const Vector<ch_dbl>                     & a_grids,
                                const double                             & a_coarsest_dx,
                                const Vector< ch_ldf_cell_ptr >          & a_acoef,
                                const Vector< ch_ldf_flux_ptr >          & a_bcoef,
                                const string                             & a_domBC,
                                const double                             & a_alpha,
                                const double                             & a_beta)
    {

      CH_TIME("PrChUtilities::getProtoConductivityFactory");
      typedef PrCh_AMR_Elliptic::Proto_Conductivity_Op_Factory  prch_cond_fact;
      vector<pr_dbl>                        grids_pr;
      vector<pr_lbd_cell_ptr>               acoef_pr;
      vector<pr_lbd_xfac_ptr>               bco_x_pr;
      vector<pr_lbd_yfac_ptr>               bco_y_pr;
      vector<pr_lbd_zfac_ptr>               bco_z_pr;
      vector<int>                           refrat_pr;

      protoizeConductivityArguments(grids_pr,
                                    acoef_pr,
                                    bco_x_pr,
                                    bco_y_pr,
                                    bco_z_pr,
                                    refrat_pr,
                                    a_grids,a_acoef,a_bcoef, a_ref_ratios);

      prch_cond_fact * derived_ptr =
        new prch_cond_fact(a_coarsest_dom,
                           grids_pr,
                           refrat_pr,
                           a_coarsest_dx,
                           a_domBC,
                           a_alpha,
                           a_beta,
                           acoef_pr,
                           bco_x_pr,
                           bco_y_pr,
                           bco_z_pr);           

      ch_factory* base_ptr = static_cast<ch_factory*>(derived_ptr);
      shared_ptr<ch_factory> retval(base_ptr);

      return retval;
    }

  }; //end class HelmholtzUtilties

  ///
  /**
     This misnamed class has utilities for the vector operators (resisitivity and viscous tensor).
  **/
  class TensorUtilities
  {
  public:

    typedef         LevelData<FArrayBox>           ch_ldf_cell;
    typedef         LevelData<FluxBox>             ch_ldf_flux;
    typedef         ProblemDomain                  ch_dom;
    typedef  Proto::ProblemDomain                  pr_dom;
    typedef  Proto::Point                          pr_pt;
    typedef         IntVect                        ch_iv;
    typedef  Proto::Box                            pr_box;
    typedef         Box                            ch_box;
    typedef         DisjointBoxLayout              ch_dbl;
    typedef  Proto::DisjointBoxLayout              pr_dbl;

    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   >  pr_lbd_cell_sca;
    typedef  Proto::LevelBoxData< double, DIM, Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   >  pr_lbd_cell_vec;
    typedef  Proto::LevelBoxData< double, 3  , Proto::MEMTYPE_DEFAULT, Proto::PR_CELL   >  pr_lbd_cell_three;
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_0 >      pr_lbd_xfac; 
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_1 >      pr_lbd_yfac; 
    typedef  Proto::LevelBoxData< double, 1,   Proto::MEMTYPE_DEFAULT, Proto::PR_FACE_2 >      pr_lbd_zfac;

    typedef         AMRLevelOp<pr_lbd_cell_vec>           ch_amrlevelop;
    typedef         AMRLevelOpFactory<pr_lbd_cell_vec>    ch_factory_vec;
    typedef         AMRLevelOpFactory<pr_lbd_cell_three>  ch_factory_three;
    typedef         MayDay                                ch_mayday;
    typedef         RefCountedPtr< ch_ldf_cell  >         ch_ldf_cell_ptr;
    typedef         RefCountedPtr< ch_ldf_flux  >         ch_ldf_flux_ptr;
    typedef         shared_ptr<    pr_lbd_cell_sca  >     pr_lbd_cell_sca_ptr;
    typedef         shared_ptr<    pr_lbd_cell_vec  >     pr_lbd_cell_vec_ptr;
    typedef         shared_ptr<    pr_lbd_xfac  >         pr_lbd_xfac_ptr;
    typedef         shared_ptr<    pr_lbd_yfac  >         pr_lbd_yfac_ptr;    
    typedef         shared_ptr<    pr_lbd_zfac  >         pr_lbd_zfac_ptr;    

    ///
    static void
    protoizeCoefficientsResist(vector<pr_dbl>                        & a_grids_pr,
                               vector<pr_lbd_xfac_ptr>               & a_eta_x_pr,
                               vector<pr_lbd_yfac_ptr>               & a_eta_y_pr,
                               vector<pr_lbd_zfac_ptr>               & a_eta_z_pr,
                               vector<int>                           & a_refrat_pr,
                               const Vector<ch_dbl>                  & a_grids_ch,
                               const Vector< ch_ldf_flux_ptr >       & a_eta_ch,
                               const Vector<int>                     & a_refrat_bigV)
    {
      CH_TIME("PrChUtilities::protoizeCoefficients");
      a_refrat_pr.resize(      a_refrat_bigV.size());
      for(int ilev = 0; ilev < a_refrat_bigV.size(); ilev++)
      {
        a_refrat_pr[ilev] = a_refrat_bigV[ilev];
      }

      a_grids_pr.resize(        a_grids_ch.size());
      for( int ilev = 0; ilev < a_grids_ch.size(); ilev++)
      {
        shared_ptr<pr_dbl>  lev_ptr = PrChUtilities<DIM>::getProtoLayout(a_grids_ch[ilev]);
        a_grids_pr[ilev] = *lev_ptr;
      }

      a_eta_x_pr.resize(       a_grids_pr.size());
      a_eta_y_pr.resize(       a_grids_pr.size());
      a_eta_z_pr.resize(       a_grids_pr.size());
      for(int ilev = 0; ilev < a_grids_pr.size(); ilev++)
      {
        ch_iv ghost_iv = a_eta_ch[ilev]->ghostVect();
        pr_pt ghost_pt = ProtoCh::getPoint(ghost_iv);
        a_eta_x_pr[ilev] = pr_lbd_xfac_ptr( new pr_lbd_xfac(a_grids_pr[ilev], ghost_pt)); // 
        a_eta_y_pr[ilev] = pr_lbd_yfac_ptr( new pr_lbd_yfac(a_grids_pr[ilev], ghost_pt)); //
#if DIM==3  
        a_eta_z_pr[ilev] = pr_lbd_zfac_ptr( new pr_lbd_zfac(a_grids_pr[ilev], ghost_pt)); //
#endif          
        PrChUtilities<1>::copyToDevice(a_eta_x_pr[ilev],
                                       a_eta_y_pr[ilev],
                                       a_eta_z_pr[ilev],
                                       a_eta_ch[ilev]);

      }
    }//end function protoizeCoefficientsResist
    ///
    static void
    protoizeCoefficientsVTO(vector<pr_dbl>                        & a_grids_pr,
                            vector<pr_lbd_cell_sca_ptr>           & a_acoef_pr,
                            vector<pr_lbd_xfac_ptr>               & a_eta_x_pr,
                            vector<pr_lbd_yfac_ptr>               & a_eta_y_pr,
                            vector<pr_lbd_zfac_ptr>               & a_eta_z_pr,
                            vector<pr_lbd_xfac_ptr>               & a_lam_x_pr,
                            vector<pr_lbd_yfac_ptr>               & a_lam_y_pr,
                            vector<pr_lbd_zfac_ptr>               & a_lam_z_pr,
                            vector<int>                           & a_refrat_pr,
                            const Vector<ch_dbl>                  & a_grids_ch,
                            const Vector< ch_ldf_cell_ptr >       & a_acoef_ch,
                            const Vector< ch_ldf_flux_ptr >       & a_eta_ch,
                            const Vector< ch_ldf_flux_ptr >       & a_lam_ch,
                            const Vector<int>                     & a_refrat_bigV)
    {
      CH_TIME("PrChUtilities::protoizeCoefficients");
      a_refrat_pr.resize(      a_refrat_bigV.size());
      for(int ilev = 0; ilev < a_refrat_bigV.size(); ilev++)
      {
        a_refrat_pr[ilev] = a_refrat_bigV[ilev];
      }

      a_grids_pr.resize(        a_grids_ch.size());
      for( int ilev = 0; ilev < a_grids_ch.size(); ilev++)
      {
        shared_ptr<pr_dbl>  lev_ptr = PrChUtilities<DIM>::getProtoLayout(a_grids_ch[ilev]);
        a_grids_pr[ilev] = *lev_ptr;
      }

      a_acoef_pr.resize(       a_grids_pr.size());
      a_eta_x_pr.resize(       a_grids_pr.size());
      a_eta_y_pr.resize(       a_grids_pr.size());
      a_eta_z_pr.resize(       a_grids_pr.size());
      a_lam_x_pr.resize(       a_grids_pr.size());
      a_lam_y_pr.resize(       a_grids_pr.size());
      a_lam_z_pr.resize(       a_grids_pr.size());
      for(int ilev = 0; ilev < a_grids_pr.size(); ilev++)
      {
        {
          ch_iv ghost_iv = a_acoef_ch[ilev]->ghostVect();
          pr_pt ghost_pt = ProtoCh::getPoint(ghost_iv);
          a_acoef_pr[ilev] = pr_lbd_cell_sca_ptr( new pr_lbd_cell_sca(a_grids_pr[ilev], ghost_pt));
          PrChUtilities<1>::copyToDevice(*a_acoef_pr[ilev], *a_acoef_ch[ilev]);
        }
        {
          ch_iv ghost_iv = a_eta_ch[ilev]->ghostVect();
          pr_pt ghost_pt = ProtoCh::getPoint(ghost_iv);
          a_eta_x_pr[ilev] = pr_lbd_xfac_ptr( new pr_lbd_xfac(a_grids_pr[ilev], ghost_pt)); // 
          a_eta_y_pr[ilev] = pr_lbd_yfac_ptr( new pr_lbd_yfac(a_grids_pr[ilev], ghost_pt)); //
#if DIM==3  
          a_eta_z_pr[ilev] = pr_lbd_zfac_ptr( new pr_lbd_zfac(a_grids_pr[ilev], ghost_pt)); //
#endif          

          a_lam_x_pr[ilev] = pr_lbd_xfac_ptr( new pr_lbd_xfac(a_grids_pr[ilev], ghost_pt)); // 
          a_lam_y_pr[ilev] = pr_lbd_yfac_ptr( new pr_lbd_yfac(a_grids_pr[ilev], ghost_pt)); //
#if DIM==3  
          a_lam_z_pr[ilev] = pr_lbd_zfac_ptr( new pr_lbd_zfac(a_grids_pr[ilev], ghost_pt)); //
#endif          
          PrChUtilities<1>::copyToDevice(a_eta_x_pr[ilev],
                                         a_eta_y_pr[ilev],
                                         a_eta_z_pr[ilev],
                                         a_eta_ch[ilev]);
          PrChUtilities<1>::copyToDevice(a_lam_x_pr[ilev],
                                         a_lam_y_pr[ilev],
                                         a_lam_z_pr[ilev],
                                         a_lam_ch[ilev]);
        }
      }
    }//end function protoizeCoefficientsVTO
    ///
    static shared_ptr<ch_factory_vec>
    getProtoViscousTensorOpFactory(const ch_dom                             & a_coarsest_dom,
                                   const Vector<int>                        & a_ref_ratios,
                                   const Vector<ch_dbl>                     & a_grids,
                                   const double                             & a_coarsest_dx,
                                   const Vector< ch_ldf_cell_ptr >          & a_aco,
                                   const Vector< ch_ldf_flux_ptr >          & a_eta,
                                   const Vector< ch_ldf_flux_ptr >          & a_lam,
                                   const string                             & a_domBC,
                                   const double                             & a_alpha,
                                   const double                             & a_beta)
    {
      CH_TIME("PrChUtilities::getProtoViscousTensorOpFactory");

      typedef PrCh_AMR_Elliptic::Proto_Viscous_Tensor_Op_Factory  prch_visc_fact;
      vector<pr_dbl>                        grids_pr;
      vector<pr_lbd_cell_sca_ptr>           acoef_pr;
      vector<pr_lbd_xfac_ptr>               eta_x_pr;
      vector<pr_lbd_yfac_ptr>               eta_y_pr;
      vector<pr_lbd_zfac_ptr>               eta_z_pr;
      vector<pr_lbd_xfac_ptr>               lam_x_pr;
      vector<pr_lbd_yfac_ptr>               lam_y_pr;
      vector<pr_lbd_zfac_ptr>               lam_z_pr;
      vector<int>                           refrat_pr;

      protoizeCoefficientsVTO(grids_pr,
                              acoef_pr,
                              eta_x_pr,
                              eta_y_pr,
                              eta_z_pr,
                              lam_x_pr,
                              lam_y_pr,
                              lam_z_pr,
                              refrat_pr,
                              a_grids,
                              a_aco,
                              a_eta,
                              a_lam,
                              a_ref_ratios);

      prch_visc_fact * derived_ptr =
        new prch_visc_fact(a_coarsest_dom,
                           grids_pr,
                           refrat_pr,
                           a_coarsest_dx,
                           a_domBC,
                           a_alpha,
                           a_beta,
                           acoef_pr,
                           eta_x_pr,
                           eta_y_pr,
                           eta_z_pr,
                           lam_x_pr,
                           lam_y_pr,
                           lam_z_pr);           

      ch_factory_vec* base_ptr = static_cast<ch_factory_vec*>(derived_ptr);
      shared_ptr<ch_factory_vec> retval(base_ptr);
      return retval;
    } //end function get viscous factory

    ///
    static shared_ptr<ch_factory_vec>
    getProtoFastVTOFactory(const ch_dom                             & a_coarsest_dom,
                           const Vector<int>                        & a_ref_ratios,
                           const Vector<ch_dbl>                     & a_grids,
                           const double                             & a_coarsest_dx,
                           const Vector< ch_ldf_cell_ptr >          & a_aco,
                           const Vector< ch_ldf_flux_ptr >          & a_eta,
                           const Vector< ch_ldf_flux_ptr >          & a_lam,
                           const string                             & a_domBC,
                           const double                             & a_alpha,
                           const double                             & a_beta)
    {
      CH_TIME("PrChUtilities::getProtoViscousTensorOpFactory");

      typedef PrCh_AMR_Elliptic::Proto_FastVTO_Factory  prch_visc_fact;
      vector<pr_dbl>                        grids_pr;
      vector<pr_lbd_cell_sca_ptr>           acoef_pr;
      vector<pr_lbd_xfac_ptr>               eta_x_pr;
      vector<pr_lbd_yfac_ptr>               eta_y_pr;
      vector<pr_lbd_zfac_ptr>               eta_z_pr;
      vector<pr_lbd_xfac_ptr>               lam_x_pr;
      vector<pr_lbd_yfac_ptr>               lam_y_pr;
      vector<pr_lbd_zfac_ptr>               lam_z_pr;
      vector<int>                           refrat_pr;

      protoizeCoefficientsVTO(grids_pr,
                              acoef_pr,
                              eta_x_pr,
                              eta_y_pr,
                              eta_z_pr,
                              lam_x_pr,
                              lam_y_pr,
                              lam_z_pr,
                              refrat_pr,
                              a_grids,
                              a_aco,
                              a_eta,
                              a_lam,
                              a_ref_ratios);

      prch_visc_fact * derived_ptr =
        new prch_visc_fact(a_coarsest_dom,
                           grids_pr,
                           refrat_pr,
                           a_coarsest_dx,
                           a_domBC,
                           a_alpha,
                           a_beta,
                           acoef_pr,
                           eta_x_pr,
                           eta_y_pr,
                           eta_z_pr,
                           lam_x_pr,
                           lam_y_pr,
                           lam_z_pr);           

      ch_factory_vec* base_ptr = static_cast<ch_factory_vec*>(derived_ptr);
      shared_ptr<ch_factory_vec> retval(base_ptr);
      return retval;
    } //end function get viscous factory
    
    ///
    static shared_ptr<ch_factory_three>
    getProtoResistivityOpFactory(const ch_dom                             & a_coarsest_dom,
                                 const Vector<int>                        & a_ref_ratios,
                                 const Vector<ch_dbl>                     & a_grids,
                                 const double                             & a_coarsest_dx,
                                 const Vector< ch_ldf_flux_ptr >          & a_eta,
                                 const string                             & a_domBC,
                                 const double                             & a_alpha,
                                 const double                             & a_beta)
    {
      CH_TIME("PrChUtilities::getProtoViscousTensorOpFactory");

      typedef PrCh_AMR_Elliptic::Proto_Resistivity_Op_Factory  prch_resist_fact;
      vector<pr_dbl>                        grids_pr;
      vector<pr_lbd_cell_sca_ptr>           acoef_pr;
      vector<pr_lbd_xfac_ptr>               eta_x_pr;
      vector<pr_lbd_yfac_ptr>               eta_y_pr;
      vector<pr_lbd_zfac_ptr>               eta_z_pr;
      vector<int>                           refrat_pr;

      protoizeCoefficientsResist(grids_pr,
                                 eta_x_pr,
                                 eta_y_pr,
                                 eta_z_pr,
                                 refrat_pr,
                                 a_grids,
                                 a_eta,
                                 a_ref_ratios);

      prch_resist_fact * derived_ptr =
        new prch_resist_fact(a_coarsest_dom,
                             grids_pr,
                             refrat_pr,
                             a_coarsest_dx,
                             a_domBC,
                             a_alpha,
                             a_beta,
                             eta_x_pr,
                             eta_y_pr,
                             eta_z_pr);

      ch_factory_three* base_ptr = static_cast<ch_factory_three*>(derived_ptr);
      shared_ptr<ch_factory_three> retval(base_ptr);
      return retval;
    }//end function getProtoResistivityFactory

  };  //end class TensorUtilities
}     //end namespace Chombo

#endif //include guard
