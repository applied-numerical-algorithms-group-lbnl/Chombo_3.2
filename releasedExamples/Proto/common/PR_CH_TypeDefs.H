

///bit of code I can insert so I can be using the same typedefs everywhere
typedef  Chombo::Box                                           ch_box;
typedef  Chombo::FArrayBox                                     ch_fab;
typedef  Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef  Chombo::DisjointBoxLayout                             ch_dbl;
typedef  Protoo::DisjointBoxLayout                             pr_dbl;
typedef  Chombo::AMRGrid                                       ch_amr_grid;
typedef  Chombo::ProblemDomain                                 ch_probdom;
                                                              
typedef  Proto::Point                                          pr_pt;
typedef  Proto::Box                                            pr_box;
typedef  Proto::AMRGrid                                        pr_grid;
typedef  Chombo::AMRGrid                                       ch_grid;
typedef  Proto::AMRData                                        pr_data;
typedef  Chombo::AMRData                                       ch_data;
typedef  Chombo::AMRData                                       ch_data;
typedef  Proto::ProblemDomain                                  pr_probdom;


typedef  Proto::LevelBoxData<Real>                             pr_leveldata;
typedef Chombo::AMRMultigrid<pr_leveldata>                     ch_solver;
typedef Chombo::PR_CH_AMRPoissonOpFactory                      prch_amrp_factory;

