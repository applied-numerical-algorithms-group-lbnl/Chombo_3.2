

#ifndef __PR_CH_TYPEDEFS__
#define __PR_CH_TYPEDEFS__
///bit of code I can insert so I can be using the same typedefs everywhere
typedef  Chombo::Box                                           ch_box;
typedef  Chombo::FArrayBox                                     ch_fab;
typedef  Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef  Proto::LevelBoxData<Real, 1>                          pr_lbd;
typedef  Chombo::DisjointBoxLayout                             ch_dbl;
typedef  Proto::DisjointBoxLayout                              pr_dbl;
typedef  Proto::Stencil<double>                                pr_sten;
typedef  Chombo::ProblemDomain                                 ch_probdom;
                                                              
typedef  Proto::Point                                          pr_pt;
typedef  Proto::Box                                            pr_box;
typedef  Proto::ProblemDomain                                  pr_probdom;

#endif

