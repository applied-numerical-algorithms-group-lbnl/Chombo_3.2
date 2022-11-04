
#ifndef PR_CH_TYPEDEFS__H
#define PR_CH_TYPEDEFS__H

///bit of code I can insert so I can be using the same typedefs everywhere
typedef  Chombo::Box                                           ch_box;
typedef  Chombo::FArrayBox                                     ch_fab;
typedef  Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef  Chombo::DisjointBoxLayout                             ch_dbl;
typedef  Chombo::AMRGrid                                       ch_amr_grid;
typedef  Chombo::ProblemDomain                                 ch_probdom;
                                                              
typedef  Proto::Box                                            pr_box;
typedef  Proto::AMRGrid                                        pr_grid;
typedef  Chombo::AMRGrid                                       ch_grid;
typedef  Proto::AMRData                                        pr_data;
typedef  Chombo::AMRData                                       ch_data;
typedef  Proto::ProblemDomain                                  pr_probdom;
typedef  Proto::BoxOp_Laplace<double>                          pr_poissonop;
typedef  Proto::AMRSolver_FASMultigrid<pr_poissonop, double>   pr_amr_solver;


#endif
