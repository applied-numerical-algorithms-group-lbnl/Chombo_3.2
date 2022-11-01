
///bit of code I can insert so I can be using the same typedefs everywhere
typedef  Chombo::Box                                           ch_box;
typedef  Chombo::FArrayBox                                     ch_fab;
typedef  Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef  Chombo::DisjointBoxLayout                             ch_dbl;
typedef  Chombo::AMRGrid                                       ch_amr_grid;
typedef  Chombo::ProblemDomain                                 ch_probdom;
                                                              
typedef  Proto::Box                                            pr_box;
typedef  Proto::AMRGrid                                        pr_amr_grid;
typedef  Proto::ProblemDomain                                  pr_probdom;
typedef  Proto::BoxOp_Laplace<double>                          pr_poissonop;
typedef  Proto::AMRSolver_FASMultigrid<pr_poissonop, double>   pr_amr_solver;


