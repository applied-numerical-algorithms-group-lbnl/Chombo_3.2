
///bit of code I can insert so I can be using the same typedefs everywhere/
//no guards anywhere
typedef  Chombo::Box                                           ch_box;
typedef  Proto::LevelBoxData<double, 1>                        pr_lbd1;
typedef  Chombo::SideIterator                                  ch_sit; // 
typedef  Chombo::MayDay                                        ch_mayday;
typedef  Chombo::Side                                          ch_sd;
typedef  Chombo::FArrayBox                                     ch_fab;
typedef  Chombo::LevelData<ch_fab>                             ch_leveldata;
typedef  Proto::BoxData<double, 1>                             pr_bd1;
typedef  Chombo::DisjointBoxLayout                             ch_dbl;
typedef  Proto::DisjointBoxLayout                              pr_dbl;
typedef  Proto::Stencil<double>                                pr_sten;
typedef  Chombo::ProblemDomain                                 ch_probdom;
typedef  Proto::Shift                                          pr_shift;
typedef  Proto::DataIndex<Proto::BoxPartition>                 pr_di;
typedef  Proto::Point                                          pr_pt;
typedef  Proto::Box                                            pr_box;
typedef  Proto::ProblemDomain                                  pr_probdom;



