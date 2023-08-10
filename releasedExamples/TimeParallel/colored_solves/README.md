# coloredSolves
* NProcsTotal comes in through mpiexec.
* NProcsPerColor  comes in through ParmParse;
* NColor = NProcsTotal/NProcsPerColor
* This test runs AMRMultigrid::solve.
* First we solve NColor times serially on the world communicator.
* Next we solve  NColor times simultaneously using the split communicators.
* The timers are carefully set to make the comparison obvious.
