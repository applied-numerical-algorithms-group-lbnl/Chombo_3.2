#!/usr/bin/awk -f
{if (($1 =="JFNK") && ($2 == "initial")) print  " "}
{if ($1 =="Picard") print " "}
{if (($1 =="JFNK") && ($4 == "residual")) print " "}
{if (($1 =="BiCGStab::") && ($2 == "initial")) print "   0  "  gensub(/,/,"",1,$6)}
{if (($1 =="BiCGStab::") && ($2 == "iteration")) print "   " gensub(/,/,"",1,$4) "  "  gensub(/,/,"",1,$8)}
{if (($1 =="AMRMultiGrid::") && ($2 == "iteration")) print "   " gensub(/,/,"",1,$4) " " gensub(/,/,"",1,$8)}
{if (($1 =="RelaxSolver::") && ($2 == "initial")) print "   0  "  gensub(/,/,"",1,$6)}
{if (($1 =="RelaxSolver::") && ($2 == "iteration")) print "   " gensub(/,/,"",1,$4) "  "  gensub(/,/,"",1,$8)}
{if (($1 == "0") && ($2 == "KSP")) print "  "}
{if (($2 =="KSP") ) print "   " $2 " " $5}
{if (($1 =="KSP::") ) print "   " $4 " " $8}

