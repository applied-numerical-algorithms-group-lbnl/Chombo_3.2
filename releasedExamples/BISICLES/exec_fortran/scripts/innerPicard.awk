#!/usr/bin/awk -f
{if (($1 =="JFNK") && ($2 == "initial")) print  " "}
{if ($1 =="Picard") print " "}
{if (($1 =="JFNK") && ($4 == "residual")) print " "}
{if (($1 =="BiCGStab::") && ($2 == "initial")) print "   0  "  gensub(/,/,"",1,$6)}
{if (($1 =="BiCGStab::") && ($2 == "iteration")) print "   " gensub(/,/,"",1,$4) "  "  gensub(/,/,"",1,$8)}
{if (($1 =="AMRMultiGrid::") && ($2 == "iteration")) print "   " gensub(/,/,"",1,$4) " " gensub(/,/,"",1,$8)}
{if (($2 =="GMRES") && ($3 == "residual")) print "   " gensub(/*/,"",1,gensub(/)/,"",1,$1) ) " " gensub(/,/,"",1,$5)}
{if (($1 == "0") && ($2 == "KSP")) print "  "}
{if (($2 =="KSP") ) print "   " $1 " " $5}

