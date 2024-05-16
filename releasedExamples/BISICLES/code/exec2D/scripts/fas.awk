#!/usr/bin/awk -f
{if ($1 == "Sum(rhs)") print " "}
{if ($2 == "converged") print " "}
{if (($1 =="AMRFAS::") && ($2 == "iteration") ) print gensub(/,/,"",1,$4) " " gensub(/,/,"",1,$8)}
