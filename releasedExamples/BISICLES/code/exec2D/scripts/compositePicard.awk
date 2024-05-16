#!/usr/bin/awk -f
{if (($1 =="JFNK") && ($2 == "initial")) print  "0   " ($6)}
{if ($1 =="Picard") print "\n" $3 " " ($6)}
{if (($1 =="JFNK") && ($4 == "residual")) print $3 " " ($7)}
{if (($1 =="BiCGStab::") && ($2 == "iteration")) print "   " $4 $8}
{if (($1 =="AMRMultiGrid::") && ($2 == "iteration")) print "   " $4 $8}
