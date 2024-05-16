#!/usr/bin/awk -f
{if ($1 =="Picard") print "" }
{if (($1 =="JFNK") && ($4 == "residual")) print "\n"}
{if (($1 =="AMRMultiGrid::") && ($9 == "rate")) { print "   " (gensub(/,/,"",1,$4)-1) "  " $11} }
