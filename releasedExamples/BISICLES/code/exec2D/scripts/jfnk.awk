#!/usr/bin/awk -f
{if ($1 == "Sum(rhs)") print " "}
{if ($2 == "converged") print " "}
{if ( ($1 =="Picard") && ($4 == "max(resid)")) print $3 " " ($6)}
{if (($1 =="JFNK") && ($2 == "initial")) print  "0   " ($6)}
{if (($1 =="Picard") && ($4 == "residual")) print ($3+1) " " ($7)}
{if (($1 =="JFNK") && ($4 == "residual")) print ($3+1) " " ($7)}
