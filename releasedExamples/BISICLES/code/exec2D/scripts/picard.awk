#!/usr/bin/awk -f
{if ( ($1 =="Picard") && ($2 == "iteration")) print $3 " " ($6)}
