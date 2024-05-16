#!/usr/bin/awk -f
{if (($2 =="SNES") ) print "   "}
{if (($2 =="KSP") ) print "   " $1 " " $5}

