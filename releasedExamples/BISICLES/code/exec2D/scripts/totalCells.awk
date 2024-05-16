#!/usr/bin/awk -f
{if ( ($4 =="cells") && ($5 == "advanced") ) print ($3) " " ($7)}

