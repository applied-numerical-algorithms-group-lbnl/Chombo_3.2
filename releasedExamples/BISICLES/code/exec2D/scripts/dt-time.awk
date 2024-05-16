#!/usr/bin/awk -f
{if ($11 =="dt") print $7 " " ($13)}
