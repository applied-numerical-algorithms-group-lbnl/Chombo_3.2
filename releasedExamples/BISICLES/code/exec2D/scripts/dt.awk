#!/usr/bin/awk -f
{if ($11 =="dt") print $2 " " ($13)}
