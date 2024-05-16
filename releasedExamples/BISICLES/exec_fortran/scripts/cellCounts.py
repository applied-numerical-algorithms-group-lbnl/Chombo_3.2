#!/usr/bin/env python
import sys;
print sys.argv[1];
fname = sys.argv[1];
f = open(fname);
sum = 0;
n = 0;
for line in f:
#    print line;
    words = line.rsplit();
    newSum = int(words[1]);
    sum = sum + newSum;
    n = n + 1;
    print n, sum;

avgCells = sum/n;
print "Avg: ", avgCells






