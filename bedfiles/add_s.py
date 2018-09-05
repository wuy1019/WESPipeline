#!/usr/bin/env python
import sys
if len(sys.argv) < 2:
	print "py bed"
	exit(1)

n = 1
for i in open(sys.argv[1]):
	i = i.rstrip()
	print "\t".join([i, "s"+str(n)])
	n += 1
