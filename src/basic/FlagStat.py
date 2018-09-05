#!/usr/bin/env python
import sys
if len(sys.argv) < 2:
	print 'py in.flagstat'
	exit(1)

pct = open(sys.argv[1]).readlines()[4].split('(')[-1].split(':')[0]
print "%s,%s" % (sys.argv[1].split('/')[-1].split('.')[0], pct)

