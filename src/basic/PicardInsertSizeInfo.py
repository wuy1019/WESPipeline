#!/usr/bin/env python
import sys

if len(sys.argv) < 2:
    print 'py  picard-insertsize.out '
    exit(1)

f = open(sys.argv[1]).readlines()
keys = f[7].split()
phead = """#SAMPLE,MEDIAN_INSERT_SIZE,MEDIAN_ABSOLUTE_DEVIATION,MEAN_INSERT_SIZE,STANDARD_DEVIATION"""
pout  = ",".join([sys.argv[1].split('/')[-1].split('.')[0], keys[0], keys[1], keys[4], keys[5]])
#print phead
print pout



