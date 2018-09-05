#!/usr/bin/env python
import sys

if len(sys.argv) < 2:
    print 'py  *_rmdup.metrics'
    exit(1)

f = open(sys.argv[1]).readlines()
keys = f[7].split()
phead = """#SAMPLE,READ_PAIRS_EXAMINED,READ_PAIR_DUP,PCT_DUP,OPTICAL_DUP,OPITICAL_PCT"""
#pout  = ",".join([sys.argv[1].split('/')[-1].split('.')[0], keys[2], keys[5], keys[6], keys[7]])


SAMPLE = sys.argv[1].split('/')[-1].split('.')[0]
READ_PAIRS_EXAMINED = int(keys[2])
READ_PAIR_DUPLICATES = int(keys[5])
READ_PAIR_OPTICAL_DUPLICATES = int(keys[6])
PERCENT_DUPLICATION =  float(keys[-2]) * 100
if READ_PAIR_DUPLICATES != 0:
    OPITICAL_PCT = float(READ_PAIR_OPTICAL_DUPLICATES)/READ_PAIR_DUPLICATES * 100
else:
    OPITICAL_PCT = 0
#pout = "%s,%d,%d,%f%%,%d,%f%%" % (SAMPLE, READ_PAIRS_EXAMINED, READ_PAIR_DUPLICATES, PERCENT_DUPLICATION, READ_PAIR_OPTICAL_DUPLICATES, OPITICAL_PCT)
pout = "%s,%d,%d,%f%%" % (SAMPLE, READ_PAIRS_EXAMINED, READ_PAIR_DUPLICATES, PERCENT_DUPLICATION)









#print phead
print pout



