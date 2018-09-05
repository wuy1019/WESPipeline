#!/usr/bin/env python
# File Name: scripts/comtamination.py
# Author:
# Created Time: Wed 15 Aug 2018 01:48:49 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

def get_comtamination(txt):
    sampleid = txt.split("/")[-1].split(".")[0]
    contamination = open(txt).readlines()[-1].split()[1]
    print ",".join([sampleid, contamination])




if __name__ == "__main__":
    import sys
    get_comtamination(sys.argv[1])
