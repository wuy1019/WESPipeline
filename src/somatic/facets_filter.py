#!/usr/bin/env python
# File Name: src/somatic/facets_filter.py
# Author:
# Created Time: Mon 20 Aug 2018 09:00:28 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import sys
import numpy as np


def facetsfilter(inseg, reffai, gainlog2, losslog2, outseg, avinput):
    faidic = {}
    with open(reffai) as fai:
        for i in fai:
            t = i.rstrip().split()
            faidic[t[0]] = t[1]
    faidic["23"] = faidic["X"]
    faidic["24"] = faidic["Y"]

    segfile = open(inseg).readlines()
    out = open(outseg, "w")
    avin = open(avinput, "w")
    Median = np.median(map(float, list(np.loadtxt(inseg, str, delimiter=',')[:, 4][1:])))
    out.write(segfile[0] )
    for i in segfile[1:]:
        i = i.rstrip()
        t = i.split(",")
        log2 = float(t[4]) - Median
        start, end = map(int, map(float, [t[9], t[10]]))
        if end > int(faidic[t[0]]):
            end = int(faidic[t[0]])
        size = end - start
        out.write(",".join(map(str, t[:4] + [log2] + t[5:9] + [start, end] + t[11:])) + "\n")
        if t[0] == "23":
            chrom = "X"
        elif t[0] == "24":
            chrom = "Y"
        else:
            chrom = t[0]
        #loss

        if log2 <= float(losslog2):
            avin.write("\t".join([chrom, t[9], t[10], "0", "0",
                                  "CopyNumber=%s;Size=%s;cnlr_median=%s;CNVType=loss"%(t[12], size, log2)]) + "\n")

        elif log2 >= float(gainlog2):
            avin.write("\t".join([chrom, t[9], t[10], "0", "0",
                                  "CopyNumber=%s;Size=%s;cnlr_median=%s;CNVType=gain" % (t[12], size, log2)])+"\n")
    avin.close()
    out.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "py inseg xxx"
        exit(1)
    t = sys.argv
    facetsfilter(t[1], t[2], t[3], t[4], t[5], t[6])
