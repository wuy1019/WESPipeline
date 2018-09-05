#!/usr/bin/env python
# File Name: src/somatic/combineCNV.py
# Author:
# Created Time: Wed 22 Aug 2018 06:39:25 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

from module import *

def str2dic(s):
    dic = {}
    lis = s.split(";")
    for i in lis:
        key, val = i.split("=")
        dic[key] = val
    return dic


def combine(avinput, avoutput):
    outdic, header = txt2dic(avoutput, headidx=0, keyidx=[0, 1, 2])
    with open(avinput) as fin:
        for i in fin:
            i = i.rstrip()
            t = i.split()
            dic = str2dic(t[-1])
            key = tuple(t[:3])
            outdic[key]["CNV_Size"] = dic["Size"]
            outdic[key]["cnlr_median"] = dic["cnlr_median"]
            outdic[key]["CNVType"] = dic["CNVType"]

    header = ["Chr", "Start", "End", "Gene.refGene", "cytoBand",
              "CNV_Size", "cnlr_median", "CNVType"]
    print "\t".join(header)
    for key, val in outdic.items():
        lis = []
        for h in header:
            lis.append(val[h])
        print "\t".join(lis)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print "avinput avoutput"
        exit(1)
    combine(sys.argv[1], sys.argv[2])
