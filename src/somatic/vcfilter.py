#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# File Name: src/somatic/vcfilter.py
# Author:
# Created Time: Mon 20 Aug 2018 08:54:04 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import sys

#未添加假突变数据库过滤
def mutectv1filter(invcf, outvcf):
    out = open(outvcf, "w")
    with open(invcf) as vcfin:
        for i in vcfin:
            i = i.rstrip()
            if i.startswith("#"):
                out.write(i + "\n")
                continue
            t = i.split()
            fit = t[6]
            if fit != "PASS":
                continue

            out.write(i+"\n")
    out.close()
    return outvcf

#未添加假突变数据库过滤
def scalplefilter(invcf, outvcf, tumorid, normalid):
    out = open(outvcf, "w")
    with open(invcf) as vcfin:
        for i in vcfin:
            i = i.rstrip()
            if i.startswith("##"):
                out.write(i + "\n")

            elif i.startswith("#"):
                t = i.split()
                header = '\t'.join(t[:9] + ['TUMOR', 'CONTROL'])
                out.write(header+"\n")
                if t[9] == tumorid:
                    TUMOR = 9
                    CONTROL = 10
                elif t[9] == normalid:
                    CONTROL = 9
                    TUMOR = 10
                else:
                    sys.stderr.write('ERROR:\n%s\t%s\tNot in VCF\n' % (tumorid, normalid))
                    exit(1)

            else:
                t = i.rstrip().split("\t")
                if t[6] == "PASS" and t[7].startswith("SOMATIC"):
                    out.write("\t".join(t[:9] + [t[TUMOR], t[CONTROL]])+"\n")
    out.close()
    return outvcf


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "py caller xxx"
        print "caller: mutectv1 scalpel"
        exit(1)
    t = sys.argv
    caller = sys.argv[1]

    if caller == "mutectv1":
        mutectv1filter(t[2], t[3])
    elif caller == "scalpel":
        scalplefilter(t[2], t[3], t[4], t[5])

