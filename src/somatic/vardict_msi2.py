#!/usr/bin/env python
# File Name: vardict_msi2.py
# Author:
# Created Time: Thu 30 Aug 2018 02:59:26 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

def str_modify(ref, alt):
    #print ref,alt
    pre_count=0
    min_len=min(len(ref), len(alt))
    if min_len!=0:
        for n in range(min_len):
            if ref[n]==alt[n]:
                pre_count+=1
            else:
                break
        ref=ref[pre_count:]
        alt=alt[pre_count:]
    #print min_len,pre_count
    #print ref,alt
    ext_count=0
    min_len=min(len(ref), len(alt))
    if min_len!=0:
        for n in range(0, min_len):
            if ref[-1-n]==alt[-1-n]:
                ext_count-=1
            else:
                break
        if ext_count<0:
            ref=ref[:ext_count]
            alt=alt[:ext_count]
    #print min_len,ext_count
    #print (ref,alt,pre_count,ext_count)
    return ref,alt


def indel(Ref, Alt, lseq, rseq, non_monomer=5, monomer=3):
    ref, alt = str_modify(Ref, Alt)
    seq = lseq + Ref + rseq
    base = ""
    if len(ref) > len(alt):
        base = ref
    elif len(alt) > len(ref):
        base = alt
    subseq = list(set(seq[i:i+x] for i, j in enumerate(reversed(range(len(seq)+1))) for x in range(1, j+1)))
    if len(base) == 1:
        if base*int(non_monomer) in subseq:
            return "msi%s" % non_monomer
        else:
            return "PASS"
    elif len(base) > 1:
        if base * int(monomer) in subseq:
            return "msi%s" % monomer
        else:
            return "PASS"
    else:
        return "PASS"


def snp(Ref, Alt, lseq, rseq, non_monomer=10, monomer=5):
    seq = lseq + Ref + rseq

    subseq = list(set(seq[i:i + x] for i, j in enumerate(reversed(range(len(seq) + 1))) for x in range(1, j + 1)))
    if len(Ref) == 1:
        if Ref*int(non_monomer) in subseq or Alt*int(non_monomer) in subseq:
            return "msi%s" % non_monomer
        else:
            return "PASS"
    else:
        if Ref * int(monomer) in subseq or Alt * int(monomer) in subseq:
            return "msi%s" % monomer
        else:
            return "PASS"


def info2dic(info):
    dic = {}
    infolis = info.split(";")
    for i in infolis:
        key, val = i.split("=")
        dic[key] = val
    return dic


def format2dic(FORMAT, sampleinfo):
    dic = {}
    keylis = FORMAT.split(":")
    infolis = sampleinfo.split(":")
    for i in range(len(keylis)):
        dic[keylis[i]] = infolis[i]
    return dic


def msi(vcfile, non_monomer_indel=5, monomer_indel=3, non_monomer_snp=10, monomer_snp=5):

    with open(vcfile) as vcf:
        for i in vcf:
            i = i.rstrip()
            if i.startswith("#"):
                print i
                continue
            t = i.split()
            chrom, pos, varid, ref, alt, qual, fil, info, Format, tumor, normal = t
            infodic = info2dic(info)
            sampleinfo = format2dic(Format, tumor)
            infotag = "ALD=%s;SN=%s;NM=%s;MQ=%s;PMEAN=%s;SBF=%s;QSTD=%s"%(sampleinfo["ALD"], sampleinfo["SN"], sampleinfo["NM"], sampleinfo["MQ"], sampleinfo["PMEAN"], sampleinfo["SBF"], sampleinfo["QSTD"])

            if len(ref) == len(alt):
                tag = snp(ref, alt, infodic["LSEQ"], infodic["RSEQ"], non_monomer_snp, monomer_snp)
            else:
                tag = indel(ref, alt, infodic["LSEQ"], infodic["RSEQ"], non_monomer_indel, monomer_indel)

            print "\t".join([chrom, pos, varid, ref, alt, qual, fil,info+";MSI2=%s;%s"%(tag, infotag), Format, tumor, normal])


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print "py vacile"
        print "py vcfile non_monomer_indel monomer_indel non_monomer_snp monomer_snp"
        exit(1)
    t = sys.argv
    msi(t[1])




