#!/usr/bin/env python
# File Name: vcfaddtag.py
# Author:
# Created Time: Thu 23 Aug 2018 10:07:06 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import gzip
import re
import sys
import argparse

from module import *


def vcf2dic(vcfile):
    """

    :param vcfile:file name, str
    :return: defaultdict
    """
    dic = Ddict()
    headerdic= {}
    if vcfile.endswith(".gz"):
        fin = gzip.open(vcfile, "r")
    else:
        fin = open(vcfile, "r")

    for i in fin:
        i = i.rstrip()
        if i.startswith("#"):
            try:
                infokey = re.search(r"##INFO=<ID=(.+?),.+", i).group(1)
                headerdic[infokey] = i
            except:
                pass
            continue

        chrom, pos, ID, ref, alt, qual, Filter, info = i.split()
        key = (chrom, pos, ref, alt)
        for it in info.split(";"):
            if "=" in it:
                key2, val = re.split("=", it, 1)
                dic[key][key2] = val

    return dic, headerdic


def fileopen(filepath):
    if filepath == "-":
        fin = sys.stdin.readlines()
    elif filepath.endswith(".gz"):
        fin = gzip.open(filepath).readlines()
    else:
        fin = open(filepath).readlines()
    return fin

def header(openfile, headerdic, infolis):
    fin = openfile
    head1 = []
    head2 = []
    for i in fin:
        i = i.rstrip()
        if i.startswith("##"):
            head1.append(i)
        elif i.startswith("#"):
            head2.append(i)
        else:
            break
    print "\n".join(head1)
    for itag in infolis:
        if headerdic.has_key(itag):
            print headerdic[itag]
    print "\n".join(head2)



def addtag(invcf, tagvcf, taglis, null="."):
    """

    :param invcf: in vcf file,support stand in or gzip file.
    :param tagvcf: vcf file(database)
    :param tag: add tagvcf INFOs, list
    :param null: null value string, default "."
    :return: print vcf
    """
    tagdic, headerdic = vcf2dic(tagvcf)

    fin = fileopen(invcf)
    header(fin, headerdic, taglis)
#    fin.seek(0)
    for i in fin:
        i = i.rstrip()
        if i.startswith("#"):
            continue

        t = i.split()
        chrom, pos, ID, ref, alt, qual, Filter, info = t[:8]
        key = (chrom, pos, ref, alt)
        if tagdic.has_key(key):
            for itag in taglis:
                if not tagdic[key].has_key(itag):
                    sys.stderr.write("Error tag: %s tag not exist in %s\n"%(itag, tagvcf))
                    exit(1)
                info = info + ";%s=%s"%(itag, tagdic[key][itag])
        else:
            for itag in taglis:
                info = info+";%s=%s"%(itag,null)

        print "\t".join(t[:7] + [info] + t[8:])


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="input vcf file", required=True)
    parser.add_argument("-d", "--database", help="database vcf file, add tag from this file", required=True)
    parser.add_argument("-t", "--taglist", help="Comma-delimited names of INFO fields to retain as columns in database vcf", required=True)
    parser.add_argument("-n", "--nastring", help="string to display when a score is not available. By default . string is printed in the output file.", default=".")
    return parser.parse_args()


if __name__ == "__main__":
    if len(sys.argv)==1:
        sys.argv.append("-h")
    args=parseargs()
    addtag(args.input_file, args.database, args.taglist.split(","), args.nastring)



