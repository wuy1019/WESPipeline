#!/usr/bin/env python
# File Name: mafilter.py
# Author:
# Created Time: Tue 04 Sep 2018 05:35:37 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

from module import *


def str2float(s):
    if s == ".":
        return 0
    elif s == "":
        return 0
    else:
        return float(s)


def mutype():
    return [
        "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
        "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
        "Splice_Site", "Translation_Start_Site", ""
    ]


def symbol():
    return ["HGNC"]


def biotype():
    return ["protein_coding"]


def ald2bias(ALD):
    rad, lad = map(float, ALD.split(","))
    return rad / (rad + lad)


def popfiltered(frequency, poplis):
    if len(filter(lambda x: x >= frequency, map(str2float, poplis))):
        return False
    else:
        return True


def filtermaf(inmaf, popfreq, tdp, tad, ndp, varqual, MSI, bias, pmean,
              totalAC, LAC, tad_hot, totalAC_hot, LAC_hot):

    dic, header = txt2dic(inmaf, headidx=0, keyidx=[5, 6, 11, 13])
    print "Commons\t" + "\t".join(header)
    for key, val in dic.items():
        commons = []

        mtype = val["Variant_Classification"]
        tdepth = int(val["t_depth"])
        talt = int(val["t_alt_count"])
        ndepth = int(val["n_depth"])
        symbol_source = val["SYMBOL_SOURCE"]
        Biotype = val["BIOTYPE"]
        poplis = [
            val["AF"], val["EAS_AF"], val["ExAC_AF"], val["ExAC_AF_EAS"],
            val["gnomAD_AF"], val["gnomAD_EAS_AF"]
        ]
        var_qual = str2float(val["variant_qual"])
        msi = float(val["MSI"])
        msi2 = val["MSI2"]
        Bias = ald2bias(val["ALD"])
        Pmean = float(val["PMEAN"])
        Cnt = str2float(val["CNT"])
        TotalAC = str2float(val["TOTALAC"])
        Lac = str2float(val["LAC"])
        NM = float(val["NM"])
        MQ = float(val["MQ"])

        # hotspot(cosmic CNT >= 20)
        if Cnt >= 20:
            hotspot = True
        else:
            hotspot = False

        #
        if mtype not in mutype():
            commons.append("Muttype")

        if tdepth < tdp:
            commons.append("tdepth")
        if ndepth < ndp:
            commons.append("ndepth")
        if symbol_source not in symbol():
            commons.append("Symbol")
        if Biotype not in biotype():
            commons.append("Biotype")
        if msi > MSI:
            commons.append("MSI")

        if not popfiltered(popfreq, poplis):
            commons.append("Popfreq")
        if var_qual < varqual:
            commons.append("varqual")
        if msi2 != "PASS":
            commons.append("MSI2")
        if Bias < bias or Bias > (1 - bias):
            commons.append("BIAS")
        if Pmean < pmean:
            commons.append("PMEAN")

        if (MQ < 55 and NM > 1) or (MQ < 60 and NM > 2):
            commons.append("MQ_NM")

        if hotspot:
            if talt < tad_hot:
                commons.append("tad")
            if TotalAC > totalAC_hot or Lac > LAC_hot:
                commons.append("Common_Variant")
        else:
            if talt < tad:
                commons.append("tad")
            if TotalAC > totalAC or Lac > LAC:
                commons.append("Common_Variant")

        if not commons:
            commons.append("PASS")

        lis = []
        lis.append(";".join(commons))
        for h in header:
            lis.append(val[h])
        print "\t".join(lis)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print "py inmaf"
        exit(1)
    filtermaf(sys.argv[1], 0.01, 20, 4, 15, 70, 8, 0.1, 25, 20, 10, 3, 50, 20)
