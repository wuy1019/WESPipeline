#!/usr/bin/env python
# File Name: Sporadic_Run_201_CombineInfos.py
# Author: Wuy
# Created Time: Fri 02 Sep 2016 09:44:22 AM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description: Combine Infos_TMP Infos
#########################################################################

from collections import defaultdict
from basic import *
def Ddict():
    return defaultdict(dict)

def data2dic(file, title=["SAMPLE", "READS", "LENGTH", "BASES", "GC", "Q20", "Q30", "PPM"]):
    """
    Infos
    """
    dicR1 = Ddict()
    dicR2 = Ddict()
    with open(file) as f:
        for i in f:
            i = i.rstrip()
            if i.startswith("#"):
                continue
            t = i.split(",")
            if len(t) != len(title):
                raise ValueError("title not match all infos")

            key = t[0].split("_")[0]
            pair = t[0].split(".")[0].split("_")[1]
            if pair == "R1":
                for n in range(1, len(t)):
                    dicR1[key][title[n]] = t[n]

            elif pair == "R2":
                for n in range(1, len(t)):
                    dicR2[key][title[n]] = t[n]
            else:
                raise ValueError("title not match all infos")
    return dicR1, dicR2

def file2dic(file, title=[]):
    """
    RawData Infos and CleanData Infos
    """
    dic = Ddict()
    with open(file) as f:
        for i in f:
            i = i.rstrip()
            t = i.split(",")
            if len(t) != len(title):
                raise ValueError("title not match all infos")
            key = t[0]
            for n in range(1, len(t)):
                dic[key][title[n]] = t[n]
    return dic

def comR1R2(dicR1, dicR2):
    """
    combine R1 and R2 infos
    """
    dic = Ddict()
    def ave(val1, val2):
        return round((float(val1)+float(val2))/2, 2)
    def com(val1, val2):
        return round(float(val1)+float(val2), 2)
    for key in dicR1:
        dic[key]["READS"]  = dicR1[key]["READS"]
        dic[key]["LENGTH"] = ave(dicR1[key]["LENGTH"], dicR2[key]["LENGTH"])
        dic[key]["BASES"]  = com(dicR1[key]["BASES"], dicR2[key]["BASES"])
        dic[key]["GC"]     = ave(dicR1[key]["GC"], dicR2[key]["GC"])
        dic[key]["Q20"]    = ave(dicR1[key]["Q20"], dicR2[key]["Q20"])
        dic[key]["Q30"]    = ave(dicR1[key]["Q30"], dicR2[key]["Q30"])
        dic[key]["PPM"]    = com(dicR1[key]["PPM"], dicR2[key]["PPM"])
    return dic

def combine(RawInfos, CleanInfos, CovDepOriginal, CovDepDedup,
        OnTargetraw, OnTargetadd, Duplicate, FlagStat, InsertSize, Contamination):
    RawData = comR1R2(data2dic(RawInfos)[0], data2dic(RawInfos)[1])
    CleanData = comR1R2(data2dic(CleanInfos)[0], data2dic(CleanInfos)[1])
    CovDepOri = file2dic(CovDepOriginal, title=["Sample", "MEAN_DEPTH", "1X_COVERAGE(%)",
        "10X_COVERAGE(%)", "20X_COVERAGE(%)", "50X_COVERAGE(%)", "20%MEAN_COVERAGE(%)"])
    CovDepDup = file2dic(CovDepDedup, title=["Sample", "MEAN_DEPTH_DEDUP", "1X_COVERAGE_DEDUP(%)",
        "10X_COVERAGE_DEDUP(%)", "20X_COVERAGE_DEDUP(%)", "50X_COVERAGE_DEDUP(%)",
        "20%MEAN_COVERAGE_DEDUP(%)"])
    OnTarget_raw = file2dic(OnTargetraw, title=["Sample", "CLEANBASE1",
        "MAPPEDBASE1", "ON_TARGET_CORE(%)"])
    OnTarget_add = file2dic(OnTargetadd, title=["Sample", "CLEANBASE1",
        "MAPPEDBASE2", "ON_TARGET_EXT(%)"])
    DuplicateDic = file2dic(Duplicate, title=["Sample", "TOTAL", "UNPAIRED", "DUPLICATE(%)"])
    Mapping = file2dic(FlagStat, title=["Sample", "RATIO_OF_MAPPED(%)"])
    Insert  = file2dic(InsertSize, title=["Sample", "INSERT", "INSERT_SD", "MEAN_INSERT", "MEAN_INSERT_SD"])
    DataInfo = Ddict()
    Conta = file2dic(Contamination, title=["Sample", "Contamination"])
    for key in RawData:
        DataInfo[key]["LENGTH"]             = RawData[key]["LENGTH"]
        DataInfo[key]["GC(%)"]              = RawData[key]["GC"]
        DataInfo[key]["N(ppm)"]             = RawData[key]["PPM"]
        DataInfo[key]["Q20(%)"]             = RawData[key]["Q20"]
        DataInfo[key]["Q30(%)"]             = RawData[key]["Q30"]
        DataInfo[key]["PF_READS"]           = RawData[key]["READS"]
        DataInfo[key]["CLEAN_READS"]        = CleanData[key]["READS"]
        DataInfo[key]["RATIO_OF_READS(%)"]  = round(float(DataInfo[key]["CLEAN_READS"])
                /float(DataInfo[key]["PF_READS"]), 2)
        DataInfo[key]["PF_BASES"]           = RawData[key]["BASES"]
        DataInfo[key]["CLEAN_BASES"]        = CleanData[key]["BASES"]
        DataInfo[key]["RATIO_OF_BASES(%)"]  = round(float(DataInfo[key]["CLEAN_BASES"])
                /float(DataInfo[key]["PF_BASES"]), 2)

        DataInfo[key].update(CovDepOri[key])
        DataInfo[key].update(CovDepDup[key])
        DataInfo[key].update(OnTarget_raw[key])
        DataInfo[key].update(OnTarget_add[key])
        DataInfo[key].update(DuplicateDic[key])
        DataInfo[key].update(Mapping[key])
        DataInfo[key].update(Insert[key])
        DataInfo[key].update(Conta[key])
    return DataInfo


def qcfilter(inid, qcdic, Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION, DEPTH=0):
    dic = qcdic[inid]
    lis = []
    if dic["Q30(%)"] < float(Q30):
        lis.append("Q30 Fail")
    if dic["DUPLICATE(%)"] > float(DUPLICATE):
        lis.append("Duplicate Fail")
    if float(dic["RATIO_OF_MAPPED(%)"].replace("%", "")) < float(MAPPED):
        lis.append("MAPPED Fail")
    if float(dic["1X_COVERAGE(%)"]) < float(COVERAGE):
        lis.append("Coverage Fail")
    if float(dic["INSERT"]) > float(INSERTHIGH) or float(dic["INSERT"]) < float(INSERTLOW):
        lis.append("InsertSize Fail")
    if float(dic["Contamination"]) > float(CONTAMINATION):
        lis.append("Contamination Fail")
    if float(dic["MEAN_DEPTH"]) < float(DEPTH):
        lis.append("Depth Fail")
    if not lis:
        return "PASS"
    else:
        return "|".join(lis)


def main(inid, Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION, DEPTH=0, outdir="Results"):
    qcdic = combine("{outdir}/Infos/RawInfos.csv".format(outdir=outdir),
            "{outdir}/Infos/CleanInfos.csv".format(outdir=outdir),
            "{outdir}/Infos/CovDepOriginalInfos.csv".format(outdir=outdir),
            "{outdir}/Infos/CovDepDedupInfos.csv".format(outdir=outdir),
            "{outdir}/Infos/OnTarget.raw.csv".format(outdir=outdir),
            "{outdir}/Infos/OnTarget.add.csv".format(outdir=outdir),
            "{outdir}/Infos/Duplicate.csv".format(outdir=outdir),
            "{outdir}/Infos/FlagStat.csv".format(outdir=outdir),
            "{outdir}/Infos/InsertSize.csv".format(outdir=outdir),
            "{outdir}/Infos/Contamination.csv".format(outdir=outdir) )
    QCtag = qcfilter(inid,qcdic,Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION, DEPTH)
    title = prepare(Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION)
    dic = qcdic[inid]
    dic["#SAMPLE"] = inid
    dic["QCFlag"] = QCtag
    lis = []
    for i in title:
        lis.append(dic[i])
    print ",".join(map(str, lis))


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        print "py inid Q30 DUPLICATE INSERTLOW INSERTHIGH MAPPED COVERAGE CONTAMINATION DEPTH"
        exit(1)

    t = sys.argv
    main(t[1], t[2], t[3], t[4],t[5], t[6], t[7], t[8], t[9])

