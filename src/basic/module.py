#!/usr/bin/env python
# File Name: basic.py
# Author:
# Created Time: Thu 09 Aug 2018 02:54:42 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

import os

from ctDNAv3_adapter_cutter_cln2 import *


def prepare(Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION, outdir = "Results"):
    info_annotation = """##Infos\n##Q30 threshold:%s\n##DUPLICATE threshold:%s\n##INSERT Size threshold:%s-%s\
            \n##MAPPING RATIO threshold:%s\n##COVERAGE threshold:%s\n##CONTAMINATION threshold:%s\
            \n"""%(Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION)


    infofile = "%s/Infos.csv"%outdir
    title = ["#SAMPLE", "QCFlag", "LENGTH", "GC(%)", "N(ppm)", "Q20(%)", "Q30(%)", "PF_READS", "CLEAN_READS",
         "RATIO_OF_READS(%)", "PF_BASES", "CLEAN_BASES", "RATIO_OF_BASES(%)", "INSERT",
         "DUPLICATE(%)", "ON_TARGET_CORE(%)", "ON_TARGET_EXT(%)", "RATIO_OF_MAPPED(%)",
         "MEAN_DEPTH", "1X_COVERAGE(%)", "10X_COVERAGE(%)", "20X_COVERAGE(%)", "50X_COVERAGE(%)",
         "20%MEAN_COVERAGE(%)", "MEAN_DEPTH_DEDUP", "1X_COVERAGE_DEDUP(%)", "10X_COVERAGE_DEDUP(%)",
         "20X_COVERAGE_DEDUP(%)", "50X_COVERAGE_DEDUP(%)", "20%MEAN_COVERAGE_DEDUP(%)", "Contamination"]
    if os.path.exists(infofile):
        return title
    else:
        os.system("mkdir -p TMP raw RawData %s/CleanData %s/BAM %s/Infos/InsertSize"%(outdir, outdir, outdir))
        w = open(infofile, "w")
        w.write(info_annotation)
        w.write(",".join(title) + "\n")
        w.close()


def combinedata(inid, rawdir = "raw", enddir = "RawData"):
    rawlis = os.popen("ls %s/*gz"%rawdir).read().strip().split("\n")
    lis = []
    for s in rawlis:
        ID = s.split("_")[0].split("/")[-1]
        if inid == ID and "_R1_" in s:
            lis.append(s)
    cmd1 = ""
    cmd2 = ""
    R1 = "%s/%s_R1.fastq.gz"%(enddir, inid)
    R2 = "%s/%s_R2.fastq.gz"%(enddir, inid)
    UMIbarcode=False
    for s in lis:
        if index_check(s) == True:
            UMIbarcode=True
            break

    if not lis:
        sys.stderr.write("%s data is not found"% inid)

    elif len(lis) >= 2:
        cmd1 = "cat %s > %s/%s_R1.fastq.gz" % (
                " ".join(lis), enddir, inid)
        cmd2 = cmd1.replace("_R1", "_R2")
    else:
        cmd1 = "ln -s  ../%s %s/%s_R1.fastq.gz" % (lis[0], enddir, inid)
        cmd2 = cmd1.replace("_R1", "_R2")
    cmds = "\n".join([cmd1, cmd2])
    return cmds, R1, R2, UMIbarcode


def qc(inid, java, trimmomatic, javatmp_dir, adapters4qc, scriptspath, UMIbarcode=False, outdir="Results"):
    r1 = 'RawData/%s_R1.fastq.gz' % inid
    r2 = 'RawData/%s_R2.fastq.gz' % inid
    h1 = '%s/CleanData/%s_R1.fastq.gz'%(outdir, inid)
    h2 = '%s/CleanData/%s_R2.fastq.gz'%(outdir, inid)
    tmp_r1 = 'TMP/%s_R1_trimm.fastq.gz' %inid
    tmp_r2 = 'TMP/%s_R2_trimm.fastq.gz' %inid
    qc_raw = '%s/FastqInfos.py %s >> %s/Infos/RawInfos.csv &&\
            %s/FastqInfos.py %s >> %s/Infos/RawInfos.csv &\
            '%(scriptspath, r1, outdir, scriptspath, r2, outdir)


    if UMIbarcode == False:
        qc_cmd = "%s -d64 -Djava.io.tmpdir=%s -Xmx30g -jar %s PE -phred33 -threads  8 %s %s %s \
                TMP/%s_R1.unpaired.fastq.gz %s TMP/%s_R2.unpaired.fastq.gz ILLUMINACLIP:%s:2:20:10:1:true \
                LEADING:15 TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:36"%(java, javatmp_dir, trimmomatic,
                        r1, r2, h1, inid, h2, inid, adapters4qc)
    else:
        qc_cmd = "%s -d64 -Djava.io.tmpdir=%s -Xmx30g -jar %s PE -phred33 -threads  8 %s %s %s \
                TMP/%s_R1.unpaired.fastq.gz %s TMP/%s_R2.unpaired.fastq.gz ILLUMINACLIP:%s:2:20:10:1:true \
                TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:36 && \
                %s/ctDNAv3_adapter_cutter_cln2.py %s %s %s %s" %(java, javatmp_dir, trimmomatic, r1, r2, tmp_r1,
                        inid, tmp_r2, inid, adapters4qc, scriptspath, tmp_r1, tmp_r2, h1, h2)


    qc_hq = '%s/FastqInfos.py %s >> %s/Infos/CleanInfos.csv &&\
            %s/FastqInfos.py %s >> %s/Infos/CleanInfos.csv &\
            '%(scriptspath, h1, outdir, scriptspath, h2, outdir)

    cmds = "\n".join([qc_raw, qc_cmd, qc_hq])

    return cmds, h1, h2


def bwa_mem(inid, bwa, ref, r1, r2, outdir="Results"):
    samfile = 'TMP/%s.sam.gz'%inid
    cmd = '%s mem -M -t 8 -R "@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:illumina" \
            %s  %s %s | gzip > %s'%(bwa, inid, inid, inid, ref, r1, r2, samfile)
    return cmd, samfile


def sort_bam(inid, samfile, java, javatmp_dir, picard, ref ):
    sortbam = 'TMP/%s.sorted.bam'%inid
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s SortSam \
            INPUT=%s \
            OUTPUT=%s \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true'%(java, javatmp_dir, picard, samfile, sortbam)
    return cmd, sortbam

def rmdup(inid, inbam, java, javatmp_dir, picard, ref ):
    rmdupbam = 'TMP/%s.sorted.rmdup.bam'%inid
    metrics = 'TMP/%s.sorted.rmdup.metrics'%inid
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s MarkDuplicates \
            INPUT=%s \
            OUTPUT=%s \
            METRICS_FILE=%s \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            CREATE_INDEX=true'%(java, javatmp_dir, picard, inbam, rmdupbam, metrics)
    return cmd, rmdupbam, metrics


def realign(inid, inbam, java, javatmp_dir, gatk, ref,
        indel_gsdb, indel_1000g, indel_mills, target):
    intervals = 'TMP/%s.sorted.rmdup.realigner.intervals'%inid
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s \
            -T RealignerTargetCreator \
            -R %s \
            -I %s \
            -known %s \
            -known %s \
            -known %s \
            -L %s \
            -o %s \
            --disable_auto_index_creation_and_locking_when_reading_rods \
            '%(java, javatmp_dir, gatk, ref, inbam, indel_gsdb,
                    indel_1000g, indel_mills, target, intervals)
    return cmd, intervals


def indelrealigner(inid, inbam, intervals, java, javatmp_dir, gatk,
        ref, indel_gsdb, indel_1000g, indel_mills):
    realignedbam = 'TMP/%s.sorted.rmdup.realigned.bam'%inid
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s \
            -T IndelRealigner \
            -R %s \
            -I %s \
            -known %s \
            -known %s \
            -known %s \
            -targetIntervals %s \
            -o %s \
            --disable_auto_index_creation_and_locking_when_reading_rods \
            '%(java, javatmp_dir, gatk, ref, inbam, indel_gsdb, indel_1000g,
                    indel_mills, intervals, realignedbam)
    return cmd, realignedbam


def baserecalibrator(inid, inbam, java, javatmp_dir, gatk, ref,
        indel_gsdb, indel_1000g, indel_mills,snp_gsdb, snp_dbsnp, target):
    table = 'TMP/%s.sorted.rmdup.recal_data.table'%inid
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s \
            -T BaseRecalibrator \
            -R %s \
            -I %s \
            -L %s \
            -nct 8 \
            -knownSites %s \
            -knownSites %s \
            -knownSites %s \
            -knownSites %s \
            -knownSites %s \
            -o %s \
            --disable_auto_index_creation_and_locking_when_reading_rods \
            '%(java, javatmp_dir, gatk, ref, inbam, target, indel_gsdb,
                    indel_1000g, indel_mills,snp_gsdb, snp_dbsnp, table)
    return cmd, table


def printReads(inid, inbam, table, java, javatmp_dir, picard, ref, outdir="Results"):
    outbam = '%s/BAM/%s.sorted.rmdup.realigned.recal.bam'%(outdir, inid)
    cmd = '%s -Djava.io.tmpdir=%s \
            -jar %s \
            -T PrintReads \
            -R %s \
            -nct 8 \
            -I %s \
            -BQSR %s \
            -o %s \
            --disable_auto_index_creation_and_locking_when_reading_rods \
            '%(java, javatmp_dir, picard, ref, inbam, table, outbam)
    return cmd, outbam


def run_contamination(inid, inbam, gatk4, exacVCF, outdir="Results"):
    output = 'TMP/%s.CalculateContamination.txt'%(inid)
    cmd1  = '%s GetPileupSummaries \
            -I %s \
            -V %s \
            -O TMP/%s.GetPileupSummaries.txt \
            '%(gatk4, inbam, exacVCF, inid)
    cmd2 = '%s CalculateContamination \
            -I TMP/%s.GetPileupSummaries.txt \
            -O %s \
            '%(gatk4, inid, output)
    cmds = "&&".join([cmd1, cmd2])
    return cmds, output


def infos(inid, java, javatmp_dir, picard, samtools, bedtools, gatk4, exacVCF,
        scriptspath, sortbam, rmdupbam, rawbed, ext100bed, outdir="Results"):
    inf1 = '%s coverage \
            -d -abam %s \
            -b %s \
            > TMP/%s.sorted.raw.base'%(bedtools, sortbam, rawbed,inid)
    inf2 = '%s/InfoExtracter.py \
            TMP/%s.sorted.raw.base \
            %s/Infos/CleanInfos.csv'%(scriptspath, inid, outdir)
    inf3 = "%s coverage \
            -d -abam %s \
            -b %s > TMP/%s.dedup.raw.base \
            " % (bedtools, rmdupbam, rawbed, inid)
    inf4 = "%s/InfoExtracter.py \
            TMP/%s.dedup.raw.base %s/Infos/CleanInfos.csv \
            " % (scriptspath, inid, outdir)
    inf5 = "%s coverage \
            -d -abam %s \
            -b %s > TMP/%s.sorted.add.base \
            " %(bedtools, sortbam, ext100bed, inid)
    inf6 = "%s/InfoExtracter.py \
            TMP/%s.sorted.add.base %s/Infos/CleanInfos.csv \
            " % (scriptspath, inid, outdir)
    inf7 = "%s flagstat \
            %s \
            > TMP/%s.flagstat \
            " % (samtools, sortbam, inid)
    inf8 = "%s/FlagStat.py \
            TMP/%s.flagstat \
            >> %s/Infos/FlagStat.csv \
            " % (scriptspath , inid, outdir)
    inf9 = "%s/Duplicate.py \
            TMP/%s.sorted.rmdup.metrics \
            >> %s/Infos/Duplicate.csv" % (scriptspath, inid, outdir)
    inf10 = "%s -Djava.io.tmpdir=%s \
            -jar %s CollectInsertSizeMetrics \
            I=%s \
            O=TMP/%s.insertsize.txt \
            HISTOGRAM_FILE=%s/Infos/InsertSize/%s.pdf \
            " % (java, javatmp_dir, picard, sortbam, inid, outdir, inid)
    inf11 = "%s/PicardInsertSizeInfo.py \
            TMP/%s.insertsize.txt \
            >> %s/Infos/InsertSize.csv" % (scriptspath, inid, outdir)
    inf12, ctout = run_contamination(inid,
            rmdupbam,
            gatk4,
            exacVCF)
    inf13 = "%s/Contamination.py %s \
            >> %s/Infos/Contamination.csv \
            "%(scriptspath, ctout, outdir)

    cmds = " && ".join([inf1, inf2, inf3, inf4, inf5, inf6, inf7, inf8, inf9, inf10, inf11, inf12, inf13])
    return cmds


def combineinfo(inid, scriptspath, Q30, DUPLICATE, INSERTLOW, INSERTHIGH, MAPPED, COVERAGE, CONTAMINATION, DEPTH, outdir="Results"):
    cmd = "%s/CombineInfos.py %s \
            %s  %s %s %s %s %s %s %s \
            >> %s/Infos.csv"%(scriptspath, inid,Q30, DUPLICATE,
                    INSERTLOW, INSERTHIGH, MAPPED,
                    COVERAGE, CONTAMINATION, DEPTH, outdir)
    return cmd
