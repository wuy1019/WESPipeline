#!/usr/bin/env python
# File Name: WEStage1.py
# Author:
# Created Time: Thu 09 Aug 2018 05:54:16 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import ConfigParser as configparser

from basic.module import *
def readconfig(inifile):
    config = configparser.ConfigParser()
    config.read(inifile)
    config.read(config.get("common", "config"))
    return config


def qc2mapping(inid, targetini, depth):
    config = readconfig(targetini)
    #mkdir dir
    prepare(config.get("QCparameter","Q30"),
             config.get("QCparameter","DUPLICATE"),
             config.get("QCparameter","INSERTLOW"),
             config.get("QCparameter","INSERTHIGH"),
             config.get("QCparameter","MAPPED"),
             config.get("QCparameter","COVERAGE"),
             config.get("QCparameter","CONTAMINATION"))
    # combine or link data
    print 'echo "==== Combine or Link Data ===="'
    combine, R1, R2, UMI = combinedata(inid)
    print combine
    print 'echo "==== fastq QC ===="'
    qc_cmd, h1, h2 = qc(inid,
            config.get("software", "JAVA"),
            config.get("software", "TRIMMOMATIC"),
            config.get("JavaTmp","JAVATMP"),
            config.get("adapters4qc","ADAPTERS4QC"),
            config.get("ScriptPath","BASIC"),
            UMIbarcode=UMI)
    print qc_cmd
    print 'echo "==== Mapping ===="'
    bwa_cmd, samfile = bwa_mem(inid,
            config.get("software", "BWA"),
            config.get("ReferenceFasta", "REFFASTA"),
            h1, h2 )
    print bwa_cmd
    print 'echo "==== Sort BAM ===="'
    sort_cmd, sortbam = sort_bam(inid,
            samfile,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "PICARD"),
            config.get("ReferenceFasta", "REFFASTA"))
    print sort_cmd

    print 'echo "==== Remove Duplicate ===="'
    rmdup_cmd, rmdupbam, metrics = rmdup(inid,
            sortbam,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "PICARD"),
            config.get("ReferenceFasta", "REFFASTA"))
    print rmdup_cmd

    print 'echo "==== RealignerTargetCreator ===="'
    realign_cmd, intervals = realign(inid,
            rmdupbam,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "GATK"),
            config.get("ReferenceFasta", "REFFASTA"),
            config.get("CorrectVCF", "INDEL_GSDB"),
            config.get("CorrectVCF", "INDEL_1000G"),
            config.get("CorrectVCF", "INDEL_MILLS"),
            config.get("Target", "EXT100"))
    print realign_cmd

    print 'echo "==== IndelRealigner ===="'
    indelrealigner_cmd, realignedbam = indelrealigner(inid,
            rmdupbam,
            intervals,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "GATK"),
            config.get("ReferenceFasta", "REFFASTA"),
            config.get("CorrectVCF", "INDEL_GSDB"),
            config.get("CorrectVCF", "INDEL_1000G"),
            config.get("CorrectVCF", "INDEL_MILLS") )
    print indelrealigner_cmd

    print 'echo "==== BaseRecalibrator ===="'
    baserecalibrator_cmd, table = baserecalibrator(inid,
            realignedbam,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "GATK"),
            config.get("ReferenceFasta", "REFFASTA"),
            config.get("CorrectVCF", "INDEL_GSDB"),
            config.get("CorrectVCF", "INDEL_1000G"),
            config.get("CorrectVCF", "INDEL_MILLS"),
            config.get("CorrectVCF", "SNP_DBSNP"),
            config.get("CorrectVCF", "SNP_GSDB"),
            config.get("Target", "EXT100"))
    print baserecalibrator_cmd

    print 'echo "==== PrintReads ===="'
    printReads_cmd, outbam = printReads(inid,
            realignedbam,
            table,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "GATK"),
            config.get("ReferenceFasta", "REFFASTA"))
    print printReads_cmd

    print 'echo "==== Infos ===="'
    info_cmds = infos(inid,
            config.get("software", "JAVA"),
            config.get("JavaTmp","JAVATMP"),
            config.get("software", "PICARD"),
            config.get("software", "SAMTOOLS"),
            config.get("software", "BEDTOOLS"),
            config.get("software", "GATK4"),
            config.get("CorrectVCF", "EXACVCF"),
            config.get("ScriptPath","BASIC"),
            sortbam,
            rmdupbam,
            config.get("Target", "RAW"),
            config.get("Target", "EXT100"))
    print info_cmds
    print 'date'

    combine_info = combineinfo(inid,
            config.get("ScriptPath","BASIC"),
            config.get("QCparameter","Q30"),
            config.get("QCparameter","DUPLICATE"),
            config.get("QCparameter","INSERTLOW"),
            config.get("QCparameter","INSERTHIGH"),
            config.get("QCparameter","MAPPED"),
            config.get("QCparameter","COVERAGE"),
            config.get("QCparameter","CONTAMINATION"),
            depth)
    print combine_info


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print "py inid panelconfig depth"
        exit(1)
    qc2mapping(sys.argv[1], sys.argv[2], sys.argv[3])
