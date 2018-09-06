#!/usr/bin/env python
# File Name: common.py
# Author: wuy
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import ConfigParser as configparser

import combineanno
from module import *


def printcmd(cmd, echostr):
    print 'date'
    print 'echo "==== Start Run %s ===="' % echostr
    print cmd


def bam2anno(inbam, targetini):

    config = readconfig(targetini)
    sample = file2sampleid(inbam)
    tmpdir = "TMP/%s" % sample
    outdir = "Results/ComDis/raw"

    makedir(tmpdir, outdir)

    hc_cmd, hcvcf = hccmd(inbam, config.get("software", "JAVA"),
                          config.get("JavaTmp", "JAVATMP"),
                          config.get("software", "GATK"),
                          config.get("ReferenceFasta", "REFFASTA"),
                          config.get("Target", "EXT10"))
    printcmd(hc_cmd, "HaplotypeCaller ")

    av_cmd, avinput, avoutput = annovarcmd(
        hcvcf, config.get("Annovar", "CONVERT2ANNOVAR"),
        config.get("Annovar", "TABLE_ANNOVAR"),
        config.get("Annovar", "AnnovarDB"))
    printcmd(av_cmd, "AnnoVar")
    inter_cmd, interout = intervar(hcvcf, config.get("InterVar", "InterVar"))
    printcmd(inter_cmd, "InterVar")
    vep_cmd, vepcsv = vepcmd(hcvcf, config.get("vep83", "vep"),
                             config.get("ReferenceFasta", "REFFASTA"),
                             config.get("ScriptPath", "COMMONDISEASE"))
    printcmd(vep_cmd, "VEP83")
    #combine vcf, vep, intervar, annovar
    anno1 = "%s/%s.anno1.tsv" % (outdir, sample)
    combine_cmd, combine_tsv = combineanno(
        config.get("ScriptPath", "COMMONDISEASE"), hcvcf, avoutput,
        vepcsv, interout, config.get("genetic", "omim"),
        config.get("genetic", "hgmd"), anno1)
    printcmd(combine_cmd, "Combine anno files")
