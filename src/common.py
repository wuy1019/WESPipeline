#!/usr/bin/env python
# File Name: common.py
# Author: wuy
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import ConfigParser as configparser

from commondisease.module import *


def bam2anno(inbam, targetini):

    config = readconfig(targetini)
    sample = file2sampleid(inbam)
    tmpdir = "TMP/%s" % sample
    if not checkFile(tmpdir):
        makedir(tmpdir)

    hc_cmd, hcvcf = hccmd(inbam, config.get("software", "JAVA"),
                          config.get("JavaTmp", "JAVATMP"),
                          config.get("software", "GATK"),
                          config.get("ReferenceFasta", "REFFASTA"),
                          config.get("Target", "EXT10"))

    av_cmd, avinput, avoutput = annovarcmd(hcvcf, config.get("Annovar", "CONVERT2ANNOVAR"),
                               config.get("Annovar", "TABLE_ANNOVAR"),
                               config.get("Annovar", "AnnovarDB"))
    

    pass