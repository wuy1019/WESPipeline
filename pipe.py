#!/usr/bin/env python
# File Name: src/pipe.py
# Author: wuy
# Created Time: Thu 16 Aug 2018 10:45:34 AM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import os

from module import *


def mutectv1snv(tumorbam, normalbam, inifile, scriptout, outdir="Results/mutectv1"):
    config = readconfig(inifile)
    script = open(scriptout, "w")
    tumorid, normalid = map(bam2id, [tumorbam, normalbam])
    tmpdir = "TMP/%s/mutectv1" % (tumorid)
    #call
    cmd1, rawvcf = mutectv1(tumorbam, normalbam,
                            config.get("software","JAVA7"),
                            config.get("JavaTmp","JAVATMP"),
                            config.get("software", "MUTECT"),
                            config.get("ReferenceFasta", "REFFASTA"),
                            config.get("Parameters2MuTectv1", "fraction_contamination"),
                            config.get("Target","EXT10"),
                            config.get("Parameters2MuTectv1","min_qscore"),
                            config.get("AnnoVCF", "Cosmic"),
                            config.get("AnnoVCF", "dbSNP"),
                            tmpdir)

    #filter
    filtervcf = "%s/%s-%s.vcf"%(outdir, tumorid, normalid)
    cmd2 = "%s/vcfilter.py mutectv1 %s %s"%(config.get("ScriptPath", "SOMATIC"),
                                            rawvcf, filtervcf)

    #anno
    maf = "%s/%s-%s.maf"%(outdir, tumorid, normalid)
    cmd3, maf = vcfanno(filtervcf, tumorid, normalid,
                   config.get("software", "VCF2MAF"),
                   config.get("Vcf2maf", "VEPPATH"),
                   config.get("Vcf2maf", "VEPDATA"),
                   config.get("ReferenceFasta", "REFFASTA"),
                   config.get("AnnoVCF", "ExAC"),
                   config.get("Vcf2maf", "ENST"),
                   maf)

    print >> script,"\n".join([cmd1, cmd2, cmd3])
    script.close()


def scalpelindel(tumorbam, normalbam, inifile, scriptout, outdir="Results/scalpel"):
    config = readconfig(inifile)
    script = open(scriptout, "w")
    tumorid, normalid = map(bam2id, [tumorbam, normalbam])
    tmpdir = "TMP/%s/scalpel" % (tumorid)
    tumorid, normalid = map(bam2id, [tumorbam, normalbam])
    cmd1, rawvcf = scalpel(tumorbam, normalbam,
                           config.get("software", "SCALPEL"),
                           config.get("Target", "EXT10"),
                           config.get("ReferenceFasta", "REFFASTA"),
                           config.get("Parameters2Scalpel", "numprocs"),
                           config.get("Parameters2Scalpel", "max-ins-size"),
                           config.get("Parameters2Scalpel", "max-del-size"),
                           config.get("Parameters2Scalpel", "min-alt-count-tumor"),
                           config.get("Parameters2Scalpel", "min-coverage-tumor"),
                           config.get("Parameters2Scalpel", "max-coverage-tumor"),
                           config.get("Parameters2Scalpel", "min-coverage-normal"),
                           config.get("Parameters2Scalpel", "max-coverage-normal"),
                           config.get("Parameters2Scalpel", "min-phred-fisher"),
                           config.get("Parameters2Scalpel", "min-vaf-tumor"),
                           tmpdir)


    #filter
    filtervcf = "%s/%s.vcf"%(outdir, tumorid)
    cmd2 = "%s/vcfilter.py scalpel %s %s %s %s "%(config.get("ScriptPath", "SOMATIC"),
                                            rawvcf, filtervcf, tumorid, normalid )

    #anno
    maf = "%s/%s.maf"%(outdir, tumorid)
    cmd3, maf = vcfanno(filtervcf, tumorid, normalid,
                   config.get("software", "VCF2MAF"),
                   config.get("Vcf2maf", "VEPPATH"),
                   config.get("Vcf2maf", "VEPDATA"),
                   config.get("ReferenceFasta", "REFFASTA"),
                   config.get("AnnoVCF", "ExAC"),
                   config.get("Vcf2maf", "ENST"),
                   maf)

    print >> script, "\n".join([cmd1, cmd2, cmd3])
    script.close()


def facetscnv(tumorbam, normalbam, inifile, scriptout, outdir="Results/facets"):
    config = readconfig(inifile)
    tumorid = tumorbam.split("/")[-1].split(".")[0]
    normalid = normalbam.split("/")[-1].split(".")[0]
    tmpdir = "TMP/%s/facets" % (tumorid)
    ugvcf = "%s/%s.ug.vcf" % (tmpdir, normalid)
    outtmp = "%s/raw" % outdir
    makedir(tmpdir, outdir, outtmp)

    script = open(scriptout, "w")

    cmd1, normalvcf = ug(normalbam,
                         config.get("software", "JAVA"),
                         config.get("JavaTmp", "JAVATMP"),
                         config.get("software", "GATK"),
                         config.get("ReferenceFasta", "REFFASTA"),
                         config.get("Target", "EXT10"),
                         config.get("Parameters2ug", "nct"),
                         ugvcf)
    snp = "%s/%s.snp" % (tmpdir, tumorid)
    cmd2, snpfile = snp_pileup(config.get("software", "SNP_PILEUP"),
                               normalvcf, normalbam, tumorbam,
                               config.get("Parameters2Facets", "min-map-quality"),
                               config.get("Parameters2Facets", "min-base-quality"),
                               config.get("Parameters2Facets", "pseudo-snps"),
                               config.get("Parameters2Facets", "min-read-counts"),
                               snp)
    cmd3, cnvseg = facets(config.get("ScriptPath", "somatic") + "/facets_purity.R",
                          tmpdir + "/", tumorid, tmpdir)

    seg = "%s/%s.seg" % (outdir, tumorid)
    avinput = "%s/%s.avinput" % (outtmp, tumorid)
    cmd4 = "%s/facets_filter.py %s %s.fai %s %s %s %s " % (
        config.get("ScriptPath", "SOMATIC"),
        cnvseg,
        config.get("ReferenceFasta", "REFFASTA"),
        config.get("Parameters2FilterFacets", "gainlog2"),
        config.get("Parameters2FilterFacets", "losslog2"),
        seg, avinput
    )

    avoutput = "%s/%s" %(outtmp, tumorid)
    cmd5 = cnvannovar(avinput, avoutput,
                      config.get("Annovar","TABLE_ANNOVAR"),
                      config.get("Annovar", "AnnovarDB"))

    cnvanno = "%s/%s.tsv"%(outdir, tumorid)
    cmd6 = "%s/combineCNV.py " \
            "%s %s.hg19_multianno.txt " \
            "> %s "%(config.get("ScriptPath", "SOMATIC"),
                    avinput, avoutput, cnvanno)

    cmd7 = "perl %s/circos.pl " \
           "-seg %s " \
           "-tumor %s " \
           "-outdir %s " \
           "-max %s " \
           "-min %s " \
           "-circos %s"%(
        config.get("ScriptPath", "SOMATIC"),
        seg, tumorid, outtmp,
        config.get("Parameters2FilterFacets", "gainlog2"),
        config.get("Parameters2FilterFacets", "losslog2"),
        config.get("software", "CIRCOS")
    )

    cmds = "\n".join([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7])
    print >> script, cmds
    script.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print "py tumorbam normalbam inifile output"
        exit()
    t = sys.argv
    snv = t[4] + ".snv.sh"
    indel = t[4] + ".indel.sh"
    cnv = t[4] + ".cnv.sh"
    mutectv1snv(t[1], t[2], t[3], snv)
    scalpelindel(t[1],t[2],t[3], indel)
    facetscnv(t[1],t[2],t[3], cnv)





