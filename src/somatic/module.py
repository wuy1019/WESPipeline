#!/usr/bin/env python
# File Name: src/somatic/module.py
# Author: wuy
# Created Time: Thu 16 Aug 2018 10:45:34 AM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import os

import numpy as np
import ConfigParser as configparser
from collections import defaultdict


def Ddict():
    return defaultdict(dict)


def makedir(*dirs):
    for path in dirs:
        if not os.path.exists(path):
            os.makedirs(path)


def checkFile(fileToCheck):
    if os.path.isfile(fileToCheck):
        return True
    else:
        return False


def readconfig(inifile):
    config = configparser.ConfigParser()
    config.read(inifile)
    config.read(config.get("common", "config"))
    return config


def bam2id(bam):
    return bam.split("/")[-1].split(".")[0]


def mutectv1(tumorbam,
             normalbam,
             java,
             javatmp,
             mutect,
             ref,
             fraction_contamination,
             intervals,
             min_qscore,
             cosmic,
             dbsnp,
             tmpdir,
             outdir="Results/mutectv1"):
    tumorid = bam2id(tumorbam)

    ouDir = "%s/raw" % outdir
    result = "%s/%s.mutectv1.vcf" % (ouDir, tumorid)
    makedir(tmpdir, ouDir)
    cmd = "%s -Djava.io.tmpdir=%s " \
          "-jar %s --analysis_type MuTect " \
          "--reference_sequence %s " \
          "--input_file:tumor %s " \
          "--input_file:normal %s " \
          "--out %s/%s.muetctv1.txt " \
          "--vcf %s/%s.mutectv1.vcf " \
          "--coverage_file %s/%s.muetctv1.wig.txt " \
          "--fraction_contamination %s " \
          "--intervals %s " \
          "--min_qscore %s " \
          "--cosmic %s " \
          "--dbsnp %s " \
          "--disable_auto_index_creation_and_locking_when_reading_rods " \
          "--tumor_sample_name TUMOR --normal_sample_name CONTROL " % (
              java, javatmp, mutect, ref, tumorbam, normalbam, tmpdir, tumorid,
              ouDir, tumorid, tmpdir, tumorid,  fraction_contamination,
              intervals, min_qscore, cosmic, dbsnp
          )

    return cmd, result


def scalpel(tumorbam,
            normalbam,
            scalpeldir,
            intervals,
            ref,
            numprocs,
            max_ins_size,
            max_del_size,
            min_alt_count_tumor,
            min_coverage_tumor,
            max_coverage_tumor,
            min_coverage_normal,
            max_coverage_normal,
            min_phred_fisher,
            min_vaf_tumor,
            tmpdir,
            outdir="Results"):
    tumorid = tumorbam.split("/")[-1].split(".")[0]

    ouDir = "%s/scalpel/raw" % outdir
    makedir(tmpdir, ouDir)

    result = "%s/%s.scalpel.vcf" % (ouDir, tumorid)
    cmd1 =  "%s/scalpel-discovery --somatic " \
            "--tumor %s " \
            "--normal %s " \
            "--bed %s " \
            "--ref %s " \
            "--numprocs %s " \
            "--two-pass " \
            "--dir %s "%(
        scalpeldir, tumorbam, normalbam,intervals, ref, numprocs, tmpdir
    )
    cmd2 =  "%s/scalpel-export --somatic "\
            "--db %s/twopass/somatic.db "\
            "--bed  %s "\
            "--ref %s "\
            "--output-format vcf "\
            "--max-ins-size %s "\
            "--max-del-size %s "\
            "--min-alt-count-tumor %s "\
            "--min-coverage-tumor %s "\
            "--max-coverage-tumor %s "\
            "--min-coverage-normal %s "\
            "--max-coverage-normal %s "\
            "--min-phred-fisher %s "\
            "--min-vaf-tumor %s " \
            "> %s "%(
        scalpeldir, tmpdir, intervals, ref, max_ins_size, max_del_size,
        min_alt_count_tumor, min_coverage_tumor, max_coverage_tumor,
        min_coverage_normal, max_coverage_normal, min_phred_fisher,
        min_vaf_tumor, result
    )

    cmds = "&&".join([cmd1, cmd2])
    return cmds, result


def vardict(tumorbam,
            normalbam,
            vardictjava,
            vardictperlpath,
            ref,
            af,
            intervals,
            th,
            outdir="Results/mutectv1"):

    tumorid, normalid = map(bam2id, [tumorbam, normalbam])
    outraw = "%s/raw" % outdir
    makedir(outraw)
    result = "{dir}/{tumorid}.vardict.vcf".format(tumorid=tumorid, dir=outraw)
    cmd = "{vardictjava} -G {ref} -f {af} -N {tumorid} " \
           "-b \"{tumorbam}|{normalbam}\" -z -c 1 -S 2 -E 3 " \
           "-th {th} {intervals} | {vardictperlpath}/testsomatic.R " \
           "| {vardictperlpath}/var2vcf_paired.pl -N \"{tumorid}|{normalid}\" " \
           "-f {af} > {dir}/{tumorid}.vardict.vcf" \
           "".format(vardictjava=vardictjava, ref=ref, af=af, tumorid=tumorid,
                     normalid=normalid, tumorbam=tumorbam, normalbam=normalbam,
                     th=th, intervals=intervals, vardictperlpath=vardictperlpath,
                     dir = outraw)
    return cmd, result


def ug(inbam, java, javatmp, gatk, ref, intervals, nct, result):
    cmd = "%s -Djava.io.tmpdir=%s -jar %s -T UnifiedGenotyper " \
          "-R %s -I %s -glm BOTH --min_indel_count_for_genotyping 3 " \
          "--min_indel_fraction_per_sample 0.05 -stand_emit_conf 1 " \
          "-stand_call_conf 1 -L %s -o %s -ploidy 2 -nt %s " \
          "--disable_auto_index_creation_and_locking_when_reading_rods "%(
        java, javatmp, gatk, ref, inbam, intervals, result, nct
    )
    return cmd, result


def snp_pileup(snp_pileup, normalvcf, normalbam, tumorbam, min_map_quality,
               min_base_qualit, pseudo_snps, min_read_counts, output):

    cmd = "%s -q %s -Q %s -P %s -r %s " \
           "%s %s %s %s "%(
        snp_pileup, min_map_quality, min_base_qualit, pseudo_snps,
        min_read_counts, normalvcf, output, normalbam, tumorbam
    )
    return cmd, output


def facets(facetsR, snpdir, tumorid, outdir):
    makedir(outdir)
    result = "%s/%s.seg" % (outdir, tumorid)
    cmd = "Rscript %s %s %s %s" % (facetsR, snpdir, tumorid, outdir + "/")
    return cmd, result


def vcfanno(invcf, tumorid, normalid, vcf2maf, veppath, vepdata, ref, exacvcf,
            custom_enst, outmaf):

    cmd = "perl %s --input-vcf %s --output-maf %s " \
          "--tumor-id %s --normal-id %s " \
          "--vcf-tumor-id TUMOR --vcf-normal-id CONTROL " \
          "--retain-info --any-allele " \
          "--vep-path %s --vep-data %s " \
          "--ref-fasta %s " \
          "--filter-vcf %s --custom-enst %s "%(
        vcf2maf, invcf, outmaf, tumorid, normalid, veppath, vepdata, ref, exacvcf, custom_enst
    )

    return cmd, outmaf


def facetsfilter(inseg, reffai, gainlog2, losslog2, outseg):
    faidic = {}
    with open(reffai) as fai:
        for i in fai:
            t = i.rstrip().split()
            faidic[t[0]] = t[1]
    faidic["23"] = faidic["X"]
    faidic["24"] = faidic["Y"]

    segfile = open(inseg).readlines()
    out = open(outseg, "w")
    Median = np.median(
        map(float, list(np.loadtxt(inseg, str, delimiter=',')[:, 4][1:])))
    out.write(segfile[0] + "\n")
    for i in segfile[1:]:
        t = i.split(",")
        log2 = float(t[4]) - Median
        start, end = map(int, map(float, [t[9], t[10]]))
        if end > faidic[t[0]]:
            end = faidic[t[0]]
        if log2 <= losslog2 or log2 >= gainlog2:
            out.write(",".join(
                map(str, t[:4] + [log2] + t[5:9] + [start, end] + t[11:])) +
                      "\n")
    out.close()
    return outseg


def cnvannovar(avinput, avoutput, table_annovar, annovardb):

    cmd = "perl %s " \
          "%s %s -buildver hg19 " \
          "-out %s " \
          "-remove -protocol " \
          "refGene,cytoBand " \
          "-operation " \
          "g,r " \
          "-nastring . "%(
        table_annovar, avinput, annovardb, avoutput
    )
    return cmd


###### maf combine and filter module ######


def txt2dic(txt, headidx=0, keyidx=[1, 2, 3, 4]):
    dic = Ddict()
    header = open(txt).readlines()[headidx].rstrip().split("\t")
    for i in open(txt).readlines()[headidx + 1:]:
        t = i.rstrip().split("\t")
        key = tuple([t[x] for x in keyidx])
        for s in range(len(header)):
            dic[key][header[s]] = t[s]
    return dic, header


def str2float(s):
    if not s:
        return 0
    else:
        return float(s)


def popfiltered(frequency, poplis):
    if len(filter(lambda x: x >= frequency, map(float, poplis))):
        return False
    else:
        return True


def maffilter(mafdic, type_header, type_list, freq_header, freq):
    """
    :param maf: file
    :param type_header: str
    :param type_list: list
    :param freq_header: list
    :param freq: float
    :return: defaultdict, filter maf
    """
    dic = mafdic
    end = Ddict()
    for key, val in dic.items():
        #filter gene symbol source, only reserve HGNC
        if val["SYMBOL_SOURCE"] != "HGNC":
            continue
        if val[type_header] not in type_list:
            continue
        popfreq = [str2float(val[x] for x in freq_header)]
        if not popfiltered(freq, popfreq):
            continue
        end[key] = dic[key]
    return end


def combineMAF(maf1, maf2, output):
    version = open(maf1).readlines()[0].rstrip()
    header = open(maf1).readlines()[1].rstrip().split("\t")
    maf1raw, header1 = txt2dic(maf1, headidx=1, keyidx=[4, 5, 10, 12])
    maf2raw, header2 = txt2dic(maf2, headidx=1, keyidx=[4, 5, 10, 12])
    mafraw = maf1raw.update(maf2raw)
    out = open(output, "w")
    out.write(version + "\n")
    out.write("\t".join(header) + "\n")
    for key, val in mafraw.items():
        lis = []
        for h in header:
            lis.append(val[h])
        out.write("\t".join(lis) + "\n")
    out.close()
    return output


def combineFilterMAF(maf, type_header, type_list, freq_header, freq, output):
    version = open(maf).readlines()[0].rstrip()
    header = open(maf).readlines()[1].rstrip().split("\t")
    mafraw, header = txt2dic(maf, headidx=1, keyidx=[4, 5, 10, 12])
    maf_filter = maffilter(mafraw, type_header, type_list, freq_header, freq)
    out = open(output, "w")
    out.write(version + "\n")
    out.write("\t".join(header) + "\n")
    for key, val in maf_filter.items():
        lis = []
        for h in header:
            lis.append(val[h])
        out.write("\t".join(lis) + "\n")
    out.close()
    return output


## covert format
def covertMAF(inmaf, outmaf):
    header = [
        "Hugo_Symbol", "Gene.ID", "HGVS", "Chr.start", "Chr.end",
        "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
        "Variant_Classification", "dbSNP_RS", "t-AF", "t_DP(ref:alt)", "n-AF",
        "n_DP(ref:alt)", "Entrez_Gene_Id", "NCBI_Build", "Strand",
        "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
        "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
        "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
        "Verification_Status", "Validation_Status", "Mutation_Status",
        "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score",
        "BAM_File", "Sequencer", "Tumor_Sample_UUID",
        "Matched_Norm_Sample_UUID", "HGVSp", "Transcript_ID", "all_effects",
        "Allele", "Gene", "Feature", "Feature_type", "Consequence",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND_VEP",
        "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "CANONICAL", "CCDS",
        "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "SIFT", "PolyPhen", "INTRON",
        "DOMAINS", "AF", "AFR_AF", "AMR_AF", "ASN_AF", "EAS_AF", "EUR_AF",
        "SAS_AF", "AA_AF", "EA_AF", "CLIN_SIG", "SOMATIC", "PUBMED",
        "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
        "IMPACT", "PICK", "VARIANT_CLASS", "TSL", "HGVS_OFFSET", "PHENO",
        "MINIMISED", "ExAC_AF", "ExAC_AF_AFR", "ExAC_AF_AMR", "ExAC_AF_EAS",
        "ExAC_AF_FIN", "ExAC_AF_NFE", "ExAC_AF_OTH", "ExAC_AF_SAS",
        "GENE_PHENO", "FILTER", "flanking_bps", "variant_id", "variant_qual",
        "ExAC_AF_Adj", "ExAC_AC_AN_Adj", "ExAC_AC_AN", "ExAC_AC_AN_AFR",
        "ExAC_AC_AN_AMR", "ExAC_AC_AN_EAS", "ExAC_AC_AN_FIN", "ExAC_AC_AN_NFE",
        "ExAC_AC_AN_OTH", "ExAC_AC_AN_SAS", "ExAC_FILTER", "gnomAD_AF",
        "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF",
        "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"
    ]
    mafdic, header = txt2dic(inmaf, headidx=1, keyidx=[4, 5, 10, 12])
    out = open(outmaf, "w")
    out.write("#version 2.4\n")
    out.write("\t".join(header) + "\n")

    for key, val in mafdic.items():
        val["Gene.ID"] = ", ".join([val["RefSeq"], val["EXON"]])
        if val["HGVSp_Short"]:
            val["HGVS"] = "{}({})".format(val["HGVSc"], val["HGVSp_Short"])
        else:
            val["HGVS"] = val["HGVSc"]
        val["Chr.start"] = "{}:{}".format(val["Chromosome"],
                                          val["Start_Position"])
        val["Chr.end"] = "{}:{}".format(val["Chromosome"], val["End_Position"])
        val["t-AF"] = round(
            float(val["t_alt_count"]) / float(val["t_depth"]), 2)
        val["n-AF"] = round(
            float(val["n_alt_count"]) / float(val["n_depth"]), 2)
        val["t_DP(ref:alt)"] = "{}({}:{})".format(
            val["t_depth"], val["t_ref_count"], val["t_alt_count"])
        val["n_DP(ref:alt)"] = "{}({}:{})".format(
            val["n_depth"], val["n_ref_count"], val["n_alt_count"])

        lis = []
        for h in header:
            lis.append(val[h])
        out.write("\t".join(lis) + "\n")
    out.close()
    return outmaf
