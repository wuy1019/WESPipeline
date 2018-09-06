#!/usr/bin/env python
# File Name: GSGP_Script_301_Anno.py
# Author: Wuy
# Created Time: Sat 01 Apr 2017 11:21:14 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################

import module
from module import cut


def ref_alt(ref, alt):
    """
    """
    Ref = cut(ref, alt)[0]
    Alt = cut(ref, alt)[1]
    if Ref == "":
        Ref = "-"
    if Alt == "":
        Alt = "-"
    return Ref, Alt


def genetic_anno_combine(vcf_file, annovar, vep, intervar, omim, hgmd):
    #def genetic_anno_combine(vcf_file, annovar, vep, intervar, omim):
    dic = module.Ddict()
    header = module.header()
    annovardic = module.tsvparse2(
        annovar,
        keyidx=[0, 1, 3,
                4])  #chrom with chr string; key = "chrom,start,end,ref,alt"
    rawvep = module.vep(
        vep)  #chrom with chr string; key = "chrom,start,ref,alt,hgvs"
    rawinter = module.parseintervar(
        intervar)  #chrom without chr string; key = "chrom,start,ref,alt"
    hgmd_snp, hgmd_indel = module.parsehgmd(hgmd)  #
    hgmd_gene2tran = module.hgmdgene2trans(hgmd)  # key: gene, val: trans
    rawomim = module.parsegenemap2ClnSynopsis(
        omim, key="Symbol")  # key = gene,
    #print rawomim
    vcfdic_raw = module.hcvcf2dic(
        vcf_file
    )  #chrom with chr string; key = (chrom, pos, ref, alt), vcf file raw value

    #parse vcf file
    vcfdic = module.Ddict()
    for key in vcfdic_raw:
        chrom, pos, ref, alt = key
        Ref, Alt = ref_alt(ref, alt)
        #snp
        if len(Ref) == len(Alt) and Ref != "-" and Alt != "-":
            Pos = pos
        #ins
        elif len(Ref) < len(Alt) or Ref == "-":
            Pos = pos
        #del
        elif len(Ref) > len(Alt) or Alt == "-":
            Pos = int(pos) + (len(ref) - len(Ref))
        #del
    # elif len(Ref) == 1 and Alt == "-":

    #    Pos = int(pos) + (len(ref)-len(Ref))
    #other
        else:
            Pos = pos
        vcfkey = ",".join(map(str, [chrom, Pos, Ref, Alt]))
        vcfdic[vcfkey] = vcfdic_raw[key]
        #if int(Pos) == 16327915:
        #    print vcfkey, vcfdic[vcfkey]
    #parse vep file
    vepdic = module.Ddict()
    vepkeydic = {}
    for key in rawvep:
        tmp = ",".join(key.split(",")[:4])
        if tmp not in vepkeydic:

            vepkeydic[tmp] = [key]
        else:
            vepkeydic[tmp] += [key]

    for key, val in vepkeydic.items():
        genes = []
        hgvss = []
        types = []
        exons = []
        for v in val:
            genes.append(rawvep[v][0])
            hgvss.append(rawvep[v][2])
            types.append(rawvep[v][1])
            exons.append(rawvep[v][3])
        vepdic[key]["Gene"] = ";".join(module.uniqlist(genes))
        vepdic[key]["HGVS"] = ";".join(hgvss)
        vepdic[key]["Mutation_type"] = ";".join(types)
        vepdic[key]["ExonNum"] = ";".join(exons)

    # init
    for key in annovardic.keys():
        for item in header:
            dic[key][item] = "."

    #parse annovar file
    for key, val in annovardic.items():

        for i in val:
            dic[key][i] = val[i]

        dic[key]["Chr.Start"] = ":".join([val["Chr"], val["Start"]])
        dic[key]["Chr.End"] = ":".join([val["Chr"], val["End"]])
        dic[key]["Clinvar"] = val["CLINSIG"]

        # add vcf info
        if key in vcfdic:
            dic[key]["VAF"] = vcfdic[key]["AF"]
            dic[key]["Format"] = ":".join(
                map(str, [
                    vcfdic[key]["DP"], vcfdic[key]["AD"],
                    vcfdic[key]["DP"] - vcfdic[key]["AD"]
                ]))
            if vcfdic[key]["DP"] - vcfdic[key]["AD"] > 1:
                dic[key]["Zygosity"] = "Het"
            else:
                dic[key]["Zygosity"] = "Hom"

        #add vep info
        if key in vepdic:
            dic[key]["Gene"] = vepdic[key]["Gene"]
            dic[key]["HGVS"] = vepdic[key]["HGVS"]
            dic[key]["Mutation_type"] = vepdic[key]["Mutation_type"]
            dic[key]["ExonNum"] = vepdic[key]["ExonNum"]
            #ad omim
            genelis = vepdic[key]["Gene"].split(";")
            diseases = []
            PhenotypeMIMs = []
            geneMIMs = []
            clinicalSynopsiss = []
            Inheritances = []
            hgmd_trans = []
            for gene in genelis:
                if gene in rawomim:
                    diseases += rawomim[gene]["Phenotype"]
                    PhenotypeMIMs += rawomim[gene]["PhenotypeMIM"]
                    geneMIMs += rawomim[gene]["geneMIM"]
                    clinicalSynopsiss += rawomim[gene]["clinicalSynopsis"]
                    Inheritances += rawomim[gene]["Inheritance"]
                if gene in hgmd_gene2tran:
                    hgmd_trans += [hgmd_gene2tran[gene]]
            if diseases:
                dic[key]["Disease"] = ";".join(diseases)
            if PhenotypeMIMs:
                dic[key]["PhenotypeMIM"] = ";".join(PhenotypeMIMs)
            if geneMIMs:
                dic[key]["geneMIM"] = ";".join(geneMIMs)
            if Inheritances:
                dic[key]["Inheritance"] = ";".join(
                    module.uniqlist(map(module.capital, Inheritances)))
            if clinicalSynopsiss:
                tmp = ""
                for i in range(1, len(clinicalSynopsiss) + 1):
                    tmp += "$$$%s.%s" % (i, clinicalSynopsiss[i - 1])
                dic[key]["clinicalSynopsis"] = tmp
            if hgmd_trans:
                dic[key]["HGMD_trans"] = ";".join(hgmd_trans)

        #add intervar
        rmchrkey = key.lstrip("chr")
        if rmchrkey in rawinter:
            dic[key]["intervar"] = rawinter[rmchrkey]["intervar"]
            dic[key]["ACMG"] = rawinter[rmchrkey]["ACMG"]
            dic[key]["Clinvar_intervar"] = rawinter[rmchrkey]["clinvar"]
            dic[key]["OrphaNumber"] = rawinter[rmchrkey]["OrphaNumber"]
            dic[key]["Orpha"] = rawinter[rmchrkey]["Orpha"]

        #add hgmd
        hgmd_chr = key.split(",")[0].lstrip("chr")
        hgmd_pos = key.split(",")[1]
        dic[key]["HGMD_mutation_type"], dic[key]["HGMD_variant_class"], dic[
            key]["HGMD_disease"], dic[key]["HGMD_pubmed_ID"] = module.hgmdpos(
                hgmd_chr, hgmd_pos, "0", hgmd_snp)
        dic[key]["HGMD_indel_type"], dic[key]["HGMD_Indel_class"], dic[key][
            "HGMD_indel_disease"], dic[
                key]["HGMD_indel_PMID"] = module.hgmdpos(
                    hgmd_chr, hgmd_pos, "0", hgmd_indel)

        #distance 1
        snp_dis1 = ",".join(module.hgmdpos(chrom, pos, "1", hgmd_snp))
        snp_dis_1 = ",".join(module.hgmdpos(chrom, pos, "-1", hgmd_snp))
        indel_dis1 = ",".join(module.hgmdpos(chrom, pos, "1", hgmd_indel))
        indel_dis_1 = ",".join(module.hgmdpos(chrom, pos, "-1", hgmd_indel))

        #distance 2
        snp_dis2 = ",".join(module.hgmdpos(chrom, pos, "2", hgmd_snp))
        snp_dis_2 = ",".join(module.hgmdpos(chrom, pos, "-2", hgmd_snp))
        indel_dis2 = ",".join(module.hgmdpos(chrom, pos, "2", hgmd_indel))
        indel_dis_2 = ",".join(module.hgmdpos(chrom, pos, "-2", hgmd_indel))

        #distance 3
        snp_dis3 = ",".join(module.hgmdpos(chrom, pos, "3", hgmd_snp))
        snp_dis_3 = ",".join(module.hgmdpos(chrom, pos, "-3", hgmd_snp))
        indel_dis3 = ",".join(module.hgmdpos(chrom, pos, "3", hgmd_indel))
        indel_dis_3 = ",".join(module.hgmdpos(chrom, pos, "-3", hgmd_indel))

        dis1 = module.uniqlist([snp_dis1, snp_dis_1, indel_dis1, indel_dis_1])
        dis2 = module.uniqlist([snp_dis2, snp_dis_2, indel_dis2, indel_dis_2])
        dis3 = module.uniqlist([snp_dis3, snp_dis_3, indel_dis3, indel_dis_3])
        dis1.remove(".,.,.,.")
        dis2.remove(".,.,.,.")
        dis3.remove(".,.,.,.")

        if dis1:
            dic[key]["HGMD_nearby_alleles"] = ";".join(dis1)
        elif dis2:
            dic[key]["HGMD_nearby_alleles"] = ";".join(dis2)
        elif dis3:
            dic[key]["HGMD_nearby_alleles"] = ";".join(dis3)
        else:
            dic[key]["HGMD_nearby_alleles"] = "."

    return dic


def writeanno(vcf_file, annovar, vep, intervar, omim, hgmd, output):
    w = open(output, "w")
    dic = genetic_anno_combine(vcf_file, annovar, vep, intervar, omim, hgmd)
    header = module.header()
    w.write("\t".join(header) + "\n")
    for key, val in dic.items():
        lis = []
        for i in header:
            lis.append(val[i])
        w.write("\t".join(map(str, lis) + "\n"))
    w.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print "py vcf annovar vep intervar omim hgmd"
        exit(1)
    t = sys.argv
    dic = genetic_anno_combine(t[1], t[2], t[3], t[4], t[5], t[6])

    header = module.header()
    print "\t".join(header)
    for key, val in dic.items():
        lis = []
        for i in header:
            lis.append(val[i])
        print "\t".join(map(str, lis))
