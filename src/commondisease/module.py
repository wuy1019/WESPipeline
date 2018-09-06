#!/usr/bin/env python
# File Name: module.py
# Author: wuy
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import csv
import re
from collections import Counter, defaultdict

import vcf

#import os

def Ddict():
    return defaultdict(dict)

def capital(s):
    return s.capitalize()

def conf2dic(infile):
        conf_dic = {}
        for conf in open(infile):
                tmp = conf.split()
                key = tmp[1]
                val = tmp[2]
                conf_dic[key] = val
        return conf_dic


def uniqlist(rawlis):
    func = lambda x,y:x if y in x else x + [y]
    return reduce(func, [[], ] + rawlis)


def header():
    lis = ['intervar', 'Clinvar', 'Clinvar_intervar', 'HGMD_variant_class', 'HGMD_Indel_class', 'Chr.Start', 'Chr.End',
            'Ref', 'Alt', 'Zygosity', 'VAF', 'Format', 'Gene', 'geneMIM', 'Mutation_type', 'ExonNum', 'HGVS',
            'HGMD_trans', 'PhenotypeMIM', 'Disease', 'Inheritance', 'HGMD_mutation_type', 'HGMD_disease',
            'HGMD_pubmed_ID', 'HGMD_indel_type', 'HGMD_indel_disease', 'HGMD_indel_PMID',
            'HGMD_nearby_alleles', 'snp138', '1000G_ALL', '1000G_EAS', 'ExAC_ALL', 'ExAC_EAS',
            'gnomAD_exome_ALL', 'gnomAD_exome_EAS', 'gnomAD_genome_ALL', 'gnomAD_genome_EAS',
            'ESP6500siv2_ALL', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred',
            'CADD_raw', 'CADD_phred', 'GERP++_RS', 'ACMG', 'Orpha', 'OrphaNumber', 'clinicalSynopsis',
            'Chinese_Name', 'English_Name', 'SC_GRS_HGMD', 'SC_GRS_GWAS', 'OR_GRS_GWAS', 'In_HGMD', 'In_GWAS',
            'total_num_pos_in_hgmd', 'total_num_pos_in_gwas', 'risk_num_pos_in_hgmd', 'risk_num_pos_in_gwas', 'final_state']

    return lis

def tsvparse(tsv):
    dic = Ddict()

    head = open(tsv).readline().rstrip().split("\t")
    chrix = head.index("CHR")
    posix = head.index("POS")
    refix = head.index("REF")
    altix = head.index("ALT")

    for i in open(tsv).readlines()[1:]:
        i = i.rstrip()
        t = i.split("\t")
        chrom = t[chrix]
        pos   = t[posix]
        ref   = t[refix]
        alt   = t[altix]
        hgvs  = t[head.index("HGVS nomenclature")]
        key = ",".join([chrom, pos, ref, alt, hgvs])

        for s in range(len(head)):
            dic[key][head[s]] = t[s]
    return dic


def tsvparse2(tsv, keyidx=[0,1,2,3]):
    """
    key: join in index of headline
    example: keyidx = [0,1,2], key = ",".join([headline[0], headline[1], headline[2]])
    """
    dic = Ddict()
    head = open(tsv).readline().rstrip().split("\t")
    for i in open(tsv).readlines()[1:]:
        i = i.rstrip()
        t = i.split("\t")
        key = ",".join([t[x] for x in keyidx])
        for s in range(len(head)):
            dic[key][head[s]] = t[s]
    return dic

def polarchange(hgvs):
    dic = {
            "Ala": "NP",
            "Val": "NP",
            "Leu": "NP",
            "Ile": "NP",
            "Pro": "NP",
            "Phe": "NP",
            "Trp": "NP",
            "Met": "NP",
            "Gly": "P0",
            "Ser": "P0",
            "Thr": "P0",
            "Cys": "P0",
            "Tyr": "P0",
            "Asn": "P0",
            "Gln": "P0",
            "Asp": "P-",
            "Glu": "P-",
            "Lys": "P+",
            "Arg": "P+",
            "His": "P+",
            "Ter": "."
            }
    polar = "."
    p = re.findall("\(p\.(\w{3})\d+(\w{3})\)", hgvs)
    if p:
        if p[0][0] in dic and p[0][1] in dic:
            p1 = dic[p[0][0]]
            p2 = dic[p[0][1]]
            polar = "=>".join([p1, p2])
    return polar

# one gene can output multui transcripts
def vep(vepcsv):
    """key:"chrom,pos,ref,alt,hgvs"
       val:gene, mut_type, hgvs, exonnumber, geneflag"""
    dic_vep = Ddict()
    for i in open(vepcsv):
        if i.startswith("#"):
            continue
        i = i.rstrip()
        t = i.split(",")
        chrom = t[4].split(":")[0]
        pos = t[4].split(":")[-1]
        ref = t[6]
        alt = t[7]
        gene = t[1]
        mut_type = t[8].split(";")[0]
        hgvs = "."
        if "dup" not in t[-1].split("(")[0]:
            hgvs = t[-1]

#        trans = "-"
#        if "dup" not in hgvs.split("(")[0] and hgvs != ".":
#            trans = hgvs.split("(")[0]
#        else:
#            trans = "-"
        exonnumber = "."
        geneflag   = "NO"
        if hgvs != ".":
            if "exon" in t[2] or "intron" in t[2]:
                exonnumber = t[2].split(":")[-2]+":"+t[2].split(":")[-1]
                geneflag = "YES"
            #else:
            #    geneflag = "NO"
            #    exonnumber = "."

        key1 = ",".join([chrom, pos, ref, alt])
        dic_vep[key1][hgvs] = [gene, mut_type, hgvs, exonnumber, geneflag]

    dic = {}
    for key in dic_vep:
        NMlis  = []
        other  = []
        dotlis = []
        translis = dic_vep[key].keys()
        for i in translis:
            if i.startswith("NM"):
                NMlis.append(i)
            elif i == ".":
                dotlis.append(i)
            else:
                other.append(i)
        if NMlis:
            for n in NMlis:
                tmp = key+","+n
                dic[tmp] = dic_vep[key][n]
        elif other:
            tmp = key+","+ other[0]
            dic[tmp] = dic_vep[key][other[0]]
        else:
            tmp = key+","+"."
            dic[tmp] = dic_vep[key]["."]
    return dic



def printtsv(dic):
    head = header()
    print "\t".join(head)
    for key, val in dic.items():
        lis = []
        for i in head:
            lis.append(val[i])
        print "\t".join(lis)


def parsegenemap2ClnSynopsis(genemap2ClnSynopsis, key="EntrezID"):

    #key can set for EntrezID or Gene Symbol
    dic = Ddict()
    def rmlis(lis):
        return [x for x in lis if x]
    with open(genemap2ClnSynopsis) as fin:
        for i in fin:
            i = i.rstrip("\n")
            if i.startswith("#"):
                continue
            t = i.split("\t")

            symbol = t[8].strip()
            geneid = t[9].strip()
            if key == "EntrezID":
                idkey = geneid
            elif key == "Symbol":
                idkey = symbol
            else:
                raise ValueError("parsegenemap2ClnSynopsis: parameter of key must be EntrezID or Symbol")

            if geneid:
                if geneid not in dic:
                    dic[idkey]["geneMIM"] = rmlis([t[5].strip()])
                    dic[idkey]["Symbol"] = rmlis([t[8].strip()])
                    dic[idkey]["Phenotype"] = rmlis(t[14].strip().split(";"))
                    dic[idkey]["PhenotypeMIM"] = rmlis(t[15].strip().split(";"))
                    dic[idkey]["PhenotypeKey"] = rmlis(t[16].strip().split(";"))
                    dic[idkey]["Inheritance"] = rmlis(t[17].strip().split(";"))
                    dic[idkey]["clinicalSynopsis"] = rmlis(t[18].lstrip("$").strip().split("$$$") )
                   #dic[idkey]["PhenoLink"] = link

                else:
                    dic[idkey]["geneMIM"] += rmlis([t[5].strip()])
                    dic[idkey]["Symbol"] += rmlis([t[8].strip()])
                    dic[idkey]["Phenotype"] += rmlis(t[14].strip().split(";"))
                    dic[idkey]["PhenotypeMIM"] += rmlis(t[15].strip().split(";"))
                    dic[idkey]["PhenotypeKey"] += rmlis(t[16].strip().split(";"))
                    dic[idkey]["Inheritance"] += rmlis(t[17].strip().split(";"))
                    dic[idkey]["clinicalSynopsis"] += rmlis(t[18].lstrip("$").strip().split("$$$"))
                  # dic[idkey]["PhenoLink"] += link
    links = []
    for key, val in dic.items():
        links = []
        if val["PhenotypeMIM"]:
            for i in val["PhenotypeMIM"]:
                links.append("https://www.omim.org/entry/"+i)
        dic[key]["PhenoLink"] = links

    return dic
def parsegenemap2(genemap2):
    dic = {}
    #key:gene
    #val:pheno, type is list
    with open(genemap2) as fin:
        for i in fin:
            i = i.rstrip("\n")
            if i.startswith("#"):
                continue
            t = i.split("\t")
            geneid, pheno = t[9].strip(), t[12].strip()
            if geneid and pheno:
                if geneid not in dic:
                    dic[geneid] = [pheno]
                else:
                    dic[geneid] += [pheno]
    return dic

def hgmdgene2trans(hgmd19):
    """
    key: gene
    val: transcript
    """
    dic = {}
    with open(hgmd19, "rb") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene = row["gene"]
            trans = row["refseq"]
            dic[gene] = trans
    return dic

def parsehgmd(hgmd19):
    dic_snp = Ddict()
    dic_indel = Ddict()
    # deltion of hgmd, the postion must be minus 1 can match our tsv postion
    with open(hgmd19, "rb") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["hg19_chr"] and row["hg19_start"]:

#               key = (row["hg19_chr"], row["hg19_start"])
                hgvsctmp = row["DNA"]
                hgvsptmp = row["PROT"]
                gene = row["gene"]
                hgvsc = ""
                hgvsp = ""
                NM = ""
                ref = row["hg19_ref"]
                alt = row["hg19_alt"]
                chrom, start = row["hg19_chr"], row["hg19_start"]
                if len(ref) > len(alt):
                    pos = str(int(row["hg19_start"]) - 1)
                else:
                    pos = start
                key = (chrom, pos)
                mut = "%s:%s_%s_%s"%(chrom, pos, ref, alt)
                tag = row["tag"]
                if hgvsctmp:
                    NM, hgvsc = hgvsctmp.split(":")
                if hgvsptmp and 'p.' in hgvsptmp:
                    NP, hgvsp = hgvsptmp.split(":")
                if hgvsc and hgvsp:
                    hgvs = "%s(%s):%s(%s)"%(NM, gene, hgvsc, hgvsp)
                elif hgvsc:
                    hgvs = "%s(%s):%s"%(NM, gene, hgvsc)
                else:
                    hgvs = ""

                if len(row["hg19_ref"])  == len(row["hg19_alt"]):
                    if key not in dic_snp:
                        #dic_snp[key] = Ddict()
                        dic_snp[key]["gene"] = [row["gene"]]
                        dic_snp[key]["disease"] = [row["disease"]]
                        dic_snp[key]["pmid"] = [row["pmid"]]
                        dic_snp[key]["acc_num"] = [row["acc_num"]]
                        #dic_snp[key]["hgvsc"] = [row["DNA"]]
                        #dic_snp[key]["hgvsp"] = [row["PROT"]]
                        dic_snp[key]["hgvs"] = [hgvs]
                        dic_snp[key]["mut"] = [mut]
                        dic_snp[key]["tag"] = [tag]
                    else:
                        dic_snp[key]["gene"] += [row["gene"]]
                        dic_snp[key]["disease"] += [row["disease"]]
                        dic_snp[key]["pmid"] += [row["pmid"]]
                        dic_snp[key]["acc_num"] += [row["acc_num"]]
                        #dic_snp[key]["hgvsc"] = [row["DNA"]]
                        #dic_snp[key]["hgvsp"] = [row["PROT"]]
                        dic_snp[key]["hgvs"] += [hgvs]
                        dic_snp[key]["mut"] += [mut]
                        dic_snp[key]["tag"] += [tag]
                else:
                    if key not in dic_indel:
                        #dic_snp[key] = Ddict()
                        dic_indel[key]["gene"] = [row["gene"]]
                        dic_indel[key]["disease"] = [row["disease"]]
                        dic_indel[key]["pmid"] = [row["pmid"]]
                        dic_indel[key]["acc_num"] = [row["acc_num"]]
                        #dic_indel[key]["hgvsc"] = [row["DNA"]]
                        #dic_indel[key]["hgvsp"] = [row["PROT"]]
                        dic_indel[key]["hgvs"] = [hgvs]
                        dic_indel[key]["mut"] = [mut]
                        dic_indel[key]["tag"] = [tag]
                    else:
                        dic_indel[key]["gene"] += [row["gene"]]
                        dic_indel[key]["disease"] += [row["disease"]]
                        dic_indel[key]["pmid"] += [row["pmid"]]
                        dic_indel[key]["acc_num"] += [row["acc_num"]]
                        #dic_snp[key]["hgvsc"] = [row["DNA"]]
                        #dic_snp[key]["hgvsp"] = [row["PROT"]]
                        dic_indel[key]["hgvs"] += [hgvs]
                        dic_indel[key]["mut"] += [mut]
                        dic_indel[key]["tag"] += [tag]
            else:
                pass
    return dic_snp, dic_indel

def hgmdpos(chrom, pos, distance, dic):
    """
    dic from parsehgmd def
    """
    mut_type, mut_class, mut_disease, mut_pubmed = ".", ".", ".", "."
    key = (chrom, str(int(pos)+int(distance)))
    if key in dic:
        mut_type, mut_class, mut_disease, mut_pubmed = \
                uniqlist(dic[key]["hgvs"]),\
                uniqlist(dic[key]["tag"]),\
                uniqlist(dic[key]["disease"]),\
                uniqlist(dic[key]["pmid"])
    return ",".join(mut_type), ",".join(mut_class), ",".join(mut_disease), ",".join(mut_pubmed)


def pileup2info(base):

    "yeild (alt, RDF, RDR, ADF, ADR)"

    header=re.compile("\^\S")
    insert=re.compile(r"\+[0-9]+[ACGTNacgtn]+")
    delere=re.compile(r"\-[0-9]+[ACGTNacgtn]+")
    #indelnum=re.compile(r"[\+\-][0-9]+")
    #remove read start mark
    base=header.sub("",base)
    #remove read end mark
    base=base.replace("$","")
    #Cal insert
    insertcounter=Counter(insert.findall(base))
    insertkey=set([x.upper() for x in insertcounter.keys()])
    base=insert.sub("", base)
    #Cal dele
    delecounter=Counter(delere.findall(base))
    delekey=set([x.upper() for x in delecounter.keys()])
    base=delere.sub("", base)
    #Cal base and snp
    basecounter=Counter(base)
    #chrm, pos, ref, alt, depth, forward, reverse
    rf=basecounter.get(".",0)
    rr=basecounter.get(",",0)
    hasmut=False
    for n in insertkey:
        hasmut=True
        #sys.stderr.write(n+"\n")
        yield(n, rf, rr, insertcounter.get(n, 0), insertcounter.get(n.lower(), 0))
    for n in delekey:
        hasmut=True
        #sys.stderr.write(n+"\n")
        yield(n, rf, rr, delecounter.get(n, 0), delecounter.get(n.lower(), 0))
    for n in "AGCT":
        nf=basecounter.get(n, 0)
        nr=basecounter.get(n.lower(), 0)
        if nf!=0 or nr!=0:
            hasmut=True
            #sys.stderr.write(n+"\n")
            yield(n, rf, rr, nf, nr)
    if not hasmut:
        yield(".",rf, rr, 0, 0)

def pileupcontrol(pileup):
    rawlis = [x for x in pileup2info(pileup)]
    def findmaxmut(info):
    #    alt, rdf, rdr, adf, adr = infos

    #for i in rawlis:

        alt, rdf, rdr, adf, adr = info
        AD = adf + adr
        RD = rdf + rdr
        if (AD+RD) != 0:
            AF = float(AD)/(AD+RD)
        else:
            AF = 0
        #if AD <= 10:
        if AD != 0:
            strand = float(adf)/AD
        else:
            strand = 0
        return AF, AD, strand
    lis = [findmaxmut(x) for x in rawlis]
    AF, AD, strand = sorted(lis,key=lambda x:x[0],reverse=True)[0]
    if AD <= 10:
        return "PASS"
    elif 0.1 <= strand <= 0.9:
        return "PASS"
    else:
        return "StrandFailed"

def parseintervar(intervar):
    dic = Ddict()
    with open(intervar) as fin:
        for i in fin:
            i = i.rstrip("\n")
            t = i.split("\t")
            if i.startswith("#"):
                continue
            chrom, start, end, ref, alt = t[:5]
            clinvar_raw = t[12]
            inter_raw = t[13]
            MIM, pheno_MIM, OrphaNumber, Orpha = t[-4:]
            clinvar = clinvar_raw.split()[-1]
            inter = inter_raw.split("PVS1")[0].split(": ")[-1]
            ACMG  = "PVS1" + inter_raw.split("PVS1")[-1]

            key = ",".join([chrom, start, ref, alt])
            dic[key]["intervar"] = inter
            dic[key]["ACMG"] = ACMG
            dic[key]["clinvar"] = clinvar
            dic[key]["MIM"] = MIM
            dic[key]["pheno_MIM"] = pheno_MIM
            dic[key]["OrphaNumber"] = OrphaNumber
            dic[key]["Orpha"] = Orpha

    return dic

def hcvcf2dic(vcf_file, DPfilter=0):
    """
    key: (chrom, pos, ref, alt) # vcf file raw value
    key2: AF, AD, DP
    """
    dic = Ddict()
    with open(vcf_file) as fin:
        for record in vcf.Reader(fin):
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            try:
                dp = record.samples[0]["DP"]
            except:
                dp = 0
            if dp < DPfilter:
                continue
            for altidx in range(len(record.ALT)):
                alt = str(record.ALT[altidx])
                try:
                    addp = record.samples[0]["AD"][altidx+1]
                except:
                    addp = 0
                if dp != 0:
                    af = round(float(addp)*100 / dp, 2)
                else:
                    af = 0
                key = (chrom, pos, ref, alt)
                dic[key]["DP"] = dp
                dic[key]["AD"] = addp
                dic[key]["AF"] = af
    return dic

def cut(str1, str2):
	lis1	= list(str1)
	lis2	= list(str2)
	len1	= len(lis1)
	len2	= len(lis2)

	##ins
	if len1 < len2:
		preStr2		= "".join(lis2[:len1])
		postStr2	= "".join(lis2[0-len1:])
		if preStr2 == str1:
			Str1 = ""
			Str2 = "".join(lis2[len1:])
		elif len1 > 1:
			s = 0
			for i in range(len1):
				if lis1[i] == lis2[i]:
					s += 1
				else:
					break
			Str1 = "".join(lis1[s:])
			Str2 = "".join(lis2[s:])

		else:
			Str1 = str1
			Str2 = str2
	##del
	elif len1 > len2:
		preStr1         = "".join(lis1[:len2])
		postStr1	= "".join(lis1[0-len2:])
		if preStr1 == str2:
			Str2 = ""
			Str1 = "".join(lis1[len2:])
		elif len2 > 1:
			s = 0
			for i in range(len2):
				if lis1[i] == lis2[i]:
					s += 1
				else:
					break
			Str1 = "".join(lis1[s:])
			Str2 = "".join(lis2[s:])

		else:
			Str1 = str1
			Str2 = str2
	else:
		if str1 == str2:
			Str1, Str2 = "", ""
		else:
			Str1, Str2 = "", ""
			for n, m in enumerate(str1):
				if m != str2[n]:
					Str2 += str2[n]
					Str1 += m

#			Str1 = str1
#			Str2 = str2
	return Str1, Str2
