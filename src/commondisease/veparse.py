#!/usr/bin/env python
import sys

#import hgvscDup2Ins
#All Extra terms list
Extra_terms = [
    "IMPACT", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE", "STRAND", "ENSP",
    "SWISSPROT", "TREMBL", "UNIPARC", "HGVSc", "HGVSp", "HGVS_OFFSET", "SIFT",
    "PolyPhen", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE", "CELL_TYPE", "CANONICAL", "CCDS", "INTRON", "EXON",
    "DOMAINS", "DISTANCE", "IND", "ZYG", "SV", "FREQS", "GMAF", "AFR_MAF",
    "AMR_MAF", "ASN_MAF", "EUR_MAF", "EAS_MAF", "SAS_MAF", "AA_MAF", "EA_MAF",
    "CLIN_SIG", "BIOTYPE", "TSL", "PUBMED", "SOMATIC", "PHENO", "ALLELE_NUM",
    "MINIMISED", "PICK", "REFSEQ_MATCH", "ExAC_MAF", "ExAC_Adj_MAF",
    "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", "ExAC_FIN_MAF",
    "ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF"
]

#Amino acids corresponding table
AA_dic = {
    "A": "Ala",
    "L": "Leu",
    "R": "Arg",
    "K": "Lys",
    "N": "Asn",
    "M": "Met",
    "D": "Asp",
    "F": "Phe",
    "C": "Cys",
    "P": "Pro",
    "Q": "Gln",
    "S": "Ser",
    "E": "Glu",
    "T": "Thr",
    "G": "Gly",
    "W": "Trp",
    "H": "His",
    "Y": "Tyr",
    "I": "Ile",
    "V": "Val",
    "*": "Ter"
}

dic_AA = {
    "Ala": "A",
    "Leu": "L",
    "Arg": "R",
    "Lys": "K",
    "Asn": "N",
    "Met": "M",
    "Asp": "D",
    "Phe": "F",
    "Cys": "C",
    "Pro": "P",
    "Gln": "Q",
    "Ser": "S",
    "Glu": "E",
    "Thr": "T",
    "Gly": "G",
    "Trp": "W",
    "His": "H",
    "Tyr": "Y",
    "Ile": "I",
    "Val": "V",
    "Ter": "*"
}
#Parse Extra, Return a dict
comple = {"A": "T", "T": "A", "C": "G", "G": "C"}


def hgvsc(hgvs):
    if "ins" in hgvs and "del" in hgvs:
        return hgvs.split("del")[0] + "delins" + hgvs.split("ins")[-1]
    elif "del" in hgvs:
        return hgvs.split("del")[0] + "del"
    elif "dup" in hgvs:
        return hgvs.split("dup")[0] + "dup"
    else:
        return hgvs


def hgvsp(hgvs):
    for key, val in dic_AA.items():
        hgvs = hgvs.replace(key, val)
    return hgvs


def Extra_Parse(Extra_str):
    """
    Extra_str is a string, vep file line_list[-1]
    """
    Extra_dic = {}
    for Extra in Extra_terms:
        Extra_dic[Extra] = "."
    Extra_list = Extra_str.split(";")
    for extra in Extra_list:
        key, val = extra.split("=")[0], extra.split("=")[1]
        Extra_dic[key] = val
    return Extra_dic


def Exist_Parse(Exist_str):
    """
    Exist_str is a string, vep file line_list[-2]; return dbsnp and cosmic
    """
    Exist_dic = {"dnSNP": ".", "COSMIC": "."}
    Exist_list = Exist_str.split(",")
    dnSNP_lis = []
    COSMIC_lis = []
    for i in Exist_list:
        if i.startswith("rs"):
            dnSNP_lis.append(i)
        elif i.startswith("COSM"):
            COSMIC_lis.append(i)
        else:
            pass
    if dnSNP_lis:
        Exist_dic["dnSNP"] = "|".join(dnSNP_lis)
    if COSMIC_lis:
        Exist_dic["COSMIC"] = "|".join(COSMIC_lis)
    return Exist_dic


def exac_dic(exac):
    dic = {}
    if exac == ".":
        return dic
    else:
        exac_t = exac.split(",")
        for sss in exac_t:
            t = sss.split(":")
            dic[t[0]] = t[1]
        return dic


Bases = [x for x in "ATCG"]
##Parse vep file

head_lis = [
    "##AA.stat", "Gene", "Gene.ID", "AA.Change", "Chr.start", "Chr.end", "Ref",
    "Alt", "Func", "dbSNP", "COSMIC", "1000g", "1000gEAS", "SIFT.score",
    "SIFT.pred", "PolyPhen.score", "PolyPhen.pred", "CLIN_SIG", "CANONICAL",
    "HGVS"
]

print ",".join(head_lis)


def single2triple(aa, aadic):
    lis = list(aa)
    return "".join([aadic[x] for x in lis])


def paser_vep(vepout):
    vep = open(vepout).readlines()
    mut_dic = {}
    key_lis = []
    for i in vep:
        i = i.rstrip()
        if i.startswith("#"):
            continue
        l = i.split()

        chr = "chr" + l[0].split("_")[0]
        start = l[1].split("-")[0].split(":")[-1]
        end = l[1].split("-")[-1].split(":")[-1]
        chr_start = ":".join([chr, start])
        chr_end = ":".join([chr, end])

        Ref = l[0].split("_")[-1].split("/")[0]
        Alt = l[2]
        Feature = l[4].split(":")[0].replace("-", ".")
        Consequence = l[6].replace(",", ";")
        Protein_pos = l[9].replace("-", ".")
        CDS_position = l[8].replace("-", ".")
        Amino_acids = l[10].replace("-", ".")
        if "/" in Amino_acids:
            AA_stat = "YES"

        else:
            AA_stat = "NO"

        Exist = l[-2]
        dbSNP = Exist_Parse(Exist)["dnSNP"]
        COSMIC = Exist_Parse(Exist)["COSMIC"]

        Extra = l[-1]
        dic = Extra_Parse(Extra)
        IMPACT = dic["IMPACT"]
        VARIANT_CLASS = dic["VARIANT_CLASS"]
        SYMBOL = dic["SYMBOL"]
        CANONICAL = dic["CANONICAL"]
        ENSP = dic["ENSP"]
        SIFT = dic["SIFT"]
        PolyPhen = dic["PolyPhen"]
        STRAND = dic["STRAND"]

        ExAC_MAF_raw = dic["ExAC_MAF"]
        ExAC_EAS_MAF_raw = dic["ExAC_EAS_MAF"]
        if not ExAC_MAF_raw:
            ExAC_MAF_raw = "."
        if not ExAC_EAS_MAF_raw:
            ExAC_EAS_MAF_raw = "."

        SIFT_pred_r = SIFT.split("(")[0]
        PolyPhen_pred_r = PolyPhen.split("(")[0]
        if "(" in SIFT:
            SIFT_score = SIFT.split("(")[-1][:-1]
        else:
            SIFT_score = "."
        if "(" in PolyPhen:
            PolyPhen_score = PolyPhen.split("(")[-1][:-1]
        else:
            PolyPhen_score = "."
        if SIFT_pred_r == "deleterious":
            SIFT_pred = "D"
        elif SIFT_pred_r == "tolerated":
            SIFT_pred = "T"
        elif SIFT_pred_r == "deleterious_low_confidence":
            SIFT_pred = "DL"
        elif SIFT_pred_r == "tolerated_low_confidence":
            SIFT_pred = "TL"
        else:
            SIFT_pred = "."

        if PolyPhen_pred_r == "benign":
            PolyPhen_pred = "B"
        elif PolyPhen_pred_r == "probably_damaging":
            PolyPhen_pred = "D"
        elif PolyPhen_pred_r == "possibly_damaging":
            PolyPhen_pred = "P"
        else:
            PolyPhen_pred = "."

        EXON_raw = dic["EXON"]
        INTRON_raw = dic["INTRON"]
        if EXON_raw == ".":
            EXON = "."
        else:
            EXON = "exon" + EXON_raw.split("/")[0] + "(" + EXON_raw.split(
                "/")[1] + ")"
        if INTRON_raw == ".":
            INTRON = "."
        else:
            INTRON = "intron" + INTRON_raw.split(
                "/")[0] + "(" + INTRON_raw.split("/")[1] + ")"

        CLIN_SIG = dic["CLIN_SIG"].replace(",", ";")

        HGVSc = dic["HGVSc"]
        if "synonymous_variant" in Consequence:
            HGVSp = ENSP + ":p." + single2triple(Amino_acids,
                                                 AA_dic) + Protein_pos + "="
        else:
            HGVSp = dic["HGVSp"]
            HGVSp_s = HGVSp.split(":")[-1]
            if not HGVSp_s.startswith("p."):
                HGVSp = "."

        if HGVSp != ".":
            HGVS = Feature + "(" + SYMBOL + "):" + hgvsc(
                HGVSc.split(":")[-1]) + "(" + hgvsp(HGVSp.split(":")[-1]) + ")"
        elif HGVSc != ".":
            HGVS = Feature + "(" + SYMBOL + "):" + hgvsc(HGVSc.split(":")[-1])
        else:
            HGVS = "."

        GMAF = "."
        EAS_MAF = "."

        if EXON != ".":
            Gene_id = ":".join([SYMBOL, Feature, EXON])
        elif INTRON != ".":
            Gene_id = ":".join([SYMBOL, Feature, INTRON])
        elif SYMBOL != ".":
            Gene_id = ":".join([SYMBOL, Feature])
        elif SYMBOL == "." and Feature != ".":
            Gene_id = Feature
        else:
            Gene_id = "."

        ##AAChange(SNV; InDel; framsheft)
        #Variation Type
        Type = ""
        if Ref in Bases and Alt in Bases:
            Type = "SNV"
        elif "frameshift_variant" in Consequence:
            Type = "frameshift"
        else:
            Type = "InDel"

        # AA_C, AA_P
        AA_C, AA_P = "", ""

        AAChange = "."
        key = ",".join([chr_start, Ref, Alt])

        mut_lis = [
            AA_stat, SYMBOL, Gene_id,
            AAChange.replace("*", "X"), chr_start, chr_end, Ref, Alt,
            Consequence, dbSNP, COSMIC, GMAF, EAS_MAF, SIFT_score, SIFT_pred,
            PolyPhen_score, PolyPhen_pred, CLIN_SIG, CANONICAL, HGVS
        ]
        key = ",".join([Gene_id, chr_start, chr_end, Ref, Alt])
        mut_dic[key] = mut_lis
        key_lis.append(key)
    return mut_dic, key_lis


if __name__ == "__main__":

    Usage = """
    Usage:   py shift1.vep

    Example: py  test.shift1.vep
    """

    if len(sys.argv) < 2:
        print Usage
        exit(1)
    shift1 = sys.argv[1]

    dic_shift1 = paser_vep(shift1)[0]
    lis = paser_vep(shift1)[1]

    for key in lis:
        tmp = dic_shift1[key]

        print ",".join(tmp)