#!/usr/bin/env python
#-*- coding:utf-8 -*-
#Title: Fastq info and filter
#Author: Mao Yibo(yibo.mao@geneseeq.com)
#version: 0.5
from __future__ import print_function

import os
import sys
import csv
import argparse
from BioIO import *

def hamming_distance(seq,ref):
    #print(seq1,seq2)
    hmdistance=0
    for n in xrange(len(ref)):
        if seq[n]!=ref[n]:
            hmdistance+=1
    #print(hmdistance)s
    return hmdistance

def agct_seq(num):
    if num>1:
        sub_seq=agct_seq(num-1)
        for sseq in sub_seq:
            for n in "AGCT":
                yield sseq+n
    else:
        for n in "AGCT":
            yield n

seqs=[x for x in agct_seq(3)]
#print(seqs)
INSERTS={}
INSERTS[5]=["%sCT" %x for x in seqs]
INSERTS[6]=["%sACT" %x for x in seqs]
INSERTS[7]=["%sGACT" %x for x in seqs]
INSERTS[8]=["%sTGACT" %x for x in seqs]
#print(INSERTS)

i7_list		=	['CAGAGAAGAA', 'CTGTCCTCAA', 'ACCACCATAA', 'AAGCCTTAGA', 'TATTGTGCGA', 'CAATCAGTGA', 'ATCGCAAGAA', 'TCCTGCAGAA', 'CTCTGACGAA', 'TTGAAGCGAA', 'CACAACCGAA', 'CCACGTCGAA', 'TACGCCTGAA', 'CCTCCAACAA', 'AACGAGACAA', 'TTAGCCACAA', 'CCAGATACAA', 'ACTTGAGCAA', 'CTAACTGCAA', 'TTCACACCAA', 'TCAATGCCAA', 'CACGTATCAA', 'CTGCCAGTAA', 'TCGCGCGTAA', 'CCTATCGTAA', 'CTGTTAGAGA', 'TTGCGTCAGA', 'TCGAACTAGA', 'ATCATTCCGA', 'ACATACGTGA', 'CCGTTCTTGA', 'AACTCATGCA', 'ACGATTACCA', 'ACTCCGTCCA', 'TACCGCAGTT', 'ATGTCACCGT', 'TTCATGTAGC', 'ACCAGAAGAA', 'TCACGAAGAA', 'CCTTGAAGAA', 'CCAACAAGAA', 'TTGCCAAGAA', 'CACTCAAGAA', 'ACGAAGAGAA', 'TCTGAGAGAA', 'CCATAGAGAA', 'AACACGAGAA', 'CTTGCGAGAA', 'ACTCCGAGAA', 'TCGTCGAGAA', 'CAACTGAGAA', 'TACCACAGAA', 'CTCACCAGAA', 'CAAGCCAGAA', 'CATCTAAGAG', 'TAGAAGAGAG', 'CCTACGAGAG', 'ATGCCGAGAG', 'ACGCGAGATG', 'CTGCTGCATG', 'TCTCTTCATG', 'CCAGCCTATG', 'CAGATAAGCG', 'CTGAAGAGCG', 'ATCATCAGCG', 'TATCTTAGCG', 'CATCGCATCG', 'TCAATCATCG', 'ATTCTCATCG', 'TAATACGTCG', 'CTTGAGACGC', 'AATTGCGCGC', 'TAACCTGCGC', 'TCAAGACCGC', 'ACTACAATGC', 'TTCAAGATGC', 'TCTCCGATGC', 'TAGATGATGC', 'CATCCGTGCT', 'ATGTCGTGCT', 'CCTTCAACCT', 'CTATTCACCT', 'ATGCACGCCT', 'TTAGCCGCCT', 'ACGCCGATCT', 'ACCATGATCT', 'TCAGCTAGTT', 'CCGCGAACTT', 'ACACCGTTCG', 'ACCTCCTTCG', 'ACTTCGAGTG', 'AACGCCAGTG', 'ATTACACGTC', 'AATAGCCGTC', 'TTAGCATGTC', 'TTGAAGTGTC', 'TAGCCGTTCC', 'CAAGATTGTC', 'CAGTGTGCTC', 'CTCGAACCTC', 'TCTTCTCCTT', 'CACGCCTCTT']
i5_4k_list	=	['TGGACTCTAT', 'CGGTATCGGA', 'CAGACAGGTG', 'GGAGTTGTTA', 'AACTCAACGA', 'AAGCAGGCAA', 'TGGAGAACTT', 'TGAGCGGTTG', 'TGAGGTGATT', 'TAGCTAGATG', 'TAACATGTGG', 'CTACACTTGA', 'CGCTCGGTAA', 'CGGTGTTAGT', 'CGGCTGAAGA', 'CGGAGACTAG', 'CGATAGTCTA', 'CGAAGGAATG', 'CATGTCACTA', 'CACACGCTAG', 'CAGAAGTCGT', 'CAATCGGCTG', 'CAACGATCAA', 'GCACTTCATT', 'GGTAAGGTTA', 'GGCTAACAAT', 'GGACACGGAA', 'GGAGTATCAA', 'GGAACAATTG', 'GACATACAAG', 'GAGTCGCTGT', 'GAGCAGCGAT', 'AGAGTGGAAG', 'AGAGATGAGA', 'TTCTGAGAGG', 'TTCTGAACTG', 'TTCTATGGAG', 'TTCTACGATG', 'TTCTACACGT', 'TTCGTTGAGA', 'TTCGTCGCTA', 'TTCGTCATGT', 'TTCGTGTCGA', 'TTCGTACTGA', 'TTCGTAGCAT', 'TTCGTAGGTG', 'TTCGCTTAGG', 'TTCGCTGTTA', 'TTCGCGGTAT', 'TTCGCGACTT', 'TTCGCACTAG', 'TTCGCAAGGT', 'TTCGGTAGAA', 'TTCGGCTCAG', 'GTCTTCGCAG', 'GTCTGTAGTG', 'GTCTAGTCAT', 'GTCGATGTAT', 'GGTCGTGATT', 'GGTCATGCGA', 'GGTGCATTAA', 'GGTATCACAT', 'GATGTTACGA', 'GATACAGAGG', 'GACTGGAGGA', 'GACTATGGTT', 'GGAACATGTA', 'GGAAGCACTG', 'GACTAAGAAG', 'GACGTTGGTG', 'ATAGCGTCGG', 'ACGTTCGAAT', 'ACGTTAACAG', 'ACGTCTTGTA', 'ACGCTGGAGT', 'ACGCGTAGGT', 'ACGCGAGTAG', 'ACGATCACTT', 'CTGCAAGTGG', 'CTGGCAGCTA', 'CTATTCGCTG', 'CTATTACGAG', 'CTATGTCTTG', 'CTATGTAAGG', 'CTAGTGACTT', 'CTAGCAGTAA', 'CGTACAAGAT', 'CGCTTCGTTG', 'AAGTACACAG', 'AAGTAGCGGT', 'AAGCGAACGT', 'AAGCAGATAG', 'ACGTGTACGG', 'ACATCTAGAG', 'ACACACAGGA', 'AGCTTACGTG', 'CTATCTGTAG', 'CTAGGCATAT', 'CGTGTTATTG', 'CGTGGTCTGT', 'CGCAATATAG', 'CGGTTCTGAT']
i5_nova_list	=	['ATAGAGTCCA', 'TCCGATACCG', 'CACCTGTCTG', 'TAACAACTCC', 'TCGTTGAGTT', 'TTGCCTGCTT', 'AAGTTCTCCA', 'CAACCGCTCA', 'AATCACCTCA', 'CATCTAGCTA', 'CCACATGTTA', 'TCAAGTGTAG', 'TTACCGAGCG', 'ACTAACACCG', 'TCTTCAGCCG', 'CTAGTCTCCG', 'TAGACTATCG', 'CATTCCTTCG', 'TAGTGACATG', 'CTAGCGTGTG', 'ACGACTTCTG', 'CAGCCGATTG', 'TTGATCGTTG', 'AATGAAGTGC', 'TAACCTTACC', 'ATTGTTAGCC', 'TTCCGTGTCC', 'TTGATACTCC', 'CAATTGTTCC', 'CTTGTATGTC', 'ACAGCGACTC', 'ATCGCTGCTC', 'CTTCCACTCT', 'TCTCATCTCT', 'CCTCTCAGAA', 'CAGTTCAGAA', 'CTCCATAGAA', 'CATCGTAGAA', 'ACGTGTAGAA', 'TCTCAACGAA', 'TAGCGACGAA', 'ACATGACGAA', 'TCGACACGAA', 'TCAGTACGAA', 'ATGCTACGAA', 'CACCTACGAA', 'CCTAAGCGAA', 'TAACAGCGAA', 'ATACCGCGAA', 'AAGTCGCGAA', 'CTAGTGCGAA', 'ACCTTGCGAA', 'TTCTACCGAA', 'CTGAGCCGAA', 'CTGCGAAGAC', 'CACTACAGAC', 'ATGACTAGAC', 'ATACATCGAC', 'AATCACGACC', 'TCGCATGACC', 'TTAATGCACC', 'ATGTGATACC', 'TCGTAACATC', 'CCTCTGTATC', 'TCCTCCAGTC', 'AACCATAGTC', 'TACATGTTCC', 'CAGTGCTTCC', 'CTTCTTAGTC', 'CACCAACGTC', 'CCGACGCTAT', 'ATTCGAACGT', 'CTGTTAACGT', 'TACAAGACGT', 'ACTCCAGCGT', 'ACCTACGCGT', 'CTACTCGCGT', 'AAGTGATCGT', 'CCACTTGCAG', 'TAGCTGCCAG', 'CAGCGAATAG', 'CTCGTAATAG', 'CAAGACATAG', 'CCTTACATAG', 'AAGTCACTAG', 'TTACTGCTAG', 'ATCTTGTACG', 'CAACGAAGCG', 'CTGTGTACTT', 'ACCGCTACTT', 'ACGTTCGCTT', 'CTATCTGCTT', 'CCGTACACGT', 'CTCTAGATGT', 'TCCTGTGTGT', 'CACGTAAGCT', 'CTACAGATAG', 'ATATGCCTAG', 'CAATAACACG', 'ACAGACCACG', 'CTATATTGCG', 'ATCAGAACCG']

def index8(inlist):
    list8 = [i[:8] for i in inlist]
    return list8
i7_list.extend(index8(i7_list))
i5_4k_list.extend(index8(i5_4k_list))
i5_nova_list.extend(index8(i5_nova_list))
#print(i7_list)
#print(i5_4k_list)
#print(i5_nova_list)
    

"""
i5_nova_list=['ATAGAGTCCA', 'TCCGATACCG', 'CACCTGTCTG', 'TAACAACTCC', 'TCGTTGAGTT', 
'TTGCCTGCTT', 'AAGTTCTCCA', 'CAACCGCTCA', 'AATCACCTCA', 'CATCTAGCTA', 
'CCACATGTTA', 'TCAAGTGTAG', 'TTACCGAGCG', 'ACTAACACCG', 'TCTTCAGCCG', 
'CTAGTCTCCG', 'TAGACTATCG', 'CATTCCTTCG', 'TAGTGACATG', 'CTAGCGTGTG', 
'ACGACTTCTG', 'CAGCCGATTG', 'TTGATCGTTG', 'AATGAAGTGC', 'TAACCTTACC', 
'ATTGTTAGCC', 'TTCCGTGTCC', 'TTGATACTCC', 'CAATTGTTCC', 'CTTGTATGTC', 
'ACAGCGACTC', 'ATCGCTGCTC', 'CTTCCACTCT', 'TCTCATCTCT', 'CCTCTCAGAA', 
'CAGTTCAGAA', 'CTCCATAGAA', 'CATCGTAGAA', 'ACGTGTAGAA', 'TCTCAACGAA', 
'TAGCGACGAA', 'ACATGACGAA', 'TCGACACGAA', 'TCAGTACGAA', 'ATGCTACGAA', 
'CACCTACGAA', 'CCTAAGCGAA', 'TAACAGCGAA', 'ATACCGCGAA', 'AAGTCGCGAA', 
'CTAGTGCGAA', 'ACCTTGCGAA', 'TTCTACCGAA', 'CTGAGCCGAA']
"""
def index_check(filepath):
    index_count={}
    import gzip
    with gzip.open(filepath) as filein:
        count=0
        for line in filein:
            filein.next()
            filein.next()
            filein.next()
            index=line.strip().split(":")[-1]
            if index not in index_count:
                index_count[index]=0
            index_count[index]+=1
            count+=1
            if count >200:
                break
    index=sorted(index_count.items(), key=lambda x:x[1])[-1][0]
    if "+" not in index:
        return False
    else:
        i7, i5 =index.split("+")
        if i7 in i7_list and (i5 in i5_4k_list or i5 in i5_nova_list):
            return True
        else:
            return False


def seq_compreverse(seq):
    """
    Get the complementary reverse sequence
    reviewed
    """
    pairedbase={'A':'T','G':'C','C':'G','T':'A','N':'N'}
    outreadseq=''
    for n in xrange(len(seq)):
        outreadseq=pairedbase[seq[n]]+outreadseq
    return outreadseq

def insert_check(insertseq, dual_index=True, mismatch=False):
    #print(insertseq)
    #The first 7nt read seq
    for n in range(8,4,-1):
        if insertseq[:n] in INSERTS[n]:
            return insertseq[:n],n,0
    if mismatch!=False:
        hmdisc=defaultdict(list)
        for i in INSERTS[8]:
            hmdisc[hamming_distance(insertseq[:8], i)].append(i)
        for i in INSERTS[7]:
            hmdisc[hamming_distance(insertseq[:7], i)].append(i)
        for i in INSERTS[6]:
            hmdisc[hamming_distance(insertseq[:6], i)].append(i)
        for i in INSERTS[5]:
            hmdisc[hamming_distance(insertseq[:5], i)].append(i)
        if len(hmdisc[1])==1:
            return hmdisc[1][0],len(hmdisc[1][0]),1
    return insertseq,len(insertseq),-1

def tail_check(tail, insert):
    """if the read len is less than 150, the seq is read through, and cut the insert length
    if the length is more or equal than 150, the tail isn't seqed or half seqed.
    """
    insert_cr=seq_compreverse(insert)
    length=len(insert)
    for n in xrange(length-2):
        if hamming_distance(insert_cr[:length-n], tail[n:length]) <=1:
            return length-n, tail[n:length]
    for n in xrange(length-2, length):
        if insert_cr[:length-n]==tail[n:length]:
            return length-n, tail[n:length]
    return 0, ""

def get_mid(header):
    insert=header.split(":")[-1]
    if "+" in insert:
        i7, i5=insert.split("+")
        return i7[-2:]+i5[-2:], True
    else:
        return insert[-4:], False

def header_modify(header, mid, read1_insert, read2_insert):
    #"@{}:{}:{}".format(seqinsert,midindex,read2[0][1:])
    header_pre, header_ext= header.split()
    if read1_insert<read2_insert:
        seqinsert="{}-{}:ab".format(read1_insert, read2_insert)
    elif read1_insert>read2_insert:
        seqinsert="{}-{}:ba".format(read2_insert, read1_insert)
    else:
        seqinsert="{}-{}:cc".format(read1_insert, read2_insert)
    return "@{}:{} 1{}".format(seqinsert, header_pre[1:], header_ext[1:]), "@{}:{} 2{}".format(seqinsert, header_pre[1:], header_ext[1:])

def sub_insert_cuter(read1path, read2path, rm):
    infocount=defaultdict(int)
    fastq_in=pefastq_reader(read1path, read2path)
    read1_sucess="{}_sucess.fastq.gz".format(read1path)
    read2_sucess="{}_sucess.fastq.gz".format(read2path)
    fastq_out=PeFastqWriter(read1_sucess, read2_sucess)
    for read1, read2 in fastq_in:
        read1[1]=read1[1][1:-1]
        read2[1]=read2[1][1:-1]
        read1[3]=read1[3][1:-1]
        read2[3]=read2[3][1:-1]
        infocount["Total_Reads"]+=1
        mid, dual_index=get_mid(read1[0])
        if "N" in mid:
            infocount["MID_Failed"]+=1
            continue
        read1_insert, read1_insert_len, read1_insert_mode=insert_check(read1[1][:8], dual_index)
        read2_insert, read2_insert_len, read2_insert_mode=insert_check(read2[1][:8], dual_index)
        if read1_insert_mode==-1 or read2_insert_mode==-1:
            infocount["Insert_Failed"]+=1
            if read1_insert_mode==-1:
                infocount["Read1_Insert_Failed"]+=1
            if read2_insert_mode==-1:
                infocount["Read2_Insert_Failed"]+=1
            continue
        read1_tail_len, read1_tail=tail_check(read1[1][-read2_insert_len:], read2_insert)
        read2_tail_len, read2_tail=tail_check(read2[1][-read1_insert_len:], read1_insert)
        if rm:
             read1[0], read2[0]=header_modify(read1[0], mid, read1_insert, read2_insert)
        read1_len=len(read1[1])
        read1[1]=read1[1][read1_insert_len:read1_len-read1_tail_len]
        read1[3]=read1[3][read1_insert_len:read1_len-read1_tail_len]
        read2_len=len(read2[1])
        read2[1]=read2[1][read2_insert_len:read2_len-read2_tail_len]
        read2[3]=read2[3][read2_insert_len:read2_len-read2_tail_len]
        infocount["Success"]+=1
        fastq_out.write(read1, read2)
    fastq_out.close()
    return read1_sucess, read2_sucess, infocount

def insert_cuter(read1path, read2path, read1sucess, read2sucess,rm=None, thread=8):
    if os.path.exists(read1sucess) and os.path.exists(read2sucess):
        sys.stderr.write("Result FIle Exists\n")
        return 0
    if index_check(read1path)!=True:
        sys.stderr.write("Not a ctDNA 2.0 V3 Adapter Sample\n")
        os.symlink(os.path.abspath(read1path), os.path.abspath(read1sucess))
        os.symlink(os.path.abspath(read2path), os.path.abspath(read2sucess))
        return 0
    else:
        sys.stderr.write("ctDNA 2.0 V3 Adapter Sample\n")
    if rm==None:
        rm=True
    else:
        rm=False
    fastq_folder=os.path.split(os.path.abspath(read1path))[0]
    temp_folder=temp_folder="{}/Split_{}".format(os.path.abspath("."), "".join(random.sample(string.ascii_letters,8)))
    thread=int(thread)
    split_pool=multiprocessing.Pool(2)
    read1_split=split_pool.apply_async(split_files, (read1path, 8000000, "{}/read1".format(temp_folder)))
    read2_split=split_pool.apply_async(split_files, (read2path, 8000000, "{}/read2".format(temp_folder)))
    split_pool.close()
    split_pool.join()
    read1_subpath, read1_subfiles=read1_split.get()
    read2_subpath, read2_subfiles=read2_split.get()
    sys.stderr.write("Sub File Adapter Cuter Started!\n")
    cuter_pool=multiprocessing.Pool(thread*3)
    result_list=[]
    for num, read1_path in enumerate(read1_subfiles):
        read2_path=read2_subfiles[num]
        result_list.append(cuter_pool.apply_async(sub_insert_cuter, (read1_path, read2_path, rm)) )
    cuter_pool.close()
    cuter_pool.join()
    read1list=[]
    read2list=[]
    #import ipdb
    #ipdb.set_trace()
    info_dict=defaultdict(int)
    for n in result_list:
        sub_read1, sub_read2, sub_info=n.get()
        read1list.append(sub_read1)
        read2list.append(sub_read2)
        for key, value in sub_info.iteritems():
            info_dict[key]+=value
    read1list.sort()
    read2list.sort()
    cat_file_pool=multiprocessing.Pool(thread)
    cat_file_pool.apply_async(gzip_cat_files, (read1sucess, read1list))
    cat_file_pool.apply_async(gzip_cat_files, (read2sucess, read2list))
    cat_file_pool.close()
    cat_file_pool.join()
    try:
        shutil.rmtree(temp_folder)
    except:
        pass
    info_dict["Failed"]=info_dict["Total_Reads"]-info_dict["Success"]
    for n in ["Total_Reads","Success","Failed","MID_Failed","Insert_Failed","Read1_Insert_Failed","Read2_Insert_Failed"]:
        info_dict.setdefault(n, 0)
    sys.stdout.write("#Total_Reads,Success,Failed,MID_Failed,Insert_Failed,Read1_Insert_Failed,Read2_Insert_Failed\n")
    sys.stdout.write("{Total_Reads},{Success},{Failed},{MID_Failed},{Insert_Failed},{Read1_Insert_Failed},{Read2_Insert_Failed}\n".format(**info_dict))
    
if __name__=="__main__":
    if len(sys.argv)>=5:
        insert_cuter(*sys.argv[1:])
    else:
        print("Usage:\npython script.py read1 read2 read1sucess read2sucess [remove_insertid]")
    

    
