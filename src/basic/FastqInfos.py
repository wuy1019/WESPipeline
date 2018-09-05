#!/usr/bin/env python
#-*- coding:utf-8 -*-
#Title: Fastqinfos final version
#Author: Mao Yibo
#Email: yibo.mao@geneseeq.com
#version: 0.4
from __future__ import print_function

import os,sys,subprocess

def basecount(seqstring):
    GC=seqstring.count("G")+seqstring.count("C")
    N=seqstring.count("N")
    return len(seqstring),GC,N

def qualitycount(qualstring):
    #chr(20+33)='5' chr(30+33)='?'
    q20=0
    q30=0
    for n in qualstring:
        if n>='?':
            q20+=1
            q30+=1
        elif n>='5':
            q20+=1
    return q20,q30
    
def fastqreader(fastqinpath): 
    if fastqinpath=='-':
        filein=sys.stdin
    else:
        if fastqinpath[-9:]==".fastq.gz" or fastqinpath[-6:]==".fq.gz":
            gzfastqin=subprocess.Popen(('/GPFS01/softwares/pigz-2.3.4/unpigz','-c',fastqinpath),stdout=subprocess.PIPE)
            filein=gzfastqin.stdout
	else:
	    filein=open(fastqinpath)
        #elif fastqinpath[-6:]==".fastq" or fastqinpath[-3:]==".fq":
        #    gzfastqin=subprocess.Popen(('cat',fastqinpath),stdout=subprocess.PIPE)
        #    filein=gzfastqin.stdout
        #else:
        #    raw_input("Input Fastq file format is ERROR! %s" %fastqinpath)
        #    sys.exit()
    count=0
    for line in filein:
        if line!='':
            count+=1
            #get the second line of a read (seq line)
            seqlen,gccount,ncount=basecount(filein.next().strip())
            filein.next()
            #get the fourth line of read (qual line)
            q20,q30=qualitycount(filein.next().strip())
            yield seqlen,gccount,ncount,q20,q30
        else:
            break
    filein.close()

def fastqinfo(fastqfile,id=None):
    if id==None:
        filename=os.path.basename(fastqfile)
    else:
        filename=id
    readnum=0
    maxlen=0
    basenum=0
    GCcount=0
    Ncount=0
    Q20=0
    Q30=0
    fastqin=fastqreader(fastqfile)
    for seqlen,gccount,ncount,q20,q30 in fastqin:
        readnum+=1
        maxlen=max(maxlen,seqlen)
        basenum+=seqlen
        GCcount+=gccount
        Ncount+=ncount
        Q20+=q20
        Q30+=q30
    print('#SAMPLE,READS,LEGTH,BASES,GC,Q20,Q30,PPM')
    print(','.join([str(x) for x in [filename,readnum,maxlen,max(round(basenum/1000/1000.,2),0.01),round(GCcount*100./basenum,2),round(Q20*100./basenum,2),round(Q30*100./basenum,2),round(Ncount*1000.*1000./basenum,2)]]))
    
if __name__=="__main__":
    if len(sys.argv)==2:
        assert os.path.exists(sys.argv[1]),"File %s does not exist.\nUsage:\npython script.py filepath\npython - sampleid" %sys.argv[1] 
        fastqinfo(sys.argv[1])
    elif len(sys.argv)==3:
        assert sys.argv[1]=='-',"Input information ERROR ' %s %s '.\nUsage:\npython script.py filepath\npython - sampleid" %(sys.argv[1],sys.argv[2])
        fastqinfo(*sys.argv[1:])
    else:
        print("Usage:\npython script.py filepath\npython - sampleid")
    
    
    
    
    
