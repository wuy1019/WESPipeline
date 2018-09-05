#!/usr/bin/env python
#-*- coding:utf-8 -*-
#Title: CovDepInfos
#Author: Mao Yibo
#Email: yibo.mao@geneseeq.com
#version: 0.1
from __future__ import print_function

import os,sys,csv
from collections import defaultdict

def pause(string=""):
    raw_input("%s Any key:\n======================================" %string)

def qcinfoloader(filepath,sampleid):
    infodatadict={}
    with open(filepath) as filein:
        for line in filein:
            if line.startswith("#"):
                continue
            lineinfo=line.strip().split(",")
            infodatadict[lineinfo[0]]=float(lineinfo[3])
    sample_total_base=0
    for sid, readbase in infodatadict.iteritems():
        if sampleid==sid.split("_")[0]:
            sample_total_base+=readbase
    return sample_total_base

def basefileloader(basefilepath):
    """
    Receive a bedtools coverage -d result file and return the covery depth infors.
    """
    filename=os.path.splitext(os.path.basename(basefilepath))[0]
    sampleid=filename.split(".")[0]
    datatype=filename.split(".")[1:3]
    print("Analysising %s" %filename)
    with open(basefilepath) as filein:
        #statistics total info
        totaldepth=0
        totallength=0
        totalcovery=0
        depthdict=defaultdict(lambda:0)
        #statistics region info
        regioninfodict=defaultdict(lambda:[])
        regionkey=[]
        regionlenth=0
        regioncovery=0
        regiondepth=0
        #reading file
        for line in filein:
            if line[0]=="#":
                continue
            lineinfo=line.strip().split()
            regionid=lineinfo[:4]
            #check the id
            if regionid!=regionkey:
                if regionkey==[]:
                    regionkey=regionid
                else:
                    regionlenth=int(regionkey[2])-int(regionkey[1])
                    regioninfodict[tuple(regionkey)]=[regionlenth,regioncovery,regiondepth]
                    totallength+=regionlenth
                    totalcovery+=regioncovery
                    totaldepth+=regiondepth
                    #new region init
                    regionkey=regionid
                    regionlenth=0
                    regioncovery=0
                    regiondepth=0
            basedepth=int(lineinfo[-1])
            if basedepth!=0:
                regiondepth+=basedepth
                regioncovery+=1
                depthdict[basedepth]+=1
        #after the file read loop store the last region infor
        regionlenth=int(regionkey[2])-int(regionkey[1])
        regioninfodict[tuple(regionkey)]=[regionlenth,regioncovery,regiondepth]
        totallength+=regionlenth
        totalcovery+=regioncovery
        totaldepth+=regiondepth
    #base file info load finished
    #print(totallength,totalcovery,totaldepth)
    #pause()
    #cal covdepinfos
    #meandepth=round(totaldepth/float(totalcovery),2)
    meandepth=round(totaldepth/float(totallength),2)
    meandepth20=round(meandepth*.2,2)
    coverydict=defaultdict(lambda:0.0)
    for depth,depthcount in depthdict.items():
        if depth>=1:
            coverydict['1']+=depthcount
        if depth>=10:
            coverydict['10']+=depthcount
        if depth>=20:
            coverydict['20']+=depthcount
        if depth>=50:
            coverydict['50']+=depthcount
        if depth>=meandepth20:
            coverydict[meandepth20]+=depthcount
    return sampleid,datatype,totallength,totalcovery,totaldepth,meandepth,meandepth20,regioninfodict,coverydict

def chrcmp(chra,chrb):
    chra=chra[3:]
    chrb=chrb[3:]
    if chra.isdigit()==chrb.isdigit() and chra.isdigit():
        return cmp(int(chra),int(chrb))
    elif chra.isdigit()==chrb.isdigit() and not chra.isdigit():
        return cmp(chra,chrb)
    elif chra.isdigit() and not chrb.isdigit():
        return -1
    elif not chra.isdigit() and chrb.isdigit():
        return 1


def filtercovdep(regiondict,avedepth,savefilepath,coverage=0.7,avedepthpct=0.2):
    #AverageDeepth
    #CHR,ExonStart,ExonEnd,NAME,Reads,CoveragedSites,Range,CoveragePct,Deepth,Stat,StatCov,StatDph
    fileout=open(savefilepath,'w')
    fileout.write("#AverageDepth:%s\n" %avedepth)
    fileout.write("#CHR,ExonStart,ExonEnd,NAME,Reads,CoveragedSites,Range,CoveragePct,Depth,Stat,StatCov,StatDph\n")
    #bad fileinit
    badfilepath=savefilepath.replace(".csv",".bad.csv")
    badfileout=open(badfilepath,'w')
    badfileout.write("#AverageDepth:%s\n" %avedepth)
    badfileout.write("#CHR,ExonStart,ExonEnd,NAME,Reads,CoveragedSites,Range,CoveragePct,Depth,Stat,StatCov,StatDph\n")
    regionkeylist=regiondict.keys()
    regionkeylist=sorted(regionkeylist,lambda x,y:cmp(int(x[1]),int(y[1])))
    regionkeylist=sorted(regionkeylist,lambda x,y:chrcmp(x[0],y[0]))
    for regionkey in regionkeylist:
        regioninfo=regiondict[regionkey]
        goodcode="GOOD"
        covcode="-"
        depcode="-"
        regioncov=round(1.0*regioninfo[1]/regioninfo[0],2)
        if regioninfo[1]!=0:
            regiondph=round(1.0*regioninfo[2]/regioninfo[0],2)
        else:
            regiondph=0.0
        if regioncov<=coverage:
            covcode="LowCoverage"
            goodcode="BAD"
        if regiondph<=avedepthpct*avedepth:
            depcode="LowDepth"
            goodcode="BAD"
        outinfo=list(regionkey)+["-",regioninfo[1],regioninfo[0],regioncov,regiondph,goodcode,covcode,depcode]
        outinfo=[str(x) for x in outinfo]
        fileout.write(','.join(outinfo)+"\n")
        if goodcode=="BAD":
            badfileout.write(','.join(outinfo)+"\n")
    #print("AVE DEPTH %s" %avedepth )

def covdepinfos(sampleid,meandepth,coverydict,totallength,meandepth20,savefilepath):
    #sample,depth,1X,10X,20X,50X,20%X
    outstring="%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" %(sampleid,meandepth,coverydict['1']*100/totallength,coverydict['10']*100/totallength,coverydict['20']*100/totallength,coverydict['50']*100/totallength,coverydict[meandepth20]*100/totallength)
    #print(outstring)
    #open(savefilepath,'a').write(outstring)
    os.system("echo '%s' >>%s" %(outstring,savefilepath))

def baseontarget(sampleid,hqtotalbase,totaldepth,savefilepath):
    #need hqdata readlen * readnumber *2 equal total base number
    totalbase=round(totaldepth/1000/1000.,2)
    outstring="%s,%d,%d,%.2f" %(sampleid,hqtotalbase,totalbase,totalbase/hqtotalbase*100)
    #open(savefilepath,'a').write(outstring)
    os.system("echo '%s' >>%s" %(outstring,savefilepath))

def qcinfo(basefilepath,qcinfospath):
    #load infors from base file
    sampleid,datatype,totallength,totalcovery,totaldepth,meandepth,meandepth20,regioninfodict,coverydict=basefileloader(basefilepath)
    #Load qc file total base Cleaninfo.cav Million bases
    hqtotalbase=qcinfoloader(qcinfospath, sampleid)
    avedepth=float(totaldepth)/totalcovery
    print("Total depth:%s,total covery %s,ave depth %s" %(totaldepth,totalcovery,avedepth))
    #Saving data infor
    #print(datatype)
    if datatype==["dedup","raw"]:
        #Save FilterCovDep
        savefilepath="./TMP/%s.dedup.csv" %sampleid
        filtercovdep(regioninfodict,avedepth,savefilepath)
        #Save covdepinfos
        savefilepath="Results/Infos/CovDepDedupInfos.csv"
        covdepinfos(sampleid,meandepth,coverydict,totallength,meandepth20,savefilepath)
    elif datatype==["sorted","raw"]:
        #Save FilterCovDep
        savefilepath="./TMP/%s.original.csv" %sampleid
        filtercovdep(regioninfodict,avedepth,savefilepath)
        #Save covdepinfos
        savefilepath="Results/Infos/CovDepOriginalInfos.csv"
        covdepinfos(sampleid,meandepth,coverydict,totallength,meandepth20,savefilepath)
        #Save ontargetpct
        savefilepath="Results/Infos/OnTarget.raw.csv"
        baseontarget(sampleid,hqtotalbase,totaldepth,savefilepath)
    elif datatype==["sorted","add"]:
        #Save ontargetpct
        savefilepath="Results/Infos/OnTarget.add.csv"
        baseontarget(sampleid,hqtotalbase,totaldepth,savefilepath)

if __name__=="__main__":
    if len(sys.argv)==3:
        qcinfo(sys.argv[1],sys.argv[2])
    else:
        print("Useage:\n python script.py basefilepath qcinfospath")

