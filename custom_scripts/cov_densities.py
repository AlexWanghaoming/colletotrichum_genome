#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 23:21:55 2020

@author: Pam
"""
def make_genomefile(fasta,outfilename):
    import os
    #Requires samtools
    cmd = "samtools faidx "+fasta
    cmd2 = '''awk '{print $1"\t"$2}' '''+fasta+".fai > "+outfilename
    os.system(cmd)
    os.system(cmd2)
    
def make_windowsfile(genomefile,window,sliding=None,prefix=None):
    from pybedtools import BedTool
    if prefix !=None:
    #If a prefix is supplied, write to output file
        outputfilename = prefix+"."+str(int(window/1000))+"k."+str(int(sliding/1000))+"k.windows.bed"
    else:
    #Else return a bedtools object
        outputfilename = None
        win = BedTool().window_maker(g=genomefile,w=window,s=sliding,output=outputfilename)
        return win

def import_genegff(genegff):
    import pandas as pd
    df = pd.read_csv(genegff,header=None, sep="\t",skiprows=1)
    df.columns=["chr","program","type","start","end","score","strand","frame","description"]
    df2= pd.DataFrame(df["description"].str.split(";",0).tolist(),columns=["ID","Name"])
    df2["ID"]=df2["ID"].str.replace(r"ID=","")
    allgenesdf = df.join(df2["ID"])
    return allgenesdf

def get_goi_loc(allgenesdf,genelist):
    #Make bedfile for location of specific genes of interest in a list of genes. Input is a text file 
    #of a list of genenames and a dataframe of all genes and their locations
    import pandas as pd
    from pybedtools import BedTool
    goi = pd.read_csv(genelist,names=["ID"])
    #Clean up suffixes of transcripts sometimes included
    goi["ID"]=goi["ID"].str.replace(".t1","")
    goidf = pd.merge(allgenesdf,goi,on=["ID"],how="inner")
    goidf2 = goidf[["chr","start","end"]]
    bed = BedTool.from_dataframe(goidf2)
    return bed

def genomewide_cov(genomefile,window,feature,sliding=None,prefix=None):
    #Feature can be a bedtools object from get_goi_loc or a gff or bed filename
    import pandas as pd
    from pybedtools import BedTool
    if sliding == None:
        sliding = window
    if isinstance(feature,BedTool):
        featurebed = feature
    if isinstance(feature,str):
        featurebed = BedTool(feature)
    win = make_windowsfile(genomefile,window,sliding=sliding)
    cov = win.coverage(featurebed)
    df = pd.read_csv(cov.fn,sep="\t",names=["chr","start","end","count","bases","len","cov"])
    df["mid"]=(df["start"]+df["end"])/2
    return df

