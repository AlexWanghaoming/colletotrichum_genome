#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd

parser = argparse.ArgumentParser(description="""Getfasta from gff using bedtools getfasta. Will merge features that are split.
Example usage to get CDS with two separators:
parse_gff_to_fasta2.py -gff test.gff -feature exon -num 1 -s orig_protein_id= -s \;orig_transcript -fasta ref.fa -outfile out.fasta""",
formatter_class=RawTextHelpFormatter)
parser.add_argument('-gff', metavar='gff', type=str,
                    help='gff file', required=True)
parser.add_argument('-feature', metavar='feature', type=str, help='feature type e.g. mRNA, exon', required=True)
parser.add_argument('-num', metavar='number', type=int,
                    help='field number in the description field of the gff')
parser.add_argument('-s',metavar="separators",action="append",type=str,help="separators default is ;")
parser.add_argument('-fasta', metavar='fasta', type=str,
                    help='reference fasta', required=True)
parser.add_argument('-outfile', metavar='outfile', type=str,
                    help='output filename', required=True)

args = parser.parse_args()
feature = args.feature
gff = args.gff
fieldnum = args.num
reffasta = args.fasta
outfile = args.outfile
seps = args.s

if seps is None:
    seps = [";"]


#from stackoverflow:https://stackoverflow.com/questions/4697006/python-split-string-by-list-of-separators
def split(txt, seps):
    default_sep = seps[0]

    # we skip seps[0] because that's the default seperator
    for sep in seps[1:]:
        txt = txt.replace(sep, default_sep)
    return [i.strip() for i in txt.split(default_sep)]

def import_gff_create_tmp(gff,feature,fieldnum):
    """Parse gff file so that we only have feature gff lines"""
    tmpout = open(gff+".tmp","w+")
    with open(gff) as fh:
        for lines in fh:
            lines = lines.strip()
            if lines.startswith("#"):
                pass
            else:
                line = lines.split("\t")
                if line[2]==feature:
                    field0 = split(line[8],seps)
                    field = field0[fieldnum]
                    l = lines+"\t"+field+"\n"
                    tmpout.write(l)
def get_featurename(gff):
    """Get names that will be used for fastaheaders and get start sites relative to start of document"""
    df = pd.read_csv(gff+".tmp",sep="\t",header=None)
    #seps = sep1+"|"+sep2
    #df2 = df[8].str.split(seps,expand=True)
    #df2 = df[8].str.split("orig_protein_id=|;orig_transcript_id",expand=True)
    df["ID"] = df[9]
    df["len"] = df[4]-df[3]+1
    df3 = df.groupby(["ID"]).min().reset_index()
    df3 = df3[["ID",3]]
    df3.columns = ["ID","CDSstart"]
    df4 = df.groupby(["ID"]).max().reset_index()
    df4 = df4[["ID",4]]
    df4.columns = ["ID","CDSend"]
    mdf = pd.merge(df,df3,on="ID")
    mdf = pd.merge(mdf,df4,on="ID")
    #mdf["plusfeaturestart"] = mdf[3]-mdf["CDSstart"]
    mdf["blockstart"] = mdf[3]-mdf["CDSstart"]
    #mdf["minusfeaturestart"] = mdf["CDSend"]-mdf[4]
    sdf = mdf.sort_values(["ID",3])
    sdf["bedstart"] = sdf[3]-1
    df5 = df.groupby(["ID"]).count().reset_index()
    df5 = df5[["ID",0]]
    df5.columns = ["ID","count"]
    sdf2 = pd.merge(df5,sdf,on="ID")
    sdf2["bedCDSstart"] = sdf2["CDSstart"]-1
    sdf2 = sdf2[["ID","count","bedCDSstart",0,6,"len","CDSstart","CDSend","blockstart","bedstart"]]
    sdf2 = sdf2.rename(columns={0:"seq",6:"strand"})
    #plusseqlist = sdf2["ID"][sdf2["strand"]=="+"].drop_duplicates().tolist()
    #minusseqlist = sdf2["ID"][sdf2["strand"]=="-"].drop_duplicates().tolist()
    seqlist = sdf2["ID"].drop_duplicates().tolist()
    return sdf2,seqlist


def sql_queries2(sdf2,seqlist):
    from sqlalchemy import create_engine
    from sqlalchemy.sql import text
    engine = create_engine('sqlite://', echo=False)

    sdf2.to_sql(name="gfftable",con=engine,if_exists="append")

    lengths = []

    for i in seqlist:
        s = text("SELECT len FROM gfftable WHERE ID=:fastaname")
        out = engine.execute(s,fastaname=i).fetchall()
        out2 = [str(x[0]) for x in out]
        out3 = (',').join(out2)
        out4 = (i,out3+",")
        lengths.append(out4)


    ldf = pd.DataFrame(lengths)
    ldf.columns = ["ID","block_lengths"]

    blockstarts = []

    for i in seqlist:
        s = text("SELECT blockstart FROM gfftable WHERE ID=:fastaname")
        out = engine.execute(s,fastaname=i).fetchall()
        out2 = [str(x[0]) for x in out]
        out3 = (',').join(out2)
        out4 = (i,out3+",")
        blockstarts.append(out4)


    bsdf = pd.DataFrame(blockstarts)
    bsdf.columns = ["ID","block_starts"]


    return ldf,bsdf

def merge_it_all(ldf,bsdf,sdf2,gff):
    """merge to get bed12 format"""
    sdf2["phase"] = 0
    sdf3 = sdf2[["seq","bedCDSstart","CDSend","ID","phase","strand","bedCDSstart","bedCDSstart","phase","count"]]
    sdf3 = sdf3.drop_duplicates()
    mergedf = pd.merge(sdf3,ldf,on="ID")
    mergedf = pd.merge(mergedf,bsdf,on="ID")
    mergedf.to_csv(gff+"tmp.bed12",sep="\t",index=None,header=None)


def get_fasta(reffasta,outfile,gff):
    import os
    cmd = ("bedtools getfasta -fi "+reffasta+" -bed "+gff+"tmp.bed12 -split -name -s > "+outfile+".tmp")
    cmd2 = ("fasta_formatter -i "+outfile+".tmp -o "+outfile+" -w 60")
    cmd3 = ("rm "+outfile+".tmp")
    cmd4 = ("rm "+gff+"tmp.bed12")
    os.system(cmd)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

import_gff_create_tmp(gff,feature,fieldnum)
sdf2,seqlist = get_featurename(gff)
ldf,bsdf = sql_queries2(sdf2,seqlist)
merge_it_all(ldf,bsdf,sdf2,gff)
get_fasta(reffasta,outfile,gff)
