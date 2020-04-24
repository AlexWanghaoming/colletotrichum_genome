
def find_sequence(motif,sequence,sequencename,strand,outputname,seqtype):
    import re
    with open(outputname,"a+") as fh:
        from re import finditer
        for match in finditer(motif,sequence,flags=re.IGNORECASE):
            start,end = match.span()
            start = start+1
            end = end
            out = sequencename+"\t"+seqtype+"\t"+seqtype+"\t"+str(start)+"\t"+str(end)+"\t0\t"+strand+"\t.\tNote=(TTAGGG)\n"
            fh.write(out)

def make_telomere_gff(fasta,outputname):
    import os
    from Bio import SeqIO
    records = list(SeqIO.parse((fasta),"fasta"))
    try:
        os.remove(outputname)
    except:
        pass
    for sequences in records:
        find_sequence("TTAGGG",str(sequences.seq),str(sequences.id),"+",outputname,"telrep")
        find_sequence("CCCTAA",str(sequences.seq),str(sequences.id),"-",outputname,"telrep")


def find_telomeres(window,genomefile,telomeres_gff,percentile,outfilename):
    from pybedtools import BedTool
    import pandas as pd
    import numpy as np
    win = BedTool().window_maker(g=genomefile,w=window)
    tel = BedTool(telomeres_gff)
    cov = win.coverage(tel)
    df = pd.read_csv(cov.fn,sep="\t",names=["chr","start","end","depth","n bp","len","cov"])
    df["per_bp"] = df["depth"]/df["len"]
    a = np.percentile(df["per_bp"],percentile)
    tel = df[df["per_bp"]>=a].reset_index()
    tel["num"] = "tel"+tel.index.map(str)
    tel2=tel[["chr","start","end","num"]]
    tel2.to_csv(outfilename,sep="\t",header=None,index=None)
    

