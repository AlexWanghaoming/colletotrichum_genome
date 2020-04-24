from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser(description='Calculate dinucleotide frequencies across sequence. Output: chrom 	start 	end 	AA 	AC 	AG 	AT 	CA 	CC 	CG 	CT 	GA 	GC 	GG 	GT 	TA 	TC 	TG 	TT')
parser.add_argument('-fasta', metavar='fasta', type=str,
                    help='fasta file')
parser.add_argument('-windowsize', metavar='window size', type=int,
                    help='size of sliding windows')
parser.add_argument('-stepsize', metavar='step size', type=int,
                    help='overlap between windows')
parser.add_argument('-output', metavar='output', type=str,
                    help='output csv file with dinucleotide frequencies')

args = parser.parse_args()

multifastafile  = args.fasta
records = list(SeqIO.parse((multifastafile),"fasta"))

def count_dinucleotides(sequence):
    '''Count TpA and ApT'''

    dinucfreq = []
    count = 0
    #AA
    for match in re.findall(r"A(?=A)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"A(?=C)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"A(?=G)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"A(?=T)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"C(?=A)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"C(?=C)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"C(?=G)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"C(?=T)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"G(?=A)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"G(?=C)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"G(?=G)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"G(?=T)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"T(?=A)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"T(?=C)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"T(?=G)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    count = 0
    for match in re.findall(r"T(?=T)",str(sequence.seq),re.IGNORECASE):
        count=count+1
    dinucfreq.append(count)

    return dinucfreq[:]


def sliding_window_analysis(sequence, function, window_size=args.windowsize, step_size=args.stepsize):

    """Return an iterator that yields (start, end, property) tuples.

    Where start and end are the indices used to slice the input list
    and property is the return value of the function given the sliced
    list.
    """

    for start in range(0, len(sequence), step_size):
        end = start + window_size
        if end > len(sequence):
            break

        yield start, end, function(sequence[start:end])

###http://book.biologistsguide2computing.com/en/stable/data-analysis.html#creating-a-sliding-window-gc-content-function

with open(args.output,"w") as fh:
    header = "seq_id,start,mid,end,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT\n"
    fh.write(header)
    for sequence in records:
        for start, end, freqs in sliding_window_analysis(sequence, count_dinucleotides, window_size=args.windowsize, step_size=args.stepsize):
            AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT = freqs
            mid = (int(start)+int(end))/2
            row = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(sequence.id,start,mid,end,AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT)
            fh.write(row)
