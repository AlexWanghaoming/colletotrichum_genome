from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Calculate RIP index across sequence')
parser.add_argument('-fasta', metavar='fasta', type=str,
                    help='fasta file')
parser.add_argument('-windowsize', metavar='window size', type=int,
                    help='size of sliding windows')
parser.add_argument('-stepsize', metavar='step size', type=int,
                    help='overlap between windows')

args = parser.parse_args()

multifastafile  = args.fasta
records = list(SeqIO.parse((multifastafile),"fasta"))

def Count_TpA_and_ApT(sequence):
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

        #yield start, end, function(sequence[start:end])
        print(sequence.id,start, end, function(sequence[start:end]))
###http://book.biologistsguide2computing.com/en/stable/data-analysis.html#creating-a-sliding-window-gc-content-function


for sequence in records:
    sliding_window_analysis(sequence, Count_TpA_and_ApT, window_size=args.windowsize, step_size=args.stepsize)
