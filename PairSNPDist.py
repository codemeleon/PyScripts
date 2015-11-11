from Bio import SeqIO
import argparse
import pandas as pd
import scipy as sp

def options( ):
    parser = argparse.ArgumentParser(description="Calculates Paired SNPs Distribution. Please change the code\
                                     based on your requirement first")
    parser.add_argument("--alnfile",dest="alnfile",help="Fasta Alignemnt File",type=argparse.FileType("r"))
    parser.add_argument("--groupfile",dest="groupfile",help="Tabbed group file",type=argparse.FileType("r"))
    parser.add_argument("--ignoregap",dest="ignoregap",help = "Remove gap from the analysis",type=bool, default=False)
    return parser.parse_args( )


if __name__=='__main__':
    args = options( )
    sequences = {}
    distribution = {}
    for rec in SeqIO.parse(args.alnfile,"fasta"):
        sequences[rec.id] = list(str(rec.seq).upper( ))#Change the code for suited for your need
    sequences = pd.DataFrame.from_dict(sequences)
    groups = pd.DataFrame.from_csv(args.groupfile.read( ))
    sequence1 = sequences[[group[]]]
    sequence2 = sequences[[group[]]]
    # comparion time

    #To data frame

    # # Do t-test on fly and produce the files only for with significant
    # difference

    # Add to ignore gaps and see the difference in the result

    # Check for the consitance of the gap. begining or at the end of the
    # alignemnt

    #
    if significat:
        pd.DataFrame.from_dict( ).to_csv( )

