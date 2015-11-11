from Bio import SeqIO
import argparse

def option( ):
    opt = argparse.ArgumentParser( description="Split fasta file into sequence if file based on sequence names")
    opt.add_argument("--infasta",dest="infasta",help="Multiple fasta file",type=argparse.FileType(mode="r"))
    return opt.parse_args( )

if __name__== '__main__':
    opt = option( )
    for rec in SeqIO.parse(opt.dest,"fasta"):
        with open("%s.fasta"%rec.id,"w") as fout:
            fout.write(">%s\n%s\n"%(rec.id,rec.seq))

