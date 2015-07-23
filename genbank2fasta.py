from Bio import SeqIO
import argparse
from os import path



dTf options():
    parser = argparse.ArgumentParser(description="Extract sequences (genome,\
                                     nucleotides and aminoacids) from genbank\
                                     file")
    parser.add_argument("--inputgb", dest="inputgb", help="Input genebank file",
                        type=str)
    parser.add_argument("--outputfasta", dest="outputfasta",
                        help="Output fasta file", type=argparse.FileType("w"))
    parser.add_argument("--seq", dest="seq", help="Output sequence(s)",
                        choices = ["genome","gene","CDS","protein"],
                        default="genome",
                        type=str)

    return parser.parse_args( )

if __name__ == '__main__':
    args = options()
    try:
        if args.seq=="genome":
            SeqIO.convert(args.inputgb,"genbank",args.outputfasta,"fasta")
        else:
            record = SeqIO.read(args.inputgb,"genbank")
            for feature in record.features:
                if feature.type in ["gene","CDS"]:
                    seq = record.seq[feature.location.start:feature.location.end]
                    if feature.location.strand == -1:
                        seq = seq.reverse_complement( )
                    seq = str(seq).upper( )
                    if feature.type=="CDS" and args.seq=="protein":
                        args.outputfasta.write(">%s\t%s\t%s\n%s\n\n"%(feature.qualifiers["locus_tag"][0],
                                            feature.qualifiers["gene"][0] if "gene" in feature.qualifiers.keys( ) else "",
                                            feature.qualifiers["product"][0] if "product" in feature.qualifiers.keys( ) else "",
                                            str(feature.qualifiers["translation"][0])))
                    elif feature.type==args.seq:
                        args.outputfasta.write(">%s\t%s\t%s\n%s\n\n"%(feature.qualifiers["locus_tag"][0],
                                            feature.qualifiers["gene"][0] if "gene" in feature.qualifiers.keys( ) else "",
                                            feature.qualifiers["product"][0] if "product" in feature.qualifiers.keys( ) else "",
                                            seq))
    except IOError:
        print("%s doesn't exist or protected"%args.gbfile)


