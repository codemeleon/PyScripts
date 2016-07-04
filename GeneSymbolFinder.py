# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 14:02:06 2016

@author: devil
"""
### The program to find the possible gene symbol based on sequences
#import recipy
from Bio import SeqIO,Entrez
from os import system,remove
import click
import time
import sys

@click.command()
@click.option("--infile",type=str,help="Input fasta file",default="./unknown.fa")
@click.option("--inseq",type=str,help="Input sequence type",default="nuc")
@click.option("--outfile",type=str,help="Result Output File",default="./unknownx")
@click.option("--email",type=str,help="email for ncbi online query")
@click.option("--evalue",type=float,help="Select p-values",default = 1e-5)
def run(infile,outfile,evalue,email,inseq):
    """Gene symbool finder for given sequences.
    Requires blast suit installed and active internet connection"""
    if len(email) < 5 or '@' not in email:
        click.echo("Your email doesn't look like email")
    Entrez.email = email
    outfile = open(outfile,"w")
    for rec in SeqIO.parse(infile,'fasta'):
        with open("tmp.fa","w") as tmp_fa:
            tmp_fa.write(">%s\n%s\n"%(rec.id,rec.seq))
        if inseq == 'nuc':
            system("blastx -db nr -query tmp.fa -evalue 1E-10 -remote -out tmp.xml -outfmt 5 -max_target_seqs 20")
        elif inseq == 'aa':
            system("blastx -db nr -query tmp.fa -evalue 1E-10 -remote -out tmp.xml -outfmt 5 -max_target_seqs 20")
        else:
            click.echo("Unknown sequence format")
            sys.exit(1)
        time.sleep(5)
        gene_symbols = []
        with open("tmp.xml") as tmpxml:
            for line in tmpxml:
                if 'Hit_id' not in line:continue
                gi = line.split('|')[1]
                time.sleep(3)
                handle = Entrez.efetch(db='protein',id=gi,rettype='gb',retmode='text')
                record = SeqIO.read(handle,"genbank")
                handle.close( )
                for feature in record.features:
                    if feature.type == 'Region':
                        gene_symbols.append(feature.qualifiers['region_name'][0])
                        break
        gene_symbol_set = set(gene_symbols)
        gene_name = ''
        gene_name_count = 0
        for gene in gene_symbol_set:
            count = gene_symbols.count(gene)
            if gene_name_count < count:
                gene_name_count = count
                gene_name = gene
        outfile.write("%s\t%s\n"%(rec.id,gene_name))
    remove("tmp.fa")
    remove("tmp.xml")
if __name__=='__main__':
    run()
    
            