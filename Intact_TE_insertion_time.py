#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
import argparse
from math import log, sqrt

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

def get_args():
    ###get arguments from command line
    desc = (
        """Extract 3'LTR/5'LTR sequences from BED file, align and calculate genomic divergence (in Mya) using Kimura 2 parameters method. BED file should contain LTR coordinates, header should be "{TEid}-LTR-{3/5}-{superfamily}-{clade}" (i.e. TE_1-LTR-3-Gypsy-Tekay, TE_1-LTR-5-Gypsy-Tekay)"""
    )
    epi = """Return table with insertion time of each intact LTR-RT from file with annotations.
    
    This script depends on mafft and bedtools, Biopython package.
          """


    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("bed", action="store", help='file in BED format')
    parser.add_argument("fastafile", action="store", help='genome in ".fasta"/".fna"/".fa" format')
    parser.add_argument("outfile", action="store", help='output name and path for table')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

#function to estimate TE insertion time by Kimura 2 parameters
def estimate_insert_time(LTR3_seq, LTR5_seq):
    r = 1E-08
    with open('TE_tmp.fasta', 'w') as tmp:
        tmp.write(f'>LTR3\n{LTR3_seq}\n')
        tmp.write(f'>LTR5\n{LTR5_seq}\n')
    os.system('mafft --clustalout --quiet --thread 40 TE_tmp.fasta > TE_tmp.aln')
    with open("TE_tmp.aln", "r") as aln:
        pairs = []
        for record in SeqIO.parse(aln, "clustal"):
            if record.id == 'LTR3':
                seq1 = record.seq
            if record.id == 'LTR5':
                seq2 = record.seq
        for x in zip(seq1, seq2):
            if '-' not in x:
                pairs.append(x)
        ts_count = 0
        tv_count = 0
        length = len(pairs)
        transitions = ["ag", "ga", "ct", "tc"]
        transversions = ["ac", "ca", "at", "ta",
                         "gc", "cg", "gt", "tg"]
        for (x, y) in pairs:
            if x + y in transitions:
                ts_count += 1
            elif x + y in transversions:
                tv_count += 1
        p = float(ts_count) / length
        q = float(tv_count) / length
        try:
            d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
        except ValueError:
            print("Tried to take log of a negative number")
            return None
        insertion_time = d/(2*r*1000000)
        os.system('rm TE_tmp.aln')
        os.system('rm TE_tmp.fasta')
        return insertion_time

def main():
    print('Thanks for using xx')
    args = get_args()
    os.system('bedtools getfasta -fi {0} -bed {1} -name -fo {2}'.format(args.fastafile, args.bed, 'tmp.fasta'))
    with open(args.bed, 'r') as bed, \
    open(args.outfile, 'w') as new:
        used = []
        indx_fasta = SeqIO.index('tmp.fasta', 'fasta')
        for line in bed:
            name = line.strip().split('\t')[3]
            teid = name.split('-')[0]
            if teid not in used:
                used.append(teid)
                superf = name.split('-')[-2]
                clad = name.split('-')[-1]
                try:
                    LTR3_seq = str(indx_fasta[f'{teid}-LTR-3-{superf}-{clad}'].seq)
                    LTR5_seq = str(indx_fasta[f'{teid}-LTR-5-{superf}-{clad}'].seq)
                except:
                    continue
                insert_time = estimate_insert_time(LTR3_seq, LTR5_seq)
                if insert_time == '-0.0':
                    insert_time = '0.0'
                print(f'processing {teid}')
                new.write(f'{teid}\t{superf}\t{clad}\t{insert_time}\n')
            else:
                continue
    os.system('rm tmp.fasta')

if __name__ == "__main__":
    main()
