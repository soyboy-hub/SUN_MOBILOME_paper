#!/usr/bin/env python

import os
import sys
import argparse
import re
from itertools import groupby
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

def get_args():
    ###get arguments from command line
    desc = (
        """Extract sequence identifiers from .clstr file, conduct MSA and extract consensus"""
    )
    epi = """Return consensus sequence for each cluster.
             A typical .clstr file look like:
             >Cluster 0
             0  100nt, >seq1...*
             >Cluster 1
             0  102nt, >seq2...*
             1  101nt, >seq3...
             This script depends on mafft and Biopython package.
          """
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("clusterfile", action="store", help='file in ".clstr" format')
    parser.add_argument("fastafile", action="store", help='file in ".fasta"/".fna"/".fa" format, containing clustered sequences')
    parser.add_argument("threads", action="store", help='number of threads for MSA')
    parser.add_argument("outfile", action="store", help='output name and path for consensus file')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def main():
    print('Thanks for using xx')
    args = get_args()
    ###parcing .clsrt file, producing MSA and extracting consensus
    os.system('mkdir tmp_ParseCDHIT_align_get_consensus')
    os.system('touch {0}'.format(args.outfile))
    index_fasta = SeqIO.index(args.fastafile, 'fasta')
    ###parse clusters
    with open(args.clusterfile, "r") as cfh:
        counter = 0
        groups = {}
        groups_nested = [
            list(group)
            for key, group in groupby(cfh, lambda line: line.startswith(">Cluster"))
            if not key
        ]

        for group in groups_nested:
            if len(group) < 2:
                with open(args.outfile, 'a') as consens:
                    name = group[0].split('\t')[1].split('>')[1].split('...')[0]
                    sequence = str(index_fasta[name].seq)
                    consens.write(f'>{name}_single\n{sequence}\n')
                continue
            counter += 1
            groups[f'cluster{counter}'] = []
            for id in group:
                name = id.split('\t')[1].split('>')[1].split('...')[0]
                groups[f'cluster{counter}'].append(name)
        for group in groups:
            print(f'processing {group}')
            with open(f'tmp_ParseCDHIT_align_get_consensus/{group}.fasta', 'w') as fasta_cluster:
                for seq_name in groups[group]:
                    sequence = str(index_fasta[seq_name].seq)
                    fasta_cluster.write(f'>{seq_name}\n{sequence}\n')
            os.system('mafft --clustalout --quiet --thread {0} tmp_ParseCDHIT_align_get_consensus/{1}.fasta > tmp_ParseCDHIT_align_get_consensus/{1}.aln'.format(args.threads, group))
            alignment = AlignIO.read(open(f'tmp_ParseCDHIT_align_get_consensus/{group}.aln'), 'clustal')
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = str(summary_align.dumb_consensus(threshold=0.5, ambiguous='N').upper())
            with open(args.outfile, 'a') as consens:
                consens.write(f'>{group}_multiple\n{consensus}\n')
            os.system('rm tmp_ParseCDHIT_align_get_consensus/{0}.fasta'.format(group))
            os.system('rm tmp_ParseCDHIT_align_get_consensus/{0}.aln'.format(group))
    os.system('rm tmp_ParseCDHIT_align_get_consensus')

if __name__ == "__main__":
    main()
