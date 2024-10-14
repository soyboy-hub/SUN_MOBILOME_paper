#!/usr/bin/env python
import os
import sys
import argparse

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

def get_args():
    ###get arguments from command line
    desc = (
        """Script for parsing TEtranscript output, the TE annotation file must be the same as used in TEtranscripts"""
    )
    epi = """Return file with filtered TEs by coverage approved by two repetitions. This script depends on samtools and bedtools."""
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("sigdiff", action="store", help='path to TEtranscripts_out_sigdiff_gene_TE.txt" file')
    parser.add_argument("TE_gtf", action="store", help='path to TE annotation in ".gtf" format')
    parser.add_argument("counttable", action="store", help='path to "TEtranscripts_out.cntTable" file')
    parser.add_argument("BAM1", action="store", help='path to the first SRR alignment in ".bam" format')
    parser.add_argument("BAM2", action="store", help='path to the second SRR alignment in ".bam" format')
    parser.add_argument("--minimum_cumulative_coverage", action="store", help='minimum covered by reads TE area in persent required to recognize TEs as expressed [by default is 50.0]')
    parser.add_argument("--minimum_reads", action="store", help='number of reads in each covered TE position [by default if 2]')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

#function to calculate coverage
def calc_cover(file, threshold, minimum_coverage):
    with open(file, 'r') as cov:
        dict = {}
        approved = []
        for line in cov:
            line = line.strip().split('\t')
            name = line[3]
            if name not in dict:
                dict[name] = []
            if name in dict:
                dict[name].append(line[5])
                
        for name in dict:
            satisfied_cov = 0
            lenght = len(dict[name])
            for position in dict[name]:
                if int(position) > int(minimum_coverage):
                    satisfied_cov += 1
            coverage_per_feature = (satisfied_cov/lenght)*100
            satisfied_cov = 0
            if coverage_per_feature > float(threshold):
                approved.append(name)
    return approved

#function to sort TEs by coverage based on two repetitions 
def sort_false_pos_by_coverage(bed, treat1, treat2, threshold, minimum_coverage):
    treat1_tmp = treat1.split('.bam')[0]
    treat2_tmp = treat2.split('.bam')[0]
    print('sorting...')
    os.system('samtools sort -@ 40 -o {0}_sort.bam {1}'.format(treat1_tmp, treat1))
    os.system('samtools sort -@ 40 -o {0}_sort.bam {1}'.format(treat2_tmp, treat2))
    print('idexing...')
    os.system('samtools index -@40 {}_sort.bam'.format(treat1_tmp))
    os.system('samtools index -@40 {}_sort.bam'.format(treat2_tmp))
    print('calculating coverage...')
    os.system('bedtools coverage -a {0} -b {1}_sort.bam -d > {0}_cov_treat1.txt'.format(bed, treat1_tmp))
    os.system('bedtools coverage -a {0} -b {1}_sort.bam -d > {0}_cov_treat2.txt'.format(bed, treat2_tmp))
    var1_approved = calc_cover(f'{bed}_cov_treat1.txt', threshold, minimum_coverage)
    var2_approved = calc_cover(f'{bed}_cov_treat2.txt', threshold, minimum_coverage)
    approved_by_two = [value for value in var1_approved if value in var2_approved]
    return approved_by_two

def main():
    print('Thanks for using xx')
    args = get_args()
    dict_TEs = {}
    TEs_coordinates = {}
    #making dictionary with names and families of LTR-RT
    with open(args.TE_gtf, 'r') as gtf:
        for line in gtf:
            line = line.split('\t')
            if line[2] == 'repeat_region':
                name = line[8].split('gene_id "')[1].split('";')[0]
                class_ = line[8].split('; family_id "')[1].split('";')[0]
                dict_TEs[name] = class_
                chr = line[0]
                start = line[3]
                end = line[4]
                TEs_coordinates[name] = f'{chr}:{start}-{end}'
    
    os.system('awk "NR>1" {0} > {0}_tmp'.format(args.sigdiff))
    os.system('cut -f1,3,6,7 {0}_tmp > {0}_clear'.format(args.sigdiff))
    
    file_new = f'{args.sigdiff}_clear'.split('/')[-1].split('.txt')[0]
    path_tmp = f'{args.sigdiff}_clear'.split('/')[0:-1]
    path = '/'.join(path_tmp)
    up_TEs = {}
    with open(f'{args.sigdiff}_clear', 'r') as tmp:
        for line in tmp:
            line = line.strip().split('\t')
            if 'TE_' in line[0] or 'repeat_region_' in line[0]:
                if float(line[1]) > 0:
                    name = line[0].split(':')[0]
                    superfamily = line[0].split(':')[1]
                    family = dict_TEs[name]
                    up_TEs[name] = f'{superfamily}_{family}'

    with open(f'{path}/tmp.bed', 'w') as tmp_bed:
        for name in TEs_coordinates:
            if name in up_TEs:
                chr = TEs_coordinates[name].split(':')[0]
                start = TEs_coordinates[name].split(':')[1].split('-')[0]
                end = TEs_coordinates[name].split(':')[1].split('-')[1]
                tmp_bed.write(f'{chr}\t{start}\t{end}\t{name}\n')

    ###sort false positive results by coverage
    print('false positive correction...')
    up_TEs_approvied = {}
    threshold = 50.0
    minimum_coverage = 2
    if args.minimum_cumulative_coverage:
        threshold = args.minimum_cumulative_coverage

    if args.minimum_reads:
        minimum_coverage = args.minimum_reads

    approvied = sort_false_pos_by_coverage(f'{path}/tmp.bed', args.BAM1, args.BAM2threshold, threshold, minimum_coverage)
    for name in approvied:
        up_TEs_approvied[name] = up_TEs[name]
   
    with open(f'{args.sigdiff}_clear', 'r') as tmp, \
    open(f'{path}/{file_new}_final.txt', 'w') as new:
        new.write('name\tsuperfamily\tfamily\tLogFC2\tpvalue\tpadj\n')
        for line in tmp:
            line = line.strip().split('\t')
            if 'TE_' in line[0] or 'repeat_region_' in line[0]:
                if float(line[1]) > 0:
                    name = line[0].split(':')[0]
                    if name in up_TEs_approvied:
                        superfamily = line[0].split(':')[1]
                        family = dict_TEs[name]
                        new_line = f'{name}\t{superfamily}\t{family}\t{line[1]}\t{line[2]}\t{line[3]}\n'
                        new.write(new_line)
    
    ###normalized average coverage per feature
    print('calculating normalized average coverage...')
    file_new1 = args.counttable.split('/')[-1].split('.cntTable')[0]

    os.system('awk "NR>1" {0} > {0}_tmp'.format(args.counttable))
    os.system('cut -f2,3,4,5 {0}_tmp > {0}_clear'.format(args.counttable))

    treat1_all = 0
    treat2_all = 0
    contr1_all = 0
    contr2_all = 0
    
    with open(f'{args.counttable}_clear', 'r') as tmp_all_counts:
        for line in tmp_all_counts:
            line = line.strip().split('\t')
            treat1_all += int(line[0])
            treat2_all += int(line[1])
            contr1_all += int(line[2])
            contr2_all += int(line[3])
    
    os.system('rm {0}_tmp'.format(args.sigdiff))
    os.system('rm {0}_clear'.format(args.sigdiff))
    os.system('rm {0}_clear'.format(args.counttable))
    os.system('rm {0}_tmp'.format(args.counttable))
    os.system('rm {0}/tmp.bed'.format(path))

if __name__ == "__main__":
    main()
