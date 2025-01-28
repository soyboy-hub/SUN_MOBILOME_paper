# SUN_MOBILOME_paper
This repository contains python scripts, which was used in "Comparative Analysis of Active LTR Retrotransposons in Sunflower (Helianthus annuus L.): From Extrachromosomal Circular DNA Detection to Protein Structure Prediction" paper (https://doi.org/10.3390/ijms252413615). 
# Intact_TE_insertion_time.py
This script calculates insertion time of LTR-retrotransposons using Kimura 2 parameters method (doi: https://doi.org/10.1007/BF01731581). Requeres installed Biopython. Takes as input ".bed" file with LTRs in following format: {TEid}-LTR-{3/5}-{superfamily}-{clade} (for example TE_1-LTR-3-Gypsy-Tekay) and genome in ".fasta"/".fna"/".fa" format, gives as output TAB-delimited file with insertion time of each individual LTR-retrotransposon (in Mya).
# ParseCDHIT_align_get_consensus.py
This script extract sequence identifiers from .clstr file (obtained from CD-hit or MeShClust), conduct MSA (using mafft) and extract consensus. Requeres installed Biopython and mafft. Takes as input ".clstr" file and ".fasta" file, used for clusterization, gives as output ".fasta" file with consensus sequences of each cluster. 
