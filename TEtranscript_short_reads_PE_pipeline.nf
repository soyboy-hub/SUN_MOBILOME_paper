params.fastqc_row_output = "$projectDir/fastqc_row_results"
params.fastqc_trim_output="$projectDir/fastqc_trim_results"
params.gtf_genes = ""
params.gtf_TE = ""
params.genome = ""
params.index_output = "$projectDir/STAR_index"
params.bam_output = "$projectDir"
params.treatments = ""
params.control = ""
treat_bam_label = params.treatments.split(",")
con_bam_label = params.control.split(",")
control_map= Channel.of(con_bam_label).collate(1)
treat_map=Channel.of(treat_bam_label).collate(1)
params.TEtranscript_output= ""
stress_outdir = "$projectDir/${params.TEtranscript_output}"
params.api = "395aee5592e5c445fa836e7c037833720208"
params.accession = ""
accession=params.accession.split(",").toList()



if (params.help) {
    println """
                  TE expression short reads pipeline
    =================================================================
    Usage:
    nextflow TEtranscript_short_reads_PE_pipeline.nf --accession [SRA accessions] --genome [path to genome] 
                             --gtf_genes [path to genes GTF] --gtf_TE [path to TE GTF] --stress_outdir [output directory for TEtranscript] 
                             --treatments [sample IDs for treatments] --control [sample IDs for controls]
    
    Параметры:
    --accession             list of SRR accessions, must be comma-separated without spaces (SRR1,SRR2,SRR3,SRR4)
    --genome                path to genome sequences file (.fastq/.fa/.fna)
    --gtf_genes             path to genes GTF (.gtf)
    --gtf_TE                path to TE GTF (.gtf)
    --treatments            list of sample IDs for treatments, must be comma-separated without spaces (SRR1,SRR2)
    --control               list of sample IDs for controls, must be comma-separated without a spaces (SRR3,SRR4)
    --stress_outdir         output directory for TEtranscript
    """
    System.exit(0)
}

log.info """\
    SHORT READS PIPELINE
    =================
    reads: ${params.fastq}
    SRA accession: ${accession}
    genome: ${params.genome}
    GTF_gene: ${params.gtf_genes}
    GTF_TE: ${params.gtf_TE}
    treatments: ${treat_bam_label}
    control: ${con_bam_label}
    TEtranscript_output: ${stress_outdir}
    """
    .stripIndent()

process fastqc_row {
  publishDir params.fastqc_row_output, mode: 'move'
  tag "fastqc_row on $sample_id"
  input:
  tuple val(sample_id), path(reads_file)

  output:
  tuple val(sample_id), path ("fastqc_${sample_id}_logs")

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file[0]} -t 30 && fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file[1]} -t 30
  """
}
process trimmomatic {
  tag "trimmomatic on $sample_id"
  input:
  tuple val(sample_id), path(reads_file)

  output:
  tuple val(sample_id), path("${sample_id}_R1_paired.fastq.gz"), path("${sample_id}_R2_paired.fastq.gz"), emit: read_pairs_trimmed

  script:
  """
  TrimmomaticPE \
  -threads 85 \
  $reads_file \
  ${sample_id}_R1_paired.fastq.gz \
  ${sample_id}_R1_unpaired.fastq.gz \
  ${sample_id}_R2_paired.fastq.gz \
  ${sample_id}_R2_unpaired.fastq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10
  """
}
process fastqc_trim {
  publishDir params.fastqc_trim_output, mode: 'move'
  tag "fastqc_trim on $sample_id"
  input:
  tuple val(sample_id), path("${sample_id}_R1_paired.fastq.gz"), path("${sample_id}_R2_paired.fastq.gz")

  output:
  tuple val(sample_id), path ("fastqc_${sample_id}_logs")

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${sample_id}_R1_paired.fastq.gz -t 30 && fastqc -o fastqc_${sample_id}_logs -f fastq -q ${sample_id}_R2_paired.fastq.gz -t 30
  """
}
process STAR_genomeGenerate {
  publishDir params.index_output, mode: 'move'
  tag "index generate on $reference_genome and $reference_gtf"
  input:
  path(reference_genome)
  path(reference_gtf)

  output:
  path ("STAR_index"), emit: index_dir

  script:
  """
  mkdir STAR_index
  STAR --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles ${reference_genome} \
    --sjdbGTFfile ${reference_gtf} \
    --sjdbOverhang 50 \
    --runThreadN 80 \
    --limitGenomeGenerateRAM 171798691840
  """
}
process STAR_aligment {
  memory="420G"
  // publishDir params.bam_output, mode: 'move'
  tag "STAR aligment on $sample_id"
  input:
  each(reference_gtf_genes)
  path(index_dir)
  tuple val(sample_id), path(trim_file1), path(trim_file2)

  output:
  tuple val(sample_id), path("${sample_id}Aligned.out.bam")

  script:
  """
  STAR --genomeDir ${params.index_output}/STAR_index \
         --readFilesIn $trim_file1 $trim_file2 \
         --runThreadN 80 \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id} \
         --sjdbGTFfile ${reference_gtf_genes} \
         --outSAMtype BAM Unsorted \
         --limitBAMsortRAM 80000000000000
  """
}
process TE_transcript {
  publishDir stress_outdir, mode: 'move'
  tag "TEtranscript to treatments ($treat) and control ($con)"
  input:
  each(reference_gtf_genes)
  each(reference_gtf_TE)
  path(treat)
  path(con)
  
  output:
  path ("TEtranscripts_output")

  script:
  """
  mkdir TEtranscripts_output
  TEtranscripts -t $treat -c $con --GTF $reference_gtf_genes --TE $reference_gtf_TE --outdir TEtranscripts_output
  """
}


workflow {
  //reads=Channel.fromFilePairs(params.fastq)
  reads = Channel.fromSRA(accession, apiKey: params.api, protocol: 'ftp').view()
  fastqc_row(reads)

  trim_reads=trimmomatic(reads)
  fastqc_trim(trim_reads)

  ref_ch = Channel.fromPath(params.genome).view()
  gtf_ch = Channel.fromPath(params.gtf_genes).view()
  TE_ch = Channel.fromPath(params.gtf_TE).view()
  STAR_genomeGenerate(ref_ch, gtf_ch)
  index_ch=Channel.fromPath(params.index_output).collect().view()
  aligment=STAR_aligment(gtf_ch, index_ch, trim_reads.read_pairs_trimmed).collect().flatten().buffer(size:2).view()
  
  TE_control = control_map.combine(aligment, by: 0).collect{it[1]}.view()
  TE_treat = treat_map.combine(aligment, by: 0).collect{it[1]}.view()
  TE_transcript(gtf_ch, TE_ch, TE_treat, TE_control)
}
