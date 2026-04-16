trim5p = None
#dataset = "GSE136110-PRJNA561330"
#sample_ids = open("data/GSE136110-PRJNA561330.txt").read().strip().split("\n")
#asm_id = "GCF_000006765.1"
#strand = "forward"
#adaptor_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
#adaptor_2 = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

#dataset = "GSE163336-SRP298148"
#sample_ids = open("data/GSE163336-SRP298148.txt").read().strip().split("\n")
#asm_id = "GCF_000210855.2"
#strand = "forward"
#adaptor_1 = ""
#adaptor_2 = ""

#dataset = "SRP154464-PRJNA481981"
#sample_ids = open("data/SRP154464-PRJNA481981.txt").read().strip().split("\n")
#asm_id = "GCF_000022005.1" 
#adaptor_1 = ""
#adaptor_2 = ""
#strand = "forward"

#dataset = "GSE74425-PRJNA300395"
#asm_id = "GCF_000210855.2" 
#sample_ids = open("data/GSE74425-PRJNA300395.txt").read().strip().split("\n")
#adaptor_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#adaptor_2 = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"
#strand = "forward"

'''
dataset = "GSE129868-PRJNA533193"
asm_id = "GCF_000026965.1"
sample_ids = open("data/GSE129868-PRJNA533193.txt").read().strip().split("\n")
adaptor_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adaptor_2 = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"
strand = "forward"
'''


dataset = "GSE131520-SRP198991"
sample_ids = open("data/GSE131520-SRP198991.txt").read().strip().split("\n")
adaptor_1 = ""
adaptor_2 = ""
strand = "reverse"
asm_id = "GCF_000005845.2"


#dataset = "GSE216133-SRP403529"
#sample_ids = open("data/GSE216133-SRP403529.txt").read().strip().split("\n")
#adaptor_1 = ""
#adaptor_2 = ""
#strand = "reverse"
#asm_id = "GCF_000006765.1"

#../../RNA-RNA-interaction/bacteria-analysis/snakefiles/RIL-seq-processing.snakefile

#dataset = "GSE198671-SRP364136"
#sample_ids = open("data/GSE198671-SRP364136.txt").read().strip().split("\n")
#adaptor_1 = ""
#adaptor_2 = ""
#strand = "reverse"
#asm_id = "GCF_000006745.1"


#dataset = "GSE213005-SRP396419"
#sample_ids = open("data/GSE213005-SRP396419.txt").read().strip().split("\n")
#asm_id = "GCF_000932055.2"
#adaptor_1 = ""
#adaptor_2 = ""
#strand = "forward"

#dataset = "GSE244638-PRJNA1023942" 
#sample_ids = open("data/GSE244638-PRJNA1023942.txt").read().strip().split("\n")
#asm_id = "GCF_000016305.1"
#adaptor_1 = "AGATCGGAAGAGCACA"

#adaptor_2 = "ACCCGTCTTAGATCGG"
#sample_ids = ["SRR26281961"]
#adaptor_2  = "ACAACTCGCAGATCGG"
#sample_ids = ["SRR26281962"]
#adaptor_2  = "ACCAAGTCGAGATCGG"
#sample_ids = ["SRR26281963"]
#adaptor_2  = "AATCACTTGAGATCGG"
#sample_ids = ["SRR26281964"]
#trim5p = 9
#strand = "reverse"

'''
asm_id = "GCF_000005845.2"
dataset = "PRJEB10931-ERP012236"
adaptor_1 = ""
adaptor_2 = ""
sample_ids = open("data/PRJEB10931-ERP012236.txt").read().strip().split("\n")
strand = "reverse"

asm_id = "GCF_000742755.1"
dataset = "GSE243246-PRJNA1017475"
adaptor_1 = ""
adaptor_2 = ""
sample_ids = open("data/GSE243246-PRJNA1017475.txt").read().strip().split("\n")
strand = "reverse"
'''

#asm_id = "GCF_000210855.2"
#dataset = "GSE234792-PRJNA983241"
#adaptor_1 = ""
#adaptor_2 = ""
#sample_ids = open("data/GSE234792-PRJNA983241.txt").read().strip().split("\n")
#strand = "reverse"

window = 50
rule all:
    input:
        txt = expand("output/strands/{dataset}/{sample_id}.txt",dataset=dataset,sample_id=sample_ids),
        peaks = expand("output/peaks.annotated/{dataset}/{sample_id}.{window}.bed",dataset=dataset,sample_id=sample_ids,window=window),
        bg = expand("output/bedgraph/{dataset}/{sample_id}.{strand}.bedgraph",dataset=dataset,sample_id=sample_ids,strand = ["+","-"])

rule build_index:
    input:
        fasta = f"genomes/fasta/{asm_id}.fa"
    output:
        bt2idx = f"genomes/bowtie2-index/{asm_id}.1.bt2"
    params:
        asm_id = asm_id 
    log:
        log = f"genomes/bowtie2-index/{asm_id}"
    shell:
        """
        bowtie2-build {input.fasta} genomes/bowtie2-index/{params.asm_id} > {log.log} 2>&1
        """

rule gff2bed:
    input:
        gff = f"genomes/gff/{asm_id}.gff" 
    output:
        bed = f"genomes/bed/{asm_id}.bed"
    log:
        log = f"genomes/bed/{asm_id}.log"
    shell:
        """
        scripts/gff2bed.py --gff {input.gff}  --bed {output.bed} --feature CDS --name Name,locus_tag,gene > {log.log} 2>&1
        sort -k1,1 -k2,2n -o {output.bed}  {output.bed}
        """

rule cutadapt:
    input:
        fastq_1 = "data/{dataset}/{sample_id}_1.fastq.gz",
        fastq_2 = "data/{dataset}/{sample_id}_2.fastq.gz"
    output:
        fastq_1 = "output/trimmed/{dataset}/{sample_id}_1.fastq.gz", 
        fastq_2 = "output/trimmed/{dataset}/{sample_id}_2.fastq.gz"
    params:
        adaptor_1 = adaptor_1,
        adaptor_2 = adaptor_2,
        trim = f"--cut {trim5p}" if trim5p is not None else "--cut 0"
    log:
        log = "output/trimmed/{dataset}/{sample_id}.log"
    threads: 4
    shell:
        """
        cutadapt --trimmed-only {params.trim} -m 10 -j {threads} -a {params.adaptor_1} -A {params.adaptor_2} -o {output.fastq_1} -p {output.fastq_2} {input.fastq_1} {input.fastq_2} > {log.log} 2>&1
        """

rule mapping_paired:
    input:
        fastq_1 = "output/trimmed/{dataset}/{sample_id}_1.fastq.gz",
        fastq_2 = "output/trimmed/{dataset}/{sample_id}_2.fastq.gz",
        bt2idx = f"genomes/bowtie2-index/{asm_id}.1.bt2"
    output:
        bam = "output/bam/{dataset}/{sample_id}.bam",
        fastq_1 = "output/unmapped/{dataset}/{sample_id}_1.fastq.gz",
        fastq_2 = "output/unmapped/{dataset}/{sample_id}_2.fastq.gz"
    log:
        log = "output/bam/{dataset}/{sample_id}.log"
    params:
        bt2idx = f"genomes/bowtie2-index/{asm_id}",
    threads: 6
    shell:
        """
        bowtie2 -p {threads} -1 {input.fastq_1} -2 {input.fastq_2} -x {params.bt2idx} \
          --no-discordant --no-unal --un-conc-gz output/unmapped/{wildcards.dataset}/{wildcards.sample_id}_%.fastq.gz 2> {log.log} | samtools view --output-fmt BAM -o {output.bam}
        """


rule annotate_reads:
    input:
        bed = f"genomes/bed/{asm_id}.bed",
        bam = "output/bam/{dataset}/{sample_id}.bam"
    output:
        strand = "output/strands/{dataset}/{sample_id}.txt"
    log:
        log = "output/strands/{dataset}/{sample_id}.log"
    shell:
        """
        bash scripts/check-bam.sh {input.bam} {input.bed} > {output.strand}
        """

rule bam2bed:
    input:
        bam = "output/bam/{dataset}/{sample_id}.bam"
    output:
        bed = "output/bed/{dataset}/{sample_id}.bed"
    log:
        log = "output/bed/{dataset}/{sample_id}.log"
    params:
        strand = strand
    shell:
        """
        scripts/bam2bed.py -i {input.bam} -o {output.bed} -s {params.strand} > {log.log} 2>&1
        #sort -k1,1 -k2,2n -o {output.bed} {output.bed}
        """


rule sort:
    input:
        bed = "output/bed/{dataset}/{sample_id}.bed"
    output:
        bed = "output/bed.sorted/{dataset}/{sample_id}.bed"
    shell:
        """
        sort -k1,1 -k2,2n {input.bed} > {output.bed}
        """

rule index_genome:
    input:
        fa = f"genomes/fasta/{asm_id}.fa"
    output:
        fai = f"genomes/fasta/{asm_id}.fa.fai"
    shell:
        """
        samtools faidx {input.fa}
        """

rule get_coverage:
    input:
        bed = "output/bed.sorted/{dataset}/{sample_id}.bed",
        fai = f"genomes/fasta/{asm_id}.fa.fai"
    output:
        plus = "output/bedgraph/{dataset}/{sample_id}.+.bedgraph",
        minus = "output/bedgraph/{dataset}/{sample_id}.-.bedgraph"
    shell:
        """
        bedtools genomecov -i {input.bed} -g {input.fai} -strand + -bg > {output.plus}
        bedtools genomecov -i {input.bed} -g {input.fai} -strand - -bg > {output.minus}
        """

rule callpeak:
    input:
        bed = "output/bed.sorted/{dataset}/{sample_id}.bed"
    output:
        peak = "output/peaks/{dataset}/{sample_id}.{window}.bed"
    params:
        window = window        
    shell:
        """
        /BioII/lulab_b/jinyunfan/anaconda3/envs/peak-caller/bin/Piranha -cluster_dist 50 -s -p_threshold 0.01 -no_pval_correct \
        -background_thresh 0.9 -dist ZeroTruncatedNegativeBinomial -b {params.window} {input.bed} > {output.peak}
        """

rule annotate_peak:
    input:
        bed = "output/peaks/{dataset}/{sample_id}.{window}.bed",
        CDS = f"genomes/bed/{asm_id}.bed",
        fai = f"genomes/fasta/{asm_id}.fa.fai"
    output:
        bed = "output/peaks.annotated/{dataset}/{sample_id}.{window}.bed"
    log:
        log = "output/peaks.annotated/{dataset}/{sample_id}.{window}.log"
    shell:
        """
        sort -k1,1 -k2,2n -o {input.bed} {input.bed}
        scripts/annotate-intervals.py -g {input.CDS} -b {input.bed} -o {output.bed} -c {input.fai} > {log.log} 2>&1
        """
