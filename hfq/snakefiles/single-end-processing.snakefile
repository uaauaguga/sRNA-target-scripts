dataset = "GSE77555-PRJNA310777"
adaptor = "TGGAATTCTCGGGTGCC"
window = 50
asm_id = "GCF_000007885.1"
sample_ids = open("data/GSE77555-PRJNA310777.txt").read().strip().split("\n")
strand = "forward"

#dataset = "GSE46118-PRJNA197291"
#adaptor = "TGGAATTCTCGGGTGCC"
#window = 50
#asm_id = "GCF_000005845.2"
#sample_ids = open("data/GSE46118-PRJNA197291-MG1655.txt").read().strip().split("\n")
#strand = "forward"

#dataset = "GSE46118-PRJNA197291"
#adaptor = "TGGAATTCTCGGGTGCC"
#window = 50
#asm_id = "GCF_000008865.2"
#sample_ids = open("data/GSE46118-PRJNA197291-Sakai.txt").read().strip().split("\n")
#strand = "forward"

rule all:
    input:
        #bam = expand("output/bam/{dataset}/{sample_id}.bam",dataset="GSE136110-PRJNA561330",sample_id=sample_ids)
        #txt = expand("output/strands/{dataset}/{sample_id}.txt",dataset=dataset,sample_id=sample_ids),
        peaks = expand("output/peaks.annotated/{dataset}/{sample_id}.{window}.bed",dataset=dataset,sample_id=sample_ids,window=window),
        #gb = expand("output/bedgraph/{dataset}/{sample_id}.{strand}.bedgraph",dataset=dataset,sample_id=sample_ids,strand = ["+","-"])

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
        """

rule cutadapt:
    input:
        fastq = "data/{dataset}/{sample_id}.fastq.gz",
    output:
        fastq = "output/trimmed/{dataset}/{sample_id}.fastq.gz", 
    params:
        adaptor = adaptor
    log:
        log = "output/trimmed/{dataset}/{sample_id}.log"
    threads: 4
    shell:
        """
        cutadapt -m 10 -j {threads} -a {params.adaptor} -o {output.fastq} {input.fastq} > {log.log} 2>&1
        """

rule mapping_single:
    input:
        fastq = "output/trimmed/{dataset}/{sample_id}.fastq.gz",
        bt2idx = f"genomes/bowtie2-index/{asm_id}.1.bt2"
    output:
        bam = "output/bam/{dataset}/{sample_id}.bam",
        fastq = "output/unmapped/{dataset}/{sample_id}.fastq.gz",
    log:
        log = "output/bam/{dataset}/{sample_id}.log"
    params:
        bt2idx = f"genomes/bowtie2-index/{asm_id}",
    threads: 4
    shell:
        """
        bowtie2 -p {threads} -U {input.fastq} -x {params.bt2idx} \
          --no-discordant --no-unal --un-gz output/unmapped/{wildcards.dataset}/{wildcards.sample_id}.fastq.gz 2> {log.log} | samtools view --output-fmt BAM -o {output.bam}
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
        bash scripts/check-bam-single.sh {input.bam} {input.bed} > {output.strand}
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
        scripts/bam2bed.py --layout single -i {input.bam} -o {output.bed} -s {params.strand} > {log.log} 2>&1
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
