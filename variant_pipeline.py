#use global wildcard make sample list#
SAMPLE_LIST, = glob_wildcards("0.rawdata/{sample}_L001_R1_001.fastq.gz")
ref ="data/indexRef.fa"

rule all:
    input:
        expand("4.snp/{sample}.snp.tsv",sample=SAMPLE_LIST),
        expand("5.mpileup/{sample}.mpileup.tsv",sample=SAMPLE_LIST),
        expand("6.consensus/{sample}.fa",sample=SAMPLE_LIST)

rule qc:
    input:
        seq1="0.rawdata/{sample}_L001_R1_001.fastq.gz",
        seq2="0.rawdata/{sample}_L001_R2_001.fastq.gz"
    output:
        qc1="1.qc/{sample}.qc.1.fq.gz",
        qc2="1.qc/{sample}.qc.2.fq.gz",
        report="1.qc/{sample}.html"

    shell:
        """
        fastp -i {input.seq1} -I {input.seq2} -o {output.qc1} -O {output.qc2} -q 20 -h {output.report}
        """
rule cutadapt:
    input:
        seq1=rules.qc.output.qc1,
        seq2=rules.qc.output.qc2
    output:
        cutseq1="2.cutadapt/{sample}.cutadapt.1.fq",
        cutseq2="2.cutadapt/{sample}.cutadapt.2.fq"
    shell:
        """
        bwa mem -t 20 -a {ref} {input.seq1} {input.seq2} \
        | python ./tools/trim_primer_parts.py data/primer.bed \
        {output.cutseq1} {output.cutseq2}
        """
rule mapping:
    input:
        seq1=rules.cutadapt.output.cutseq1,
        seq2=rules.cutadapt.output.cutseq2
    output:
        bamfile="3.mapping/{sample}.bam"

    shell:
        """
        bwa mem -t 20 -a {ref} {input.seq1} {input.seq2} \
        | samtools view -@ 4 -bS -F 4 \
        | samtools sort -@ 4 -o {output.bamfile};
        """
        """
        samtools index {output.bamfile}
        """
rule mpileup:
    input:
        bamfile=rules.mapping.output.bamfile
    output:
        snptab="4.snp/{sample}.snp.tsv",
        mpileuptab="5.mpileup/{sample}.mpileup.tsv",
        consensusseq="6.consensus/{sample}.fa"

    shell:
        """
        samtools mpileup -Q 20 -d 100000 --reference {ref} {input.bamfile} \
        | ivar variants -r {ref} -p {output.snptab} -t 0;
        """
        """
        bcftools mpileup --threads 2 -Q 20 -a FORMAT/DP,FORMAT/AD,INFO/AD \
        -I -d 100000 -f {ref} {input.bamfile} > {output.mpileuptab};
        """
        """
        samtools mpileup {input.bamfile} \
        | ivar consensus -p {output.consensusseq} -q 20 -t 0 -n N -m 10
        """
