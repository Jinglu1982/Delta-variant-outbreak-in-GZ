#use global wildcard make sample list#
SAMPLE_LIST, = glob_wildcards("0.rawdata/{sample}_L001_R1_001.fastq.gz")
ref =""data/indexRef.fa"
primer="data/v3primer_indexref.bed"

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
        bamfile="2.cutadapt/{sample}.sorted.bam",
        trimbam="2.cutadapt/{sample}.trim.bam"
    params:
        prefix="2.cutadapt/{sample}.trimmed"
    shell:
        """
		bwa mem -t 20 -a {ref2} {input.seq1} {input.seq2} | samtools view -@ 4 -bS -F 4 | samtools sort -@ 4 -o {output.bamfile};
		ivar trim -b {primer} -p {params.prefix} -i {output.bamfile} -q 15 -s 4 -e -m 10;
		samtools sort {params.prefix}.bam -o {output.trimbam};
		rm {params.prefix}.bam {output.bamfile};
        """
rule mpileup:
    input:
        bamfile=rules.cutadapt.output.trimbam
    output:
        snptab="4.snp/{sample}.snp.tsv",
        mpileuptab="5.mpileup/{sample}.mpileup.tsv",
        consensusseq="6.consensus/{sample}.fa"

    shell:
        "samtools mpileup -Q 20 -d 100000 --reference {ref2} {input.bamfile} | ivar variants -r {ref2} -p {output.snptab} -t 0;"
        "bcftools mpileup --threads 2 -Q 20 -a FORMAT/DP,FORMAT/AD,INFO/AD -I -d 100000 -f {ref2} {input.bamfile} > {output.mpileuptab};"
        "samtools mpileup {input.bamfile} | ivar consensus -p {output.consensusseq} -q 20 -t 0 -n N -m 10"
