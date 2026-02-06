rule all:
    input:
        "SRR2584857-assembly.fa",
        "SRR2584857_quast",
        "SRR2584857_reads.x.SRR2584857_assembly.bam.sorted"

rule sketch:
    input:
        "SRR2584857-reads.sig.zip",
        "SRR2584857-assembly.sig.zip"

rule assemble:
    input:
        r1 = "SRR2584857_1.fastq.gz",
        r2 = "SRR2584857_2.fastq.gz"
    output:
        directory("SRR2584857_assembly")
    threads: 8
    shell: """
        megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t {threads} -o {output}
    """

rule rename_contigs:
    input: "SRR2584857_assembly"
    output: "SRR2584857-assembly.fa"
    shell: """
        cp {input}/final.contigs.fa {output}
    """

rule quast_assembly:
    input: "SRR2584857-assembly.fa"
    output: directory("SRR2584857_quast")
    threads: 8
    shell: """
        quast {input} -o {output} --threads {threads}
    """

rule prokka_assembly:
    input: "SRR2584857-assembly.fa"
    output: directory("SRR2584857_annot")
    threads: 8
    shell: """
        prokka --prefix {output} {input} --cpu {threads}
    """

rule map_to_assembly:
    input:
        ref = "SRR2584857-assembly.fa",
        r1 = "SRR2584857_1.fastq.gz",
        r2 = "SRR2584857_2.fastq.gz",
    output:
        "SRR2584857_reads.x.SRR2584857_assembly.bam",
    threads: 8
    shell: """
        minimap2 -ax sr -t {threads} {input.ref} {input.r1} {input.r2} | \
           samtools view -b - -o {output}
    """

rule sort_bam:
    input:
        "SRR2584857_reads.x.SRR2584857_assembly.bam",
    output:
        "SRR2584857_reads.x.SRR2584857_assembly.bam.sorted",
    shell: """
        samtools sort {input} -o {output}
    """

rule sketch_assembly:
    input: "SRR2584857-assembly.fa"
    output: "SRR2584857-assembly.sig.zip"
    shell: """
        sourmash sketch dna {input} -o {output} --name assembly
    """

rule sketch_reads:
    input:
        r1 = "SRR2584857_1.fastq.gz",
        r2 = "SRR2584857_2.fastq.gz",
    output: "SRR2584857-reads.sig.zip"
    shell: """
        sourmash sketch dna -p abund {input} -o {output} --name reads
    """
