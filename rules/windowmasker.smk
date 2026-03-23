def get_fasta(wildcards):
    return manifest_df.loc[(wildcards.sample, wildcards.hap), "FASTA"]


rule index_fasta:
    input:
        fasta = get_fasta
    output:
        fasta = "results/{sample}/work/windowmasker/raw_fasta/{hap}.fasta",
        fai = "results/{sample}/work/windowmasker/raw_fasta/{hap}.fasta.fai"
    threads: 1
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        tmp_fasta="{output.fasta}.tmp"
        # to prevent errors in sedef.sh due to "#" in the contig names.
        sed '/^>/ s/#/__/g' "$(readlink -f {input.fasta})" > "$tmp_fasta"
        mv $tmp_fasta {output.fasta}
        samtools faidx {output.fasta}
        """

rule dust_count:
    input:
        ref = rules.index_fasta.output.fasta,
    output:
        counts = "results/{sample}/work/windowmasker/outputs/{hap}.counts",
    resources:
        mem = 16,
        hrs = 120,
    threads: 1
    singularity:
        "docker://eichlerlab/ncbi-windowmasker:1.0"
    params:
        mem = 16384,
        smem = 2048
    shell: """
        tmp_output="{output.counts}.tmp"
        windowmasker \
        -mem {params.mem} \
        -smem {params.smem} \
        -mk_counts \
        -infmt fasta \
        -sformat obinary \
        -in {input.ref} \
        -out $tmp_output
        mv $tmp_output {output.counts}
        """

rule run_windowmasker:
    input:
        counts = rules.dust_count.output.counts,
        ref = rules.index_fasta.output.fasta,
        fai = rules.index_fasta.output.fai,
    output:
        intervals = "results/{sample}/work/windowmasker/outputs/{hap}.intervals",
    resources:
        mem = 16,
        hrs = 12,
    threads: 1
    singularity:
        "docker://eichlerlab/ncbi-windowmasker:1.0"
    shell: """
        windowmasker \
        -infmt fasta \
        -ustat {input.counts} \
        -dust T \
        -outfmt interval \
        -in {input.ref} \
        -out {output.intervals}
        """

rule run_windowmasker_bed:
    input:
        intervals = rules.run_windowmasker.output.intervals,
    output:
        bed = "results/{sample}/work/windowmasker/outputs/{hap}.dust.bed",
        sorted_bed = "results/{sample}/work/windowmasker/outputs/{hap}.mask_repeat.sorted.bed",
    resources:
        mem=16,
        hrs=12,
    params:
        script_dir=f"{SNAKEMAKE_DIR}/scripts"
    threads: 1
    singularity:
        "docker://eichlerlab/ncbi-windowmasker:1.0"
    shell: """
        tmp_out_bed="{output.bed}.tmp"
        tmp_out_sort_bed="{output.sorted_bed}.tmp"
        perl {params.script_dir}/windowmasker_to_bed.pl < {input.intervals} > $tmp_out_bed
        mv $tmp_out_bed {output.bed}
        cat {output.bed} | sort -k1,1 -k2,2n > $tmp_out_sort_bed
        mv $tmp_out_sort_bed {output.sorted_bed}
        """
