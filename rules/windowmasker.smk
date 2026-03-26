def get_fasta(wildcards):
    return asm_manifest_df.loc[(wildcards.sample, wildcards.hap), "FASTA"]



rule index_fasta:
    input:
        fasta = get_fasta
    output:
        fasta = "results/{sample}/work/windowmasker/inputs/{hap}.fasta",
        fai = "results/{sample}/work/windowmasker/inputs/{hap}.fasta.fai"
    threads: 1
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        set -euo pipefail
        # to prevent errors in sedef.sh due to "#" in the contig names.
        sed '/^>/ s/#/__/g' "$(readlink -f {input.fasta})" > {output.fasta}
        samtools faidx {output.fasta}
        """

rule make_dust_count:
    input:
        ref = rules.index_fasta.output.fasta,
    output:
        count_file = "results/{sample}/work/windowmasker/intermediates/{hap}.count",
        flag = touch("results/{sample}/work/windowmasker/flags/dust_count.{hap}.done"),
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
        set -euo pipefail
        windowmasker \
        -mem {params.mem} \
        -smem {params.smem} \
        -mk_counts \
        -infmt fasta \
        -sformat obinary \
        -in {input.ref} \
        -out {output.count_file}
        """

rule run_windowmasker:
    input:
        count_file = rules.make_dust_count.output.count_file,
        ref = rules.index_fasta.output.fasta,
        fai = rules.index_fasta.output.fai,
    output:
        intervals = "results/{sample}/work/windowmasker/intermediates/{hap}.intervals",
        flag = touch("results/{sample}/work/windowmasker/flags/run_windowmasker.{hap}.done")
    resources:
        mem = 16,
        hrs = 12,
    threads: 1
    singularity:
        "docker://eichlerlab/ncbi-windowmasker:1.0"
    shell: """
        set -euo pipefail
        windowmasker \
        -infmt fasta \
        -ustat {input.count_file} \
        -dust T \
        -outfmt interval \
        -in {input.ref} \
        -out {output.intervals}
        """

rule run_windowmasker_to_bed:
    input:
        intervals = rules.run_windowmasker.output.intervals,
    output:
        bed = "results/{sample}/work/windowmasker/intermediates/{hap}.dust.bed",
        sorted_bed = "results/{sample}/work/windowmasker/outputs/{hap}.mask_repeat.sorted.bed",
        flag = touch("results/{sample}/work/windowmasker/flags/run_windowmasker_to_bed.{hap}.done")
    resources:
        mem=16,
        hrs=12,
    params:
        script_dir=f"{SNAKEMAKE_DIR}/scripts"
    threads: 1
    singularity:
        "docker://eichlerlab/ncbi-windowmasker:1.0"
    shell: """
        set -euo pipefail
        perl {params.script_dir}/windowmasker_to_bed.pl < {input.intervals} > {output.bed}
        cat {output.bed} | sort -k1,1 -k2,2n > {output.sorted_bed}
        """
