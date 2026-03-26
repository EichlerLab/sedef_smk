def get_RM_bed(wildcards):
    return sample_manifest_df.loc[wildcards.sample, "RM"]

def get_trf_bed(wildcards):
    return sample_manifest_df.loc[wildcards.sample, "TRF"]

# The rule refine_beds is used to aviod using "#" in the query name which causes errors in sedef.
rule refine_beds:
    input:
        rm_bed = get_RM_bed,
        trf_bed = get_trf_bed,
    output:
        rm_bed = "results/{sample}/work/bed_masking/inputs/RM.refined.bed",
        trf_bed = "results/{sample}/work/bed_masking/inputs/trf.refined.bed",
        flag = touch("results/{sample}/work/bed_masking/flags/refine_beds.done")
    threads: 1,
    resources:
        mem=16,
        hrs=4,
    shell: """
        set -euo pipefail

        if [[ "{input.rm_bed}" == *.gz ]]; then
            gzip -cd "{input.rm_bed}" | sed 's/#/__/g' > {output.rm_bed}
        else
            sed 's/#/__/g' "{input.rm_bed}" > {output.rm_bed}
        fi

        if [[ "{input.trf_bed}" == *.gz ]]; then
            gzip -cd "{input.trf_bed}" | sed 's/#/__/g' > {output.trf_bed}
        else
            sed 's/#/__/g' "{input.trf_bed}" > {output.trf_bed}
        fi
    """

rule split_hap_bed:
    input:
        fai = rules.index_fasta.output.fai,
        rm_bed = rules.refine_beds.output.rm_bed,
        trf_bed = rules.refine_beds.output.trf_bed,
    output:
        hap_rm_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.RM.bed",
        hap_trf_out = "results/{sample}/work/bed_masking/intermediates/{hap}.trf.out", ## NOTE: The raw trf BEDGZ is not a BED format. 1-base / 1-base coordinate
        flag = touch("results/{sample}/work/bed_masking/flags/split_hap_bed.{hap}.done")
    threads: 1,
    resources:
        mem = 16,
        hrs = 24,
    run:
        def filter_bed_by_contigs(bed_in: str, bed_out: str, contigs: set[str]) -> None:
            bed_in_path = Path(bed_in)
            with open(bed_in_path) as fin, open(bed_out, "w") as fout:
                for line in fin:
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        fout.write(line)
                        continue
                    chrom = line.split("\t", 1)[0]
                    if chrom in contigs:
                        fout.write(line)

        hap_contigs = set()
        with open(input.fai) as fai:
            for line in fai:
                hap_contigs.add(line.strip().split("\t", 1)[0])


        # RM BED process
        filter_bed_by_contigs(input.rm_bed, output.hap_rm_bed, hap_contigs)
        # trf out process
        filter_bed_by_contigs(input.trf_bed, output.hap_trf_out, hap_contigs)


rule sort_trf_bed:
    input:
        trf_out = rules.split_hap_bed.output.hap_trf_out,
    output:
        sorted_trf_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.trf.sorted.bed",
        flag = touch("results/{sample}/work/bed_masking/flags/sort_trf_bed.{hap}.done")
    threads: 1,
    resources:
        mem = 4,
        hrs = 4,
    shell: """
        set -euo pipefail
        grep -v '^#' {input.trf_out} | \
        cut -f1-3 | \
        awk '{{print $1"\t"$2-1"\t"$3}}' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - \
        > {output.sorted_trf_bed}
        """


rule get_final_bed:
    input:
        sorted_wm_bed = rules.run_windowmasker_to_bed.output.sorted_bed,
        rm_bed = rules.split_hap_bed.output.hap_rm_bed,
        sorted_trf_bed = rules.sort_trf_bed.output.sorted_trf_bed,
        hap_fai = rules.index_fasta.output.fai,
    output:
        sorted_rm_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.RM.sorted.bed",
        merged_trf_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.trf.merged.bed",
        trf_wm_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.trf_WM.bed",
        merged_trf_wm_bed = "results/{sample}/work/bed_masking/intermediates/{hap}.merged_trf_WM.bed",
        final_bed = "results/{sample}/work/bed_masking/outputs/{hap}.final_repeat.bed",
        flag = touch("results/{sample}/work/bed_masking/flags/get_final_bed.{hap}.done")
    threads: 1,
    resources:
        mem = 16,
        hrs = 24,
    shell: """
        set -euo pipefail
        rm_out="{input.rm_bed}"
        trf_out="{input.sorted_trf_bed}"
        wm_out="{input.sorted_wm_bed}"
        tmp1="{output.sorted_rm_bed}"
        tmp2="{output.merged_trf_bed}"
        tmp3="{output.trf_wm_bed}"
        tmp4="{output.merged_trf_wm_bed}"

        # step 3
        sort -k1,1 -k2,2n $rm_out \
        > $tmp1

        # step 4
        awk '{{if ($3-$2 > 5000) print $0}}' $trf_out | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d 100 -i - | \
        bedtools slop -b 100 -g {input.hap_fai} -i - \
        > $tmp2

        # step 5
        cat $trf_out $wm_out | \
        cut -f1-3 \
        > $tmp3

        # step 6 (tmp4)
        cat $tmp2 $tmp3 | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - | \
        awk '{{if ($3-$2 > 5000) print $0}}' | \
        bedtools merge -d 100 -i - \
        > $tmp4

        # step 7 (final bed)
        cat $tmp1 $tmp2 $tmp3 $tmp4 | \
        cut -f1-3 | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - \
        > {output.final_bed}
        """

rule mask_fasta:
    input:
        fasta = rules.index_fasta.output.fasta,
        final_bed = rules.get_final_bed.output.final_bed
    output:
        masked_fasta = "results/{sample}/work/mask_fasta/intermediates/{hap}.masked.contig_renamed.fasta",
        masked_fai = "results/{sample}/work/mask_fasta/intermediates/{hap}.masked.contig_renamed.fasta.fai",
    threads: 1,
    resources:
        mem = 12,
        hrs = 4,
    singularity:
        "docker://eichlerlab/binf-basics:0.2"
    shell: """
        set -euo pipefail
        seqtk seq -l 50 -M {input.final_bed} {input.fasta} > {output.masked_fasta}
        samtools faidx {output.masked_fasta}
        """

rule restore_fasta:
    input:
        masked_fasta = rules.mask_fasta.output.masked_fasta
    output:
        masked_fasta = "results/{sample}/work/mask_fasta/outputs/{hap}.masked.fasta",
        masked_fai = "results/{sample}/work/mask_fasta/outputs/{hap}.masked.fasta.fai",
    threads: 1,
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        set -euo pipefail
        sed '/^>/ s/__/#/g' {input.masked_fasta} > {output.masked_fasta}
        samtools faidx {output.masked_fasta}
    """
