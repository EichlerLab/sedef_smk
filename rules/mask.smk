def get_RM_bed(wildcards):
    return manifest_df.loc[(wildcards.sample, wildcards.hap), "RM"]

def get_TRF_bed(wildcards):
    return manifest_df.loc[(wildcards.sample, wildcards.hap), "TRF"]

rule split_hap_bed:
    input:
        fai = rules.index_fasta.output.fai,
        rm_bed = get_RM_bed,
        trf_bed = get_TRF_bed,
    output:
        hap_rm_bed = "results/{sample}/work/split_beds/outputs/{hap}.RM.bed",
        hap_trf_out = "results/{sample}/work/split_beds/outputs/{hap}.trf.out", ## NOTE: The raw TRF BEDGZ is not a BED format. 1-base / 1-base coordinate
    threads: 1,
    resources:
        mem = 16,
        hrs = 24,
    run:
        def is_gzip(__path):
            with open(__path, "rb") as filehead:
                return filehead.read(2) == b"\x1f\x8b"

        def filter_bed_by_contigs(bed_in: str, bed_out: str, contigs: set[str]) -> None:
            bed_in_path = Path(bed_in)
            opener = gzip.open if is_gzip(bed_in_path) else open

            with opener(bed_in_path, "rt") as fin, open(bed_out, "w") as fout:
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
        # TRF out process
        filter_bed_by_contigs(input.trf_bed, output.hap_trf_out, hap_contigs)


rule sort_trf_bed:
    input:
        trf_out = rules.split_hap_bed.output.hap_trf_out,
    output:
        sorted_trf_bed = "results/{sample}/work/split_beds/outputs/{hap}.TRF.sorted.bed",
    threads: 1,
    resources:
        mem = 4,
        hrs = 4,
    shell: """
        tmp_output="{output.sorted_trf_bed}.tmp"
        grep -v '^#' {input.trf_out} | \
        cut -f1-3 | \
        awk '{{print $1"\t"$2-1"\t"$3}}' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - \
        > $tmp_output
        mv $tmp_output {output.sorted_trf_bed}
        """


rule get_final_bed:
    input:
        sorted_wm_bed = rules.run_windowmasker_bed.output.sorted_bed,
        rm_bed = rules.split_hap_bed.output.hap_rm_bed,
        sorted_trf_bed = rules.sort_trf_bed.output.sorted_trf_bed,
        hap_fai = rules.index_fasta.output.fai,
    output:
        sorted_rm_bed = "results/{sample}/work/split_beds/outputs/{hap}.RM.sorted.bed",
        merged_trf_bed = "results/{sample}/work/merge_beds/outputs/{hap}.TRF.merged.bed",
        trf_wm_bed = "results/{sample}/work/merge_beds/outputs/{hap}.TRF_WM.bed",
        merged_trf_wm_bed = "results/{sample}/work/merge_beds/outputs/{hap}.mTRF_WM.bed",
        final_bed = "results/{sample}/outputs/final_bed/{hap}.final_repeat.bed",
        flag = "results/{sample}/work/get_final_bed/flag/{hap}.final_repeat_bed.done"
    threads: 1,
    resources:
        mem = 16,
        hrs = 24,
    shell: """
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

        touch {output.flag}
        """

rule mask_fasta:
    input:
        fasta = rules.index_fasta.output.fasta,
        final_bed = rules.get_final_bed.output.final_bed
    output:
        masked_fasta = "results/{sample}/outputs/masked_fasta/{hap}.masked.fasta",
    threads: 1,
    resources:
        mem = 8,
        hrs = 4,
    singularity:
        "docker://eichlerlab/binf-basics:0.2"
    shell: """
        tmp_output="{output.masked_fasta}.tmp"
        seqtk seq -l 50 -M {input.final_bed} {input.fasta} > $tmp_output
        mv $tmp_output {output.masked_fasta}
        """


# #systemstr += " && echo '3' && cat "+repeatmask_out_file+" | sort -k1,1 -k2,2n >| "+tmp1
# systemstr += " && echo '4' && cat "+trf_out_file+" | cut -f 1-3 | awk '{if (\$3-\$2 > 5000) print \$0}' | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - | bedtools slop -b 100 -g "+species_genome_fasta+".fai -i - >| "+tmp2
# systemstr += " && echo '5' && cat "+trf_out_file+" "+windowmask_out_file+" | cut -f 1-3 >| "+tmp3
# systemstr += " && echo '6' && cat "+tmp2+" "+tmp3+" | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{if (\$3-\$2 > 5000) print \$0}' | bedtools merge -d 100 -i - >| "+tmp4
# systemstr += " && echo '7' && cat "+tmp1+" "+tmp2+" "+tmp3+" "+tmp4+" | cut -f 1-3 | sort -k1,1 -k2,2n| bedtools merge -i - >| "+final_out
# systemstr += " && module load seqtk/1.3 && seqtk seq -l 50 -M "+final_out+" "+species_genome_fasta+" > "+masked_fa+" && samtools faidx "+masked_fa
# systemstr += " && bgzip "+masked_fa+" && samtools faidx "+masked_fa+".gz"
# systemstr += '" | qsub -q eichler-short.q -N finalRM_'+species_prefix+' -pe serial 1 -l h_rt=5:00:00 -l mfree=12G'