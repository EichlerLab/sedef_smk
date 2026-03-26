rule sedef:
    input:
        masked_fasta = rules.mask_fasta.output.masked_fasta,
    output:
        sedef_final_bed = "results/{sample}/work/sedef/outputs/{hap}.sedef.bed",
        flag = touch("results/{sample}/work/sedef/flags/sedef.{hap}.done")
    threads: 32
    resources:
        hrs=72,
        mem=6,
    singularity: 
        "docker://eichlerlab/sedef:1.1"
    shell: """
        set -euo pipefail
        tmp_out="{resources.tmpdir}/{wildcards.sample}-{wildcards.hap}-sedef_out"
        sedef.sh -j {threads} -f -o $tmp_out {input.masked_fasta}
        rsync -av $tmp_out/final.bed {output.sedef_final_bed}
    """

rule vntr_filter:
    input:
        sedef_bed = rules.sedef.output.sedef_final_bed,
        sorted_trf_bed = rules.sort_trf_bed.output.sorted_trf_bed,
    output:
        sd_flanking_bed = "results/{sample}/work/sedef_filter/intermediates/{hap}.sd_flanking.bed",
        sd_flanking_vntr_bed = "results/{sample}/work/sedef_filter/intermediates/{hap}.sd_flanking.vntr.bed",
        trf_10k_bed = "results/{sample}/work/sedef_filter/intermediates/{hap}.trf_10k.bed",
        sedef_rm_vntr_bed = "results/{sample}/work/sedef_filter/intermediates/{hap}.sedef.rm_vntr.bed",
        flag = touch("results/{sample}/work/sedef_filter/flags/vntr_filter.{hap}.done")
    threads: 1
    resources:
        hrs=12,
        mem=12,
    script: f"{SNAKEMAKE_DIR}/scripts/VNTR_filter.py"

rule check_sat:
    input:
        sorted_rm_bed = rules.get_final_bed.output.sorted_rm_bed,
        sedef_rm_vntr_bed = rules.vntr_filter.output.sedef_rm_vntr_bed
    output:
        sat_count_bed = "results/{sample}/work/sedef_filter/outputs/{hap}.sat_count.bed",
        flag = touch("results/{sample}/work/sedef_filter/flags/check_sat.{hap}.done"),
    threads: 1
    resources:
        hrs=12,
        mem=12,
    shell: """
        set -euo pipefail
        grep "^#" {input.sedef_rm_vntr_bed} > {output.sat_count_bed}
        grep -e "Satellite" {input.sorted_rm_bed} \
            | cut -f 1-3 \
            | bedtools sort -i - \
            | bedtools merge -i - \
            | bedtools coverage -header -a {input.sedef_rm_vntr_bed} -b - \
            >> {output.sat_count_bed}
        sed -i '1{{s/$/\tcount_ovls\tsat_bases\ttotal_bases\tsat_coverage/}}' {output.sat_count_bed}
    """

rule check_cen:
    input:
        sorted_rm_bed = rules.get_final_bed.output.sorted_rm_bed,
        sedef_rm_vntr_bed = rules.vntr_filter.output.sedef_rm_vntr_bed
    output:
        cen_count_bed = "results/{sample}/work/sedef_filter/intermediates/{hap}.cen_count.bed",
        flag = touch("results/{sample}/work/sedef_filter/flags/check_cen.{hap}.done")
    threads: 1
    resources:
        hrs=12,
        mem=12,
    shell: """
        set -euo pipefail
        grep -e 'ALR/Alpha' -e "centr" -e "Centr" -e "Centrom" {input.sorted_rm_bed} \
        | bedtools merge -d 100 -i - \
        | awk '($3-$2)>100' \
        | bedtools merge -d 2000 -i - \
        | awk '($3-$2)>10000' \
        | bedtools merge -d 200000 -i - \
        | awk '{{print $0"\t"$3-$2}}' \
        | sort -k 1,1 -k4,4n \
        | bedtools groupby -g 1 -c 2,3 -o last,last \
        > {output.cen_count_bed}
    """

rule mock_cen_bed:
    input:
        cen_count_bed = rules.check_cen.output.cen_count_bed,
        masked_fai = rules.mask_fasta.output.masked_fai,
    output:
        cen_mock_bed = "results/{sample}/work/sedef_filter/outputs/{hap}.cen_mock.bed",
        flag = touch("results/{sample}/work/sedef_filter/flags/mock_cen_bed.{hap}.done")
    threads: 1
    resources:
        hrs=4,
        mem=6,
    run:
        dd = {}
        with open(input.cen_count_bed,'r') as fp:
            for line in fp:
                line_temp = line.strip().split('\t')
                dd[line_temp[0]] = line

        fout = open(output.cen_mock_bed,'w')
        with open(input.masked_fai,'r') as fp:
            for line in fp:
                line_temp = line.strip().split('\t')
                if line_temp[0] in dd:
                    fout.write(dd[line_temp[0]])
                else:
                    fout.write(line_temp[0]+'\t0\t1\n')
        fout.close()

rule sedef_to_bed:
    input:
        cen_mock_bed = rules.mock_cen_bed.output.cen_mock_bed,
        sat_count_bed = rules.check_sat.output.sat_count_bed,
        masked_fai = rules.mask_fasta.output.masked_fai,
    output:
        sd_bed = "results/{sample}/work/get_final_beds/intermediates/{hap}.SDs.renamed.bed",
        sd_filt_bed = "results/{sample}/work/get_final_beds/intermediates/{hap}.SDs.renamed.lowid.bed",
        flag = touch("results/{sample}/work/get_final_beds/flags/sedef_to_bed.{hap}.done"),
    params:
        script_dir=f"{SNAKEMAKE_DIR}/scripts",
        sat=0.70,
        peri=1,
        telo=1,
    threads: 1
    resources:
        hrs=12,
        mem=12,
    shell: """
        set -euo pipefail
        python {params.script_dir}/sedef_to_bed.py \
            --fai {input.masked_fai} \
            --cens {input.cen_mock_bed} \
            --sat {params.sat} \
            --peri {params.peri} \
            --telo {params.telo} \
            {input.sat_count_bed} \
            {output.sd_bed} \
            {output.sd_filt_bed}
        """

rule restore_bed:
    input:
        sd_bed = rules.sedef_to_bed.output.sd_bed,
        sd_filt_bed = rules.sedef_to_bed.output.sd_filt_bed
    output:
        sd_bed = "results/{sample}/final_outputs/beds/{hap}.SDs.bed",
        sd_filt_bed = "results/{sample}/final_outputs/beds/{hap}.SDs.lowid.bed",
        flag = touch("results/{sample}/work/get_final_beds/flags/restore_bed.{hap}.done")
    threads: 1,
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        set -euo pipefail
        sed 's/__/#/g' {input.sd_bed} > {output.sd_bed}
        sed 's/__/#/g' {input.sd_filt_bed} > {output.sd_filt_bed}
    """


rule make_bigbeds:
    input:
        masked_fai = rules.restore_fasta.output.masked_fai,
        sd_bed = rules.restore_bed.output.sd_bed,
        sd_filt_bed = rules.restore_bed.output.sd_filt_bed,
    output:
        sd_bigbed = "results/{sample}/final_outputs/bigbeds/{hap}.SDs.bb",
        sd_filt_bigbed = "results/{sample}/final_outputs/bigbeds/{hap}.SDs.lowid.bb",
        flag = touch("results/{sample}/work/get_final_beds/flags/{hap}.make_bigbeds.done")
    params:
        schema = f"{SNAKEMAKE_DIR}/schema/sedef.as"
    threads: 1
    resources:
        hrs=12,
        mem=12,
    singularity: "docker://eichlerlab/ucsc:202405"
    shell: """
        bedToBigBed -type=bed9+32 -tab -as={params.schema} {input.sd_bed} {input.masked_fai} {output.sd_bigbed}
        bedToBigBed -type=bed9+32 -tab -as={params.schema} {input.sd_filt_bed} {input.masked_fai} {output.sd_filt_bigbed}
        """