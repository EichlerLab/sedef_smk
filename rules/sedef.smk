rule sedef:
    input:
        masked_fasta = rules.mask_fasta.output.masked_fasta,
    output:
        sedef_out = directory("results/{sample}/work/sedef/final_beds/{hap}"),
        flag = "results/{sample}/work/sedef/flags/{hap}.sedef_done"
    threads: 32,
    resources:
        hrs=72,
        mem=2,
    singularity: 
        "docker://eichlerlab/sedef:1.1"
    shell: """
        set -euo pipefail
        tmp_out="{resources.tmpdir}/{wildcards.sample}-{wildcards.hap}-sedef_out"
        sedef.sh -j {threads} -f -o $tmp_out {input.masked_fasta}
        rsync -av $tmp_out/final* {output.sedef_out}/
        touch {output.flag}
    """









# module load ucsc/202403
# module load bedtools/2.31.1
# module load samtools/1.17 htslib/1.19 bcftools/1.20
# module load miniconda/4.12.0

# . ~/.bashrc

# #cp /net/eichler/vol27/projects/marmoset/nobackups/segdup/sedef_inputs/CJ1700_CJA_hap1/* ./
# mkdir -p sedef_out
# mkdir -p /tmp/sedef_input_DSA_COLO829BL_v3.hap1 && mkdir -p /tmp/sedef_output_DSA_COLO829BL_v3.hap1 && mkdir -p sedef_out && zcat DSA_COLO829BL_v3.hap1.masked.fa.gz > /tmp/sedef_input_DSA_COLO829BL_v3.hap1/DSA_COLO829BL_v3.hap1.masked.fa
# sedef.sh -j 50 -f -o /tmp/sedef_output_DSA_COLO829BL_v3.hap1 /tmp/sedef_input_DSA_COLO829BL_v3.hap1/DSA_COLO829BL_v3.hap1.masked.fa

# cp -r /tmp/sedef_output_DSA_COLO829BL_v3.hap1/final* ./sedef_out/
# rm -rf /tmp/sedef_input_DSA_COLO829BL_v3.hap1 /tmp/sedef_output_DSA_COLO829BL_v3.hap1

# python VNTR_filter.py DSA_COLO829BL_v3.hap1.trf.sorted.merged.p2k.bed
# cat <(cat DSA_COLO829BL_v3.hap1.full_mask_repeat.sorted.bed | grep -e "Satellite") | cut -f 1-3 | bedtools sort -i - | bedtools merge -i - | bedtools coverage -header -a sedef_out/final.rmVNTR.bed -b - > sedef_out/tmp_final.sat.count.bed
# grep "^#" sedef_out/final.rmVNTR.bed > sedef_out/final.sat.count.bed
# cat sedef_out/tmp_final.sat.count.bed >> sedef_out/final.sat.count.bed
# sed -i '1{{s/$/\tcount_ovls\tsat_bases\ttotal_bases\tsat_coverage/}}' sedef_out/final.sat.count.bed
# grep -e 'ALR/Alpha' -e "centr" -e "Centr" -e "Centrom" DSA_COLO829BL_v3.hap1.full_mask_repeat.sorted.bed | bedtools merge -d 100 -i - |awk '($3-$2)>100' |bedtools merge -d 2000 -i - |awk '($3-$2)>10000' |bedtools merge -d 200000 -i - |awk '{{print $0"\t"$3-$2}}' |sort -k 1,1 -k4,4n > sedef_out/tmp.cen.bed
# bedtools groupby -g 1 -c 2,3 -o last,last -i sedef_out/tmp.cen.bed > sedef_out/cen.bed
# python mock_cen_bed.py DSA_COLO829BL_v3.hap1.masked.fa.gz.fai


# ./sedef_to_bed.py --fai DSA_COLO829BL_v3.hap1.masked.fa.gz.fai --cens sedef_out/cen.mock.bed --sat 0.70 --peri 1 --telo 1 sedef_out/final.sat.count.bed DSA_COLO829BL_v3.hap1.SDs.bed DSA_COLO829BL_v3.hap1.SDs.lowid.bed
# #bedToBigBed -type=bed9+32 -tab -as=sedef.as DSA_COLO829BL_v3.hap1.SDs.bed DSA_COLO829BL_v3.hap1.masked.fa.gz.fai DSA_COLO829BL_v3.hap1.SDs.bb
# #bedToBigBed -type=bed9+32 -tab -as=sedef.as DSA_COLO829BL_v3.hap1.SDs.lowid.bed DSA_COLO829BL_v3.hap1.masked.fa.gz.fai DSA_COLO829BL_v3.hap1.SDs.lowid.bb