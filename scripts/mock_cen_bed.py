import sys
dd = {}
with open("sedef_out/cen.bed",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        dd[line_temp[0]] = line

fout = open("sedef_out/cen.mock.bed",'w')
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        if line_temp[0] in dd:
            fout.write(dd[line_temp[0]])
        else:
            fout.write(line_temp[0]+'\t0\t1\n')
fout.close()
'''
head genome.masked.1M.fa.fai
chr1_RagTag     221534370       13      50      51
chr10_RagTag    142678062       225965085       50      51
chr11_RagTag    132458108       371496723       50      51
chr12_RagTag    127527147       506604008       50      51
chr13_RagTag    123422281       636681712       50      51
chr14_RagTag    118048196       762572453       50      51
chr15_RagTag    104730676       882981627       50      51
chr16_RagTag    100538858       989806931       50      51
chr17_RagTag    79415468        1092356581      50      51
chr18_RagTag    51014245        1173360373      50      51
 hsjeong@ocelot  .../segdup/SEDEF_hap1  head -n 3 sedef_out/cen.bed
chr10_RagTag    27636988        29855310
chr11_RagTag    35915935        36824606
'''
