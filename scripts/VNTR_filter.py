import sys,os

fout = open("tmp_SD_flanking.bed",'w')
with open("sedef_out/final.bed",'r') as fp:
    for line in fp:
        if line.startswith("#"):continue
        line_temp = line.strip().split('\t')
        #print (line)
        sname = ':'.join(line_temp[0:3])
        fout.write(line_temp[0]+'\t'+line_temp[1]+'\t'+str(int(line_temp[1])+1000)+'\t'+sname+'\n')
        fout.write(line_temp[3]+'\t'+line_temp[4]+'\t'+str(int(line_temp[4])+1000)+'\t'+sname+'\n')
        if int(line_temp[2]) > 1000:
            fout.write(line_temp[0]+'\t'+str(int(line_temp[2])-1000)+'\t'+line_temp[2]+'\t'+sname+'\n')
        else:
            fout.write(line_temp[0]+'\t'+str(0)+'\t'+line_temp[2]+'\t'+sname+'\n')
        if int(line_temp[5]) > 1000:
            fout.write(line_temp[3]+'\t'+str(int(line_temp[5])-1000)+'\t'+line_temp[5]+'\t'+sname+'\n')
        else:
            fout.write(line_temp[3]+'\t'+str(0)+'\t'+line_temp[5]+'\t'+sname+'\n')
fout.close()

fout = open("trf.10k.bed",'w')
with open(sys.argv[1],'r') as fp2: #"trf.sorted.merged.bed",'r') as fp2:
    for line in fp2:
        line_temp = line.strip().split('\t')
        if int(line_temp[2]) - int(line_temp[1]) > 10000:
            fout.write(line)
fout.close()
os.system("bedtools intersect -f 1 -a tmp_SD_flanking.bed -b trf.10k.bed -wa > tmp_SD_flanking.overlapped_VNTR.bed")

dd = {}
with open("tmp_SD_flanking.overlapped_VNTR.bed",'r') as fp3:
    for line in fp3:
        line_temp = line.strip().split('\t')
        dd[line_temp[3]] = ''

fout = open("sedef_out/final.rmVNTR.bed",'w')
with open("sedef_out/final.bed",'r') as fp:
    for line in fp:
        if line.startswith("#"):
            fout.write(line)
            continue
        line_temp = line.strip().split('\t')
        sname = ':'.join(line_temp[0:3])
        if sname in dd:continue
        fout.write(line)
fout.close()




        


'''
h2tg000001l     148537  149923  h2tg000001l     151944  153275  S       20.6    +       +       1386    1429    m=10.8;g=9.9    43      98      1288    1134    154     88      66      0.880435 0.793562 0.130248        0.13109 7       509     482     459     1134    154     7       141     933M16D39M26D50M56D80M30I34M1I40M11I63M1I49M    0.875676
h2tg000001l     154749  155780  h2tg000001l     167581  168615  S       16.7    +       +       1034    1055    m=12.4;g=4.3    24      21      1010    879     131     83      48      0.870297 0.833175 0.142406        0.144022        13      568     627     461     879     131     13      45      19M4I74M1D60M1D50M2D50M11D360M2D146M6I52M4D36M7I10M1I46M2I42M1I33M3I32M 0.859238
h2tg000001l     183579  185309  h2tg000001l     11526201        11527951        S       10.1    +       +       1750    1754    m=8.6;g=1.6     24      4       1726    1576    150     102     48
        0.913094        0.898518        0.0923675       0.09322 11      969     1083    891     1576    150     11      28      32M1D119M10I89M1D378M1I146M3I155M1D170M1D190M3I69M2I172M3I94M2I112M       0.907311
bedtools subtract -A -a HG00733_hap1.SDs.bed -b <(cat trf.sorted.merged.bed | awk '{if ($3-$2 > 10000) print $0}') | bedtools merge -i - | rb bed-length -r
'''
