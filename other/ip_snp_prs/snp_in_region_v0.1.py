import HTSeq

# snp_addr  = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_clump.txt'
snp_addr  = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_lmiss0.05_clump.clumped.txt'

# ip_gene_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/IP_and_InWeb_list.txt'
# out_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_clump_SNPinIP.txt'
# out_addr  = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_lmiss0.05_clump_SNPinIP.txt'

# 20180426  NC
# ip_gene_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/NC/outputFile.tsv'  #NC
# out_addr  = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/NC/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_lmiss0.05_clump_SNPinIP_NC.txt'

#20180503  NC2
# ip_gene_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/NC2/negative_controls.txt'
# out_addr  = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/NC2/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_lmiss0.05_clump_SNPinIP_NC2.txt'

# 20181119 NEW data
ip_gene_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/bine_20181119/bine_10percent.txt'
out_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/bine_20181119/my.1000GP_Phase3_chrAll.EAS.bfile.maf0.05_lmiss0.05_clump_SNPinIP_bine_10_percent.txt'

f_ip_gene = open(ip_gene_addr)
ip_header = f_ip_gene.next()
ip_genome = HTSeq.GenomicArrayOfSets("auto", stranded=False)

for line in f_ip_gene:
    l_line = line.strip().split()
    chr_id = l_line[4]
    s_pos, e_pos = l_line[5].split('-')#  [s,e]
    s_pos = int(s_pos) - 200e3
    e_pos = int(e_pos) + 200e3
    if s_pos < 0:s_pos = 0

    ip_gene_iv = HTSeq.GenomicInterval(chr_id, int(s_pos), int(e_pos) + 1)
    ip_genome[ip_gene_iv] += '1'

f_snp_clump = open(snp_addr)
snp_header = f_snp_clump.next()

result = []
for snp_line in f_snp_clump:
    if snp_line == '\n':
        continue
    l_snp_line = snp_line.strip().split()
    snp_chr = l_snp_line[0]
    snp_pos = int(l_snp_line[3])
    snp_id = l_snp_line[2]
    snp_iv = HTSeq.GenomicPosition(snp_chr, snp_pos)
    for value in ip_genome[snp_iv]:
        result.append(snp_id)

f_out = open(out_addr,'w')
f_out.write('\n'.join(result) + '\n')

f_out.close()



