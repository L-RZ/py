import numpy as np
import pandas as pd
import HTSeq, os,sys,time

#v0.1   calculate PRS by different region from different gene groups
#v0.2   get the OR from a file which created by fix_OR.py


## for test

gt_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/test_dir/raw.txt'

or_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/test_dir/logOR.txt'

ip_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/test_dir/IP_InWeb.txt'

map_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/test_dir/tmp1.map'

out_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/test_dir/test_result1.txt'

#v0.1
# gt_addr = sys.argv[1]
# or_addr = sys.argv[2]
# ip_addr = sys.argv[3]
# map_addr = sys.argv[4]
# out_addr = sys.argv[5]

#v0.2
# gt_addr = sys.argv[1]
# or_addr = sys.argv[2]
# ip_addr = sys.argv[3]
# map_addr = sys.argv[4]
# out_addr = sys.argv[5]


def make_snp_genome(map_addr):
    f_map = open(map_addr)
    snp_genome = HTSeq.GenomicArrayOfSets("auto",stranded=False)
    all_snp_list= []
    for map_line in f_map:
        map_line_l = map_line.strip().split()
        map_chr = map_line_l[0]
        map_pos = int(map_line_l[3])
        map_snp_id = map_line_l[1]
        snp_iv = HTSeq.GenomicInterval(map_chr, map_pos, map_pos + 1)
        snp_genome[snp_iv] += map_snp_id
        all_snp_list.append(map_snp_id)

    return snp_genome, all_snp_list


def rc(seq):
    if '<' in seq:
        return seq
    c = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '+': '+'}
    if '+' in seq:
        seq_l = seq.split('+')
        seq = seq_l[0]
        ans = ''.join(map(lambda x: c[x], list(seq)))[::-1] + '+%s' % seq_l[1]
    else:
        ans = ''.join(map(lambda x: c[x], list(seq)))[::-1]
    return ans

def make_logOR(assoc_addr, gt_addr, all_snp_list):
    f_assoc = open(assoc_addr )
    assoc_header = f_assoc.next()
    assoc_dict = {}
    for line in f_assoc:
        line_l = line.split()

        snp_id = line_l[1]
        snp_A1 = line_l[3]
        snp_A2 = line_l[4]
        snp_frq = line_l[5]
        snp_OR  = line_l[8]

        # if snp_id in all_snp_list:
        #     assoc_dict[snp_id] = (snp_A1, snp_A2, snp_frq, snp_OR)
        assoc_dict[snp_id] = ' '.join([snp_A1, snp_A2, snp_frq, snp_OR])

    f_gt = open(gt_addr)
    gt_header = f_gt.next()
    gt_header_l = gt_header.split()
    logor_list = []
    result = ['\t'.join(['snp_id', 'A1', 'logOR(A1)','frq(A1)'])]
    # print gt_header_l   # test
    for gt_snp_id in gt_header_l[6:]:
        snp_id, allele = gt_snp_id.split('_')
        snp_A1, snp_A2, snp_frq, snp_OR = assoc_dict[snp_id]
        snp_logOR = np.log(float(snp_OR))
        snp_frq = float(snp_frq)
        A_T_C_G = [set(('A','T')), set(('C','G'))]
        if set((snp_A1, snp_A2)) in A_T_C_G:
            if float(snp_frq) > 0.5:   #A1 frq

                new_A1 = allele
                new_logOR = -snp_logOR
                new_frq = 1 - snp_frq
            else:
                new_A1 = allele
                new_logOR = snp_logOR
                new_frq = snp_frq
        else:
            if allele not in [snp_A1, rc(snp_A1)]:
                if allele not in [snp_A2, rc(snp_A2)]:
                    print 'error:', gt_snp_id
                else:
                    new_A1 = allele
                    new_logOR = -snp_logOR
                    new_frq = 1 - snp_frq
            else:
                new_A1 = allele
                new_logOR = snp_logOR
                new_frq = snp_frq
        logor_list.append(snp_logOR)
        result.append('\t'.join([snp_id, new_A1, str(new_logOR), str(new_frq)]))

    print '\n'.join(result)                                 # test

    return logor_list

def get_logOR(or_addr):
    f_logOR = open(or_addr)
    header = f_logOR.next()
    logor_list = [line.strip().split()[2] for line in f_logOR]
    return logor_list


print 'start',time.ctime()

snp_genome, all_snp_list = make_snp_genome(map_addr)
# logOR = make_logOR(assoc_addr, gt_addr, all_snp_list)

logOR = get_logOR(or_addr)
logOR = np.array(logOR, dtype=float)

gt_array = []
f_gt = open(gt_addr)
gt_header = f_gt.next()

phe_list = []
gt_list  = []
#

for gt_line in f_gt:
    gt_line_l = gt_line.strip().split()
    phe = gt_line_l[5]
    gt = np.array(gt_line_l[6:], dtype='float16')
    phe_list.append(phe)
    gt_list.append(gt)

df = pd.read_table(ip_addr)
bait_list = list(set(df['bait']))
bait_list.sort()

f_out = open(out_addr, 'w')
f_out.write('\t'.join(['IP\tbait\tint\tn_SNP\tsum_gene_len'] + phe_list) + '\n')
print 'done gt',time.ctime()
for bait in bait_list:
    bait_df = df[df['bait'] == bait]
    ip_list = list(set(bait_df['IP']))
    ip_list.sort()
    for ip in ip_list:
        ip_df = bait_df[bait_df['IP'] == ip]  # bait type : all_combined
        int_type = list(set(ip_df['int']))
        int_type.sort()
        for i in int_type:  # y , n
            ip_int_df = ip_df[ip_df['int'] == i]
            # ip = ' '.join(list(set(bait_int_df['IP'])))
            gene_len = 0
            snp_used = np.zeros(len(all_snp_list), dtype=bool)
            for index, row in ip_int_df.iterrows():
                snp_used_id_list = []
                gene_chr = str(row['chr'])
                gene_s, gene_e = row['pos'].split('-')
                gene_s = int(gene_s)
                gene_e = int(gene_e)

                s = gene_s - 200e3
                if s <= 0: s = 1
                e = gene_e + 200e3
                gene_iv = HTSeq.GenomicInterval(gene_chr, s, e + 1)
                gene_len += gene_e - gene_s + 1
                for iv, val in snp_genome[gene_iv].steps():
                    if len(val) != 0:
                        snp_id = list(val)[0]  # snp id
                        snp_used_id_list.append(snp_id)
                        snp_id_index = all_snp_list.index(snp_id)
                        snp_used[snp_id_index] = True
            snp_num = np.sum(snp_used)

            # for PRS
            iid_prs_list = []

            for gt in gt_list:
                gt_used = gt[snp_used]
                logOR_used = logOR[snp_used]

                iid_prs = np.sum(gt_used * logOR_used) / snp_num / gene_len * 1e6  # prs /(Mbp*SNP)
                # iid_prs = np.sum(gt_used * logOR_used) / snp_num
                iid_prs_list.append(str(iid_prs))

                print 'gt:', gt   # for test
                print 'snp_used:', snp_used
                print 'logOR_used:', logOR_used
            print ip, bait, i, snp_num, gene_len
            out_line = [ip, bait, i, str(snp_num), str(gene_len)] + iid_prs_list
            f_out.write('\t'.join(out_line) + '\n')

print 'Done', time.ctime()


