import numpy as np
import sys,os
import itertools, time
def rc(seq):
    if '<' in seq:
        return seq

    c = {'A':'T','T':'A','G':'C','C':'G','+': '+'}
    if '+' in seq:
        seq_l = seq.split('+')
        seq = seq_l[0]
        ans = ''.join(map(lambda x: c[x], list(seq)))[::-1] + '+%s'%seq_l[1]
    else:
        ans = ''.join(map(lambda x: c[x], list(seq)))[::-1]
    return ans

def make_logOR(assoc_addr, gt_addr, all_snp_list, out_addr):
    f_assoc = open(assoc_addr)
    assoc_header = f_assoc.next()
    assoc_dict = {}
    all_snp_set = set(all_snp_list)
    for line in f_assoc:
        line_l = line.split()
        snp_chr = line_l[0]
        snp_id = line_l[1]
        snp_A1 = line_l[3]
        snp_A2 = line_l[4]
        snp_frq = line_l[5]
        snp_OR  = line_l[8]

        assoc_dict[snp_id] = ' '.join([snp_A1, snp_A2, snp_frq, snp_OR])

        # if snp_id in all_snp_set:
        #     if snp_chr in assoc_dict:
        #         assoc_dict[snp_chr][snp_id] = ' '.join([snp_A1, snp_A2, snp_frq, snp_OR])
        #     else:
        #         assoc_dict[snp_chr] ={}
    print 'Done assoc_dict', time.ctime()
    f_gt = open(gt_addr)
    gt_header = f_gt.next()
    gt_header_l = gt_header.split()
    logor_list = []
    f_out = open(out_addr, 'w')
    f_out.write('\t'.join(['snp_id', 'A1', 'logOR(A1)','frq(A1)']) + '\n')
    # result = ['\t'.join(['snp_id', 'A1', 'logOR(A1)','frq(A1)'])]
    # print gt_header_l   # test
    # for gt_snp_id, gt_snp_chr in itertools.izip(gt_header_l[6:] ):
    for gt_snp_id in gt_header_l[6:]:
        snp_id, allele = gt_snp_id.split('_')
        snp_A1, snp_A2, snp_frq, snp_OR = assoc_dict[snp_id].split(' ')
        # snp_A1, snp_A2, snp_frq, snp_OR = assoc_dict[gt_snp_chr][snp_id].split(' ')
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
        # result.append('\t'.join([snp_id, new_A1, str(new_logOR), str(new_frq)]))
        f_out.write('\t'.join([snp_id, new_A1, str(new_logOR), str(new_frq)]) + '\n')
    # print '\n'.join(result)                                 # test


print ' '.join(sys.argv)
gt_addr = sys.argv[1]
assoc_addr= sys.argv[2]
map_addr = sys.argv[3]
out_addr = sys.argv[4]
print 'Start', time.ctime()
f_map = open(map_addr)
all_snp_list = []
all_snp_chr_list = []
for map_line in f_map:
    map_line_l = map_line.strip().split()
    all_snp_list.append(map_line_l[1])
    # all_snp_chr_list.append(map_line_l[0])

make_logOR(assoc_addr, gt_addr, all_snp_list, out_addr)
