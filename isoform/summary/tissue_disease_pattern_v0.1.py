# find the relation between the FM exon-skipping and FM disease

import numpy as np
import itertools
import pandas as pd
import sys, os

FM_skipping_disease_addr = sys.argv[1]
out_addr = sys.argv[2]

f_FM_skipping = open(FM_skipping_disease_addr)

# tissue: [3 : 37] # 34
#disease pval: [43 : 67] # 24

def count_FM_exon_disease(sub):
    # sub : [l_line1, l_line2,...]
    sub_res_array = np.zeros((34, 24))
    for l_line in sub:
        tissue_l = l_line[3:37]
        disease_pval_l = l_line[43:67]
        i_l = [t for t in range(tissue_l) if t != '-1']
        j_l = [d for d in range(disease_pval_l) if d != '-1']
        for index in itertools.product(i_l, j_l):
            sub_res_array[index] = 1
    sub_res_array[sub_res_array > 0] = 1
    return sub_res_array

res_array = np.zeros((34, 24))

old_id = ''
sub = []
header = f_FM_skipping.next()
header_l = header.split()
tissue_header_l = header_l[3:37]
disease_header_ = [each.split('_')[1] for each in header_l[43:67]]

for line in f_FM_skipping:
    l_line = line.split()
    type = l_line[2]
    gene_id = l_line[0]
    if type == 'p':
        if old_id == '':
            old_id = gene_id
            sub = [l_line]
        elif old_id == gene_id:
            sub.append(l_line)
        else:
            sub_res_array = count_FM_exon_disease(sub)
            res_array += sub_res_array
            old_id = gene_id
            sub = [l_line]
else:
    res_array += count_FM_exon_disease(sub)

df = pd.DataFrame(res_array)
df.to_csv(out_addr, sep='\t')

