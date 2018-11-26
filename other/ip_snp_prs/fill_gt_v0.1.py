
import sys
import numpy as np
print sys.argv

gt_addr = sys.argv[1]
out_addr = sys.argv[2]


# gt_df = pd.read_table(gt_addr)
# gt_df_filling = gt_df.fillna(gt_df.mean())
#
# pd.to_csv(out_addr, sep=' ',index =False)


f_gt = open(gt_addr)
header = f_gt.next()
n_snp = len(header.split()[6:])
n=0
sum_array = np.zeros(n_snp)
for line in f_gt:
    l_line = line.strip().split()
    gt = np.array(['0' if i=='NA' else i for i in l_line[6:]],dtype=float)
    sum_array += gt
    n += 1
    if n%1000 == 0:print n

gt_mean = sum_array / n

f_gt = open(gt_addr)

f_out = open(out_addr, 'w')
header = f_gt.next()
f_out.write(header)
for line in f_gt:
    l_line = line.strip().split()
    for i, word in enumerate(l_line[6:]):
        if word == 'NA':l_line[6+i] = str(gt_mean[i])
    f_out.write(' '.join(l_line) + '\n')
print 'Done'

