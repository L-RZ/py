import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import pandas as pd
import numpy as np

in_addr = '/Users/ruizeliu/Documents/lrz/other/sczasn_IP/result/IP_and_InWeb_scz_asiz_risk_score.txt'
f_in = open(in_addr)
header = f_in.readline()
header_l = header.split()
df = pd.read_table(in_addr,skiprows=1,header=None)
df['bait']
for line in in_addr:
    line_l = line.split()

