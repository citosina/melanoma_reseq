



import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import math
import numpy as np
import subprocess



h = pd.read_csv('h.csv', low_memory = False)

h['kg_reseq_com'] = h.gnomad_AF - h.F_U

h['col_zscore'] = (h['kg_reseq_com'] - h['kg_reseq_com'].mean())/h['kg_reseq_com'].std(ddof=0)
h.to_csv('h_test.csv')
sam = h[(h.col_zscore >= -9) & (h.col_zscore <= 9)]

sns.scatterplot(x="F_U", y='gnomad_AF', data=sam, legend = False)
plt.xlabel('AF Controls')
# Set y-axis label
plt.ylabel('AF GnomAD: Non-Finnish European')




plt.savefig('F_Uvsgnomad_AF_corrected.tif')


pasan = h[(h.col_zscore >= -9) & (h.col_zscore <= 9)]
pasan = pasan.reset_index()

print('Numero de variantes antes del filtro: ' + str(h.reset_index().SNP_reseq.unique().size))
print('Numero de variantes despues del filtro: '+ str(sam.reset_index().SNP_reseq.unique().size))
print('# de variantes que no pasan: ' + str(h.reset_index().SNP_reseq.unique().size-pasan.SNP_reseq.unique().size))


pasan = pasan[['CHROM', 'POS']].copy()

pasan.to_csv('pasan.bed', index = None, sep = '\t')
no_pasan = h[(h.col_zscore <= -9) & (h.col_zscore >= 9)]
print(no_pasan)
no_pasan.to_csv('gnomad_AF_comparison_no_pass.csv')

compare = h[['CHROM', 'POS', 'F_U', 'gnomad_AF', 'col_zscore']].copy() 

compare.to_csv('af_compare.csv')
