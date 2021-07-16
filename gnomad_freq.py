import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import math 
import numpy as np
import subprocess

info = pd.read_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/QC/final.fisher.csv')
gnomad = pd.read_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/QC/leeds_resources/gnomad_leeds_regions_nfe_AF.table', delimiter = '\t')


info['-logP'] = -(np.log10(info.P))
info['SNP_reseq'] = info.SNP +':'+ info.A2 + info.A1
info = info.drop(axis = 1, labels = 'SNP')
info = info.set_index('SNP_reseq')
info = info.drop(axis = 1, labels = {'Unnamed: 0','Unnamed: 10'})
print(info.sort_values('P').head())
print(info.info())


gnomad['SNP_reseq'] = gnomad['CHROM'] +':'+gnomad['POS'].astype(str)+':'+gnomad.REF+gnomad.ALT
gnomad['gnomad_AF'] = gnomad.AF_nfe
gnomad = gnomad.set_index('SNP_reseq')
everything = info.merge(gnomad, right_index = True, left_index = True )
print(everything)
print(everything.info())

h = everything

sns.scatterplot(x="F_U", y='gnomad_AF', data=h, legend = False)
plt.savefig('F_Uvsgnomad_AF.tif')

h.to_csv('gnomad_final.csv')

