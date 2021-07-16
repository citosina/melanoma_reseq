import pandas as pd 
import os


cwd = '/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/QC'
os.chdir(cwd)

a = pd.read_csv('mdust_overlap_val', delimiter = ' ', header = None)
a = a.rename(columns = {0:'CHROM', 1:'POS', 2:'REF', 3:'ALT', 4:'DUSTED'}).drop(axis = 1, labels = [5,6,7])
d = a[a.DUSTED != 0]
d.to_csv('mdust_vars.csv', index= None)
