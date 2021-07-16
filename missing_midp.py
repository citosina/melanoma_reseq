import pandas as pd
import os
import numpy as np
import io 


missing = pd.read_csv('cohort.missing', delimiter = "\s+")
missing['CHROM'] = missing.SNP.str.split(':').str.get(0)
missing['POS'] = missing.SNP.str.split(':').str.get(1)
missing['POS'] = missing.POS.str.split(':').str.get(0)


no_pass = missing[missing.P <0.05]

no_pass = no_pass.iloc[:, [-2,-1]]

no_pass.to_csv('missing_no_pass', index = None, header = None, sep = '\t')
