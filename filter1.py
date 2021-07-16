import pandas as pd
import numpy as np
import os
import io 



ctrl_b1 = pd.read_csv('ctrl_b1_POS.csv', delimiter = ' ', header = None)
ctrl_b2 = pd.read_csv('ctrl_b2_POS.csv', delimiter = ' ', header = None)
case_b1 = pd.read_csv('case_b1_POS.csv', delimiter = ' ', header = None)
case_b2 = pd.read_csv('case_b2_POS.csv', delimiter = ' ', header = None)

cohort = pd.read_csv('cohort_POS.csv', delimiter = ' ', header = None)

ctrl_b1['SNP'] = ctrl_b1[0]+ ':' + ctrl_b1[2].astype(str)
ctrl_b2['SNP'] = ctrl_b2[0]+ ':' + ctrl_b2[2].astype(str)
case_b1['SNP'] = case_b1[0]+ ':' + case_b1[2].astype(str)
case_b2['SNP'] = case_b2[0]+ ':' + case_b2[2].astype(str)
cohort['SNP'] = cohort[0]+ ':' + cohort[2].astype(str)

ctrl_b1 = ctrl_b1.set_index ('SNP')
ctrl_b2 = ctrl_b2.set_index('SNP')
case_b1 = case_b1.set_index('SNP')
case_b2 = case_b2.set_index('SNP')
cohort = cohort.set_index('SNP')

cohort_ctrl_b1 = cohort.merge(ctrl_b1, right_index = True, left_index = True, indicator = True)
cohort_ctrl_b2 = cohort_ctrl_b1.merge(ctrl_b2, right_index = True, left_index = True)
cohort_case_b1 = cohort_ctrl_b2.merge(case_b1, right_index = True, left_index = True)
cohort_case_b2 = cohort_case_b1.merge(cohort_case_b1, right_index = True, left_index = True)

pass_all_filters = cohort_case_b2['0_x_x']

pass_all_filters = pass_all_filters.reset_index().drop(axis = 1, labels = '0_x_x')

pass_all_filters['CHROM'] = pass_all_filters.SNP.str.split(':').str.get(0)
pass_all_filters['POS'] = pass_all_filters.SNP.str.split(':').str.get(1)

pass_all_filters = pass_all_filters.drop(axis = 1, labels = 'SNP')
pass_all_filters.to_csv('pass_filters1', index = None, header = None, sep = '\t')

