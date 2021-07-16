import pandas as pd
import io
import os
import numpy as np


ctrl_b1 = '../QC/ctrl_b1_hwe.recode.vcf'
with open(ctrl_b1, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    ctrl_b1 = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    ctrl_b1 = ctrl_b1.rename(columns={'#CHROM': 'CHROM'})
    ctrl_b1 = ctrl_b1.set_index(['CHROM', 'POS'])
    ctrl_b1_GQ = ctrl_b1.iloc[:,7:].applymap(lambda s: str(s.split(':')[3]))
    ctrl_b1_gq = ctrl_b1_GQ
    ctrl_b1_gq = ctrl_b1_gq.replace('.', np.NaN).applymap(lambda s: float(s))
    ctrl_b1_gq['avg'] = ctrl_b1_gq[ctrl_b1_gq.columns].mean(axis=1)
    
    
ctrl_b2 = '../QC/ctrl_b2_hwe.recode.vcf'
with open(ctrl_b2, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    ctrl_b2 = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    ctrl_b2 = ctrl_b2.rename(columns={'#CHROM': 'CHROM'})
    ctrl_b2 = ctrl_b2.set_index(['CHROM', 'POS'])
    ctrl_b2_GQ = ctrl_b2.iloc[:,7:].applymap(lambda s: str(s.split(':')[3]))
    ctrl_b2_gq = ctrl_b2_GQ
    ctrl_b2_gq = ctrl_b2_gq.replace('.', np.NaN).applymap(lambda s: float(s))
    ctrl_b2_gq['avg'] = ctrl_b2_gq[ctrl_b2_gq.columns].mean(axis=1)
    
ctrl_b1_gq_avg = ctrl_b1_gq['avg']
ctrl_b2_gq_avg = ctrl_b2_gq['avg']
ctrl_b1_gq_avg.to_csv('../QC/ctrl_b1_GQavg.report')
ctrl_b2_gq_avg.to_csv('../QC/ctrl_b2__GQavg.report')

ctrl_b1_gq_avgpos = ctrl_b1_gq_avg[ctrl_b1_gq_avg < 30].reset_index().drop(axis = 1, labels = 'avg').to_csv('../QC/ctrl_b1_gq_avgPASS.bed', sep = '\t', index = False)
ctrl_b2_gq_avgpos = ctrl_b2_gq_avg[ctrl_b2_gq_avg < 30].reset_index().drop(axis = 1, labels = 'avg').to_csv('../QC/ctrl_b2_gq_avgPASS.bed', sep = '\t', index = False)

case_b1 = '../archivos_fuente/case_b1.recode.vcf'
with open(case_b1, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    case_b1 = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    case_b1 = case_b1.rename(columns={'#CHROM': 'CHROM'})
    case_b1 = case_b1.set_index(['CHROM', 'POS'])
    case_b1_GQ = case_b1.iloc[:,7:].applymap(lambda s: str(s.split(':')[3]))
    case_b1_gq = case_b1_GQ
    case_b1_gq = case_b1_gq.replace('.', np.NaN).applymap(lambda s: float(s))
    case_b1_gq['avg'] = case_b1_gq[case_b1_gq.columns].mean(axis=1)
    
    
case_b2 = '../archivos_fuente/case_b2.recode.vcf'
with open(case_b2, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    case_b2 = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    case_b2 = case_b2.rename(columns={'#CHROM': 'CHROM'})
    case_b2 = case_b2.set_index(['CHROM', 'POS'])
    case_b2_GQ = case_b2.iloc[:,7:].applymap(lambda s: str(s.split(':')[3]))
    case_b2_gq = case_b2_GQ
    case_b2_gq = case_b2_gq.replace('.', np.NaN).applymap(lambda s: float(s))
    case_b2_gq['avg'] = case_b2_gq[case_b2_gq.columns].mean(axis=1)
    
case_b1_gq_avg = case_b1_gq['avg']
case_b2_gq_avg = case_b2_gq['avg']
case_b1_gq_avg.to_csv('../QC/case_b1_GQavg.report')
case_b2_gq_avg.to_csv('../QC/case_b2_GQavg.report')

case_b1_gq_avgpos = case_b1_gq_avg[case_b1_gq_avg < 30].reset_index().drop(axis = 1, labels = 'avg').to_csv('../QC/case_b1_gq_avgPASS.bed', sep = '\t', index = False)
case_b2_gq_avgpos = case_b2_gq_avg[case_b2_gq_avg < 30].reset_index().drop(axis = 1, labels = 'avg').to_csv('../QC/case_b2_gq_avgPASS.bed', sep = '\t', index = False)



