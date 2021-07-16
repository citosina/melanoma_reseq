


import pandas as pd
import numpy as np
a = pd.read_csv( 'final.fisher.csv', low_memory = False, index_col = None)
a['CHR'] = a.CHR.apply(lambda s: int(s))
a['ID_x'] = a.SNP 
a['POS'] = a.BP
a['CHROM'] = 'chr' + a.CHR.astype(str)
a = a.set_index('SNP')
a['-logP'] = -np.log10(a.P)



chr1 = a[(a['CHROM'] == 'chr1')]
chr2 = a[(a['CHROM'] == 'chr2')]
chr3 = a[(a['CHROM'] == 'chr3')]
chr4 = a[(a['CHROM'] == 'chr4')]
chr5 = a[(a['CHROM'] == 'chr5')]
chr6 = a[(a['CHROM'] == 'chr6')]
chr7 = a[(a['CHROM'] == 'chr7')]
chr9 = a[(a['CHROM'] == 'chr9')]
chr10 = a[(a['CHROM'] == 'chr10')]
chr11 = a[(a['CHROM'] == 'chr11')]
chr12 = a[(a['CHROM'] == 'chr12')]
chr14 = a[(a['CHROM'] == 'chr14')]
chr15 = a[(a['CHROM'] == 'chr15')]
chr16 = a[(a['CHROM'] == 'chr16')]
chr17 = a[(a['CHROM'] == 'chr17')]
chr18 = a[(a['CHROM'] == 'chr18')]
chr20 = a[(a['CHROM'] == 'chr20')]
chr21 = a[(a['CHROM'] == 'chr21')]
chr22 = a[(a['CHROM'] == 'chr22')]
ARNT = chr1[(chr1['BP'].astype(int) > 150754900) & (chr1['BP'].astype(int) < 151073108)].sort_values('P').rename(columns={'POS':'chr1'})
PARP1 = chr1[(chr1['BP'].astype(int) > 226258651) & (chr1['BP'].astype(int) < 226499883)].sort_values('P').rename(columns={'POS':'chr1'})
CASP8 = chr2[(chr2['BP'].astype(int) > 201047085) & (chr2['BP'].astype(int) < 201724509)].sort_values('P').rename(columns={'POS':'chr2'})
SOD3 = chr4[(chr4['BP'].astype(int) > 24692289) & (chr4['BP'].astype(int) <  25055958)].sort_values('P').rename(columns={'POS':'chr4'})
TERT = chr5[(chr5['BP'].astype(int) > 1142197) & (chr5['BP'].astype(int) < 1447788)].sort_values('P').rename(columns={'POS':'chr5'})
SLC45A2 = chr5[(chr5['BP'].astype(int) > 33840139) & (chr5['BP'].astype(int) < 34064027)].sort_values('P').rename(columns={'POS':'chr5'})
IRF4 = chr6[(chr6['BP'].astype(int) > 281588) & (chr6['BP'].astype(int) < 512178)].sort_values('P').rename(columns={'POS':'chr6'})
CYB5R4 = chr6[(chr6['BP'].astype(int) > 3948233) & (chr6['BP'].astype(int) <  84292890)].sort_values('P').rename(columns={'POS':'chr6'})
MTAP = chr9[(chr9['BP'].astype(int) > 21665209) & (chr9['BP'].astype(int) <  22151817)].sort_values('P').rename(columns={'POS':'chr9'})
CCND1 = chr11[(chr11['BP'].astype(int) > 69210414) & (chr11['BP'].astype(int) < 69803433)].sort_values('P').rename(columns={'POS':'chr11'})
TYR = chr11[(chr11['BP'].astype(int) > 89097917) & (chr11['BP'].astype(int) <  89675699)].sort_values('P').rename(columns={'POS':'chr11'})
ATM = chr11[(chr11['BP'].astype(int) > 108136364) & (chr11['BP'].astype(int) <  108468324)].sort_values('P').rename(columns={'POS':'chr11'})
OCA2 = chr15[(chr15['BP'].astype(int) > 27670078) & (chr15['BP'].astype(int) < 28337214)].sort_values('P').rename(columns={'POS':'chr15'})
FTO = chr16[(chr16['BP'].astype(int) > 53363904 ) & (chr16['BP'].astype(int) <  54384452)].sort_values('P').rename(columns={'POS':'chr16'})
CDH1 = chr16[(chr16['BP'].astype(int) > 68528729) & (chr16['BP'].astype(int) < 68923897)].sort_values('P').rename(columns={'POS':'chr16'})
MC1R = chr16[(chr16['BP'].astype(int) > 89806464) & (chr16['BP'].astype(int) < 90031424)].sort_values('P').rename(columns={'POS':'chr16'})
CARD14 = chr17[(chr17['BP'].astype(int) > 80154433) & (chr17['BP'].astype(int) < 80508392)].sort_values('P').rename(columns={'POS':'chr17'})
ASIP = chr20[(chr20['BP'].astype(int) > 33709976) & (chr20['BP'].astype(int) < 34374270)].sort_values('P').rename(columns={'POS':'chr20'})
MX2 = chr21[(chr21['BP'].astype(int) > 41194545) & (chr21['BP'].astype(int) <  41512249)].sort_values('P').rename(columns={'POS':'chr21'})
PLA2G6 = chr22[(chr22['BP'].astype(int) > 37870877) & (chr22['BP'].astype(int) <  38328731)].sort_values('P').rename(columns={'POS':'chr22'})


ARNT.sort_values('P').to_csv('../analysis/ARNT/ARNT_fisher.csv')
PARP1.sort_values('P').to_csv('../analysis/PARP1/PARP1_fisher.csv') 
CASP8.sort_values('P').to_csv('../analysis/CASP8/CASP8_fisher.csv') 
SOD3.sort_values('P').to_csv('../analysis/SOD3/SOD3_fisher.csv') 
TERT.sort_values('P').to_csv('../analysis/TERT/TERT_fisher.csv') 
SLC45A2.sort_values('P').to_csv('../analysis/SLC45A2/SLC45A2_fisher.csv') 
CYB5R4.sort_values('P').to_csv('../analysis/CYB5R4/CYB5R4_fisher.csv') 
MTAP.sort_values('P').to_csv('../analysis/MTAP/MTAP_fisher.csv') 
CCND1.sort_values('P').to_csv('../analysis/CCND1/CCND1_fisher.csv') 
TYR.sort_values('P').to_csv('../analysis/TYR/TYR_fisher.csv') 
ATM.sort_values('P').to_csv('../analysis/ATM/ATM_fisher.csv') 
OCA2.sort_values('P').to_csv('../analysis/OCA2/OCA2_fisher.csv') 
FTO.sort_values('P').to_csv('../analysis/FTO/FTO_fisher.csv') 
CDH1.sort_values('P').to_csv('../analysis/CDH1/CDH1_fisher.csv') 
MC1R.sort_values('P').to_csv('../analysis/MC1R/MC1R_fisher.csv') 
CARD14.sort_values('P').to_csv('../analysis/CARD14/CARD14_fisher.csv') 
ASIP.sort_values('P').to_csv('../analysis/ASIP/ASIP_fisher.csv') 
MX2.sort_values('P').to_csv('../analysis/MX2/MX2_fisher.csv') 
PLA2G6.sort_values('P').to_csv('../analysis/PLA2G6/PLA2G6_fisher.csv') 



