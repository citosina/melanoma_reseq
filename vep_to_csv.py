



import pandas as pd
import io
import os


VEP = 'final.vep'



with open(VEP, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    vep = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    vep = vep.rename(columns={'#CHROM': 'CHROM'})
    vep = vep.set_index(['CHROM', 'POS'])




vep = vep.reset_index()
vep['ID_x'] = vep.CHROM + ':' + vep.POS +':'+vep.REF+':'+vep.ALT
vep = vep.set_index(['CHROM', 'POS'])
IDx = vep.ID_x.to_frame()
variant_info = vep.iloc[:,:6]
variant_info = variant_info.drop(axis = 1, labels = ['ID', 'QUAL', 'FILTER'])
variant_info['samtools'] = variant_info.INFO.str.split('|').str.get(0)
variant_info['VEP'] = variant_info.INFO.str.split(';').str.get(-1)
variant_info = variant_info.reset_index().set_index(['CHROM', 'POS'])
VEP = variant_info.VEP.str.split(',', expand = True)
VEP = VEP.dropna(axis = 1, how = 'all')

stack = VEP.stack()


stack = stack.to_frame()


stack['Allele']= stack[0].str.split('|').str.get(0)
stack['Consequence']= stack[0].str.split('|').str.get(1)
stack['IMPACT']= stack[0].str.split('|').str.get(2)
stack['SYMBOL']= stack[0].str.split('|').str.get(3)
stack['Gene']= stack[0].str.split('|').str.get(4)
stack['Feature_type']= stack[0].str.split('|').str.get(5)
stack['Feature']= stack[0].str.split('|').str.get(6)
stack['BIOTYPE']= stack[0].str.split('|').str.get(7)
stack['EXON']= stack[0].str.split('|').str.get(8)
stack['INTRON']= stack[0].str.split('|').str.get(9)
stack['HGVSc']= stack[0].str.split('|').str.get(10)
stack['HGVSp']= stack[0].str.split('|').str.get(11)
stack['cDNA_position']= stack[0].str.split('|').str.get(12)
stack['CDS_position']= stack[0].str.split('|').str.get(13)
stack['Protein_position']= stack[0].str.split('|').str.get(14)
stack['Amino_acids']= stack[0].str.split('|').str.get(15)
stack['Codons']= stack[0].str.split('|').str.get(16)
stack['Existing_variation']= stack[0].str.split('|').str.get(17)
stack['DISTANCE']= stack[0].str.split('|').str.get(18)
stack['STRAND']= stack[0].str.split('|').str.get(19)
stack['FLAGS']= stack[0].str.split('|').str.get(20)
stack['VARIANT_CLASS']= stack[0].str.split('|').str.get(21)
stack['SYMBOL_SOURCE']= stack[0].str.split('|').str.get(22)
stack['HGNC_ID']= stack[0].str.split('|').str.get(23)
stack['CANONICAL']= stack[0].str.split('|').str.get(24)
stack['MANE']= stack[0].str.split('|').str.get(25)
stack['TSL']= stack[0].str.split('|').str.get(26)
stack['APPRIS']= stack[0].str.split('|').str.get(27)
stack['CCDS']= stack[0].str.split('|').str.get(28)
stack['ENSP']= stack[0].str.split('|').str.get(29)
stack['SWISSPROT']= stack[0].str.split('|').str.get(30)
stack['TREMBL']= stack[0].str.split('|').str.get(31)
stack['UNIPARC']= stack[0].str.split('|').str.get(32)
stack['GENE_PHENO']= stack[0].str.split('|').str.get(33)
stack['SIFT']= stack[0].str.split('|').str.get(34)
stack['PolyPhen']= stack[0].str.split('|').str.get(35)
stack['DOMAINS']= stack[0].str.split('|').str.get(36)
stack['miRNA']= stack[0].str.split('|').str.get(37)
stack['HGVS_OFFSET']= stack[0].str.split('|').str.get(38)
stack['AF']= stack[0].str.split('|').str.get(39)
stack['AFR_AF']= stack[0].str.split('|').str.get(40)
stack['AMR_AF']= stack[0].str.split('|').str.get(41)
stack['EAS_AF']= stack[0].str.split('|').str.get(42)
stack['EUR_AF']= stack[0].str.split('|').str.get(43)
stack['SAS_AF']= stack[0].str.split('|').str.get(44)
stack['AA_AF']= stack[0].str.split('|').str.get(45)
stack['EA_AF']= stack[0].str.split('|').str.get(46)
stack['gnomAD_AF']= stack[0].str.split('|').str.get(47)
stack['gnomAD_AFR_AF']= stack[0].str.split('|').str.get(48)
stack['gnomAD_AMR_AF']= stack[0].str.split('|').str.get(49)
stack['gnomAD_ASJ_AF']= stack[0].str.split('|').str.get(50)
stack['gnomAD_EAS_AF']= stack[0].str.split('|').str.get(51)
stack['gnomAD_FIN_AF']= stack[0].str.split('|').str.get(52)
stack['gnomAD_NFE_AF']= stack[0].str.split('|').str.get(53)
stack['gnomAD_OTH_AF']= stack[0].str.split('|').str.get(54)
stack['gnomAD_SAS_AF']= stack[0].str.split('|').str.get(55)
stack['MAX_AF']= stack[0].str.split('|').str.get(56)
stack['MAX_AF_POPS']= stack[0].str.split('|').str.get(57)
stack['CLIN_SIG']= stack[0].str.split('|').str.get(58)
stack['SOMATIC']= stack[0].str.split('|').str.get(59)
stack['PHENO']= stack[0].str.split('|').str.get(60)
stack['PUBMED']= stack[0].str.split('|').str.get(61)
stack['MOTIF_NAME']= stack[0].str.split('|').str.get(62)
stack['MOTIF_POS']= stack[0].str.split('|').str.get(63)
stack['HIGH_INF_POS']= stack[0].str.split('|').str.get(64)
stack['MOTIF_SCORE_CHANGE']= stack[0].str.split('|').str.get(65)


stack = stack.drop(axis = 1, labels = 0)

stack = stack.reset_index()

stack['SNP'] = stack.CHROM + ':' +stack.POS.astype(str)

stack = stack.set_index('CHROM', 'POS')
stack.to_csv('passvep_stack_vep.csv')

stack = stack.reset_index()
stack['BP'] = stack['POS']

chr1 = stack[(stack['CHROM'] == 'chr1')]
chr2 = stack[(stack['CHROM'] == 'chr2')]
chr3 = stack[(stack['CHROM'] == 'chr3')]
chr4 = stack[(stack['CHROM'] == 'chr4')]
chr5 = stack[(stack['CHROM'] == 'chr5')]
chr6 = stack[(stack['CHROM'] == 'chr6')]
chr7 = stack[(stack['CHROM'] == 'chr7')]
chr9 = stack[(stack['CHROM'] == 'chr9')]
chr10 = stack[(stack['CHROM'] == 'chr10')]
chr11 = stack[(stack['CHROM'] == 'chr11')]
chr12 = stack[(stack['CHROM'] == 'chr12')]
chr14 = stack[(stack['CHROM'] == 'chr14')]
chr15 = stack[(stack['CHROM'] == 'chr15')]
chr16 = stack[(stack['CHROM'] == 'chr16')]
chr17 = stack[(stack['CHROM'] == 'chr17')]
chr18 = stack[(stack['CHROM'] == 'chr18')]
chr20 = stack[(stack['CHROM'] == 'chr20')]
chr21 = stack[(stack['CHROM'] == 'chr21')]
chr22 = stack[(stack['CHROM'] == 'chr22')]
ARNT = chr1[(chr1['BP'].astype(int) > 150754900) & (chr1['BP'].astype(int) < 151073108)].rename(columns={'POS':'chr1'})
PARP1 = chr1[(chr1['BP'].astype(int) > 226258651) & (chr1['BP'].astype(int) < 226499883)].rename(columns={'POS':'chr1'})
CASP8 = chr2[(chr2['BP'].astype(int) > 201047085) & (chr2['BP'].astype(int) < 201724509)].rename(columns={'POS':'chr2'})
SOD3 = chr4[(chr4['BP'].astype(int) > 24692289) & (chr4['BP'].astype(int) <  25055958)].rename(columns={'POS':'chr4'})
TERT = chr5[(chr5['BP'].astype(int) > 1142197) & (chr5['BP'].astype(int) < 1447788)].rename(columns={'POS':'chr5'})
SLC45A2 = chr5[(chr5['BP'].astype(int) > 33840139) & (chr5['BP'].astype(int) < 34064027)].rename(columns={'POS':'chr5'})
CYB5R4 = chr6[(chr6['BP'].astype(int) > 3948233) & (chr6['BP'].astype(int) <  84292890)].rename(columns={'POS':'chr6'})
MTAP = chr9[(chr9['BP'].astype(int) > 21665209) & (chr9['BP'].astype(int) <  22151817)].rename(columns={'POS':'chr9'})
CCND1 = chr11[(chr11['BP'].astype(int) > 69210414) & (chr11['BP'].astype(int) < 69803433)].rename(columns={'POS':'chr11'})
TYR = chr11[(chr11['BP'].astype(int) > 89097917) & (chr11['BP'].astype(int) <  89675699)].rename(columns={'POS':'chr11'})
ATM = chr11[(chr11['BP'].astype(int) > 108136364) & (chr11['BP'].astype(int) <  108468324)].rename(columns={'POS':'chr11'})
OCA2 = chr15[(chr15['BP'].astype(int) > 27670078) & (chr15['BP'].astype(int) < 28337214)].rename(columns={'POS':'chr15'})
FTO = chr16[(chr16['BP'].astype(int) > 53363904 ) & (chr16['BP'].astype(int) <  54384452)].rename(columns={'POS':'chr16'})
CDH1 = chr16[(chr16['BP'].astype(int) > 68528729) & (chr16['BP'].astype(int) < 68923897)].rename(columns={'POS':'chr16'})
MC1R = chr16[(chr16['BP'].astype(int) > 89806464) & (chr16['BP'].astype(int) < 90031424)].rename(columns={'POS':'chr16'})
CARD14 = chr17[(chr17['BP'].astype(int) > 80154433) & (chr17['BP'].astype(int) < 80508392)].rename(columns={'POS':'chr17'})
ASIP = chr20[(chr20['BP'].astype(int) > 33709976) & (chr20['BP'].astype(int) < 34374270)].rename(columns={'POS':'chr20'})
MX2 = chr21[(chr21['BP'].astype(int) > 41194545) & (chr21['BP'].astype(int) <  41512249)].rename(columns={'POS':'chr21'})
PLA2G6 = chr22[(chr22['BP'].astype(int) > 37870877) & (chr22['BP'].astype(int) <  38328731)].rename(columns={'POS':'chr22'})

ARNT.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/ARNT/ARNT_vep.csv')
PARP1.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/PARP1/PARP1_vep.csv') 
CASP8.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/CASP8/CASP8_vep.csv') 
SOD3.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/SOD3/SOD3_vep.csv') 
TERT.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/TERT/TERT_vep.csv') 
SLC45A2.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/SLC45A2/SLC45A2_vep.csv') 
CYB5R4.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/CYB5R4/CYB5R4_vep.csv') 
MTAP.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/MTAP/MTAP_vep.csv') 
CCND1.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/CCND1/CCND1_vep.csv') 
TYR.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/TYR/TYR_vep.csv') 
ATM.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/ATM/ATM_vep.csv') 
OCA2.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/OCA2/OCA2_vep.csv') 
FTO.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/FTO/FTO_vep.csv') 
CDH1.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/CDH1/CDH1_vep.csv') 
MC1R.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/MC1R/MC1R_vep.csv') 
CARD14.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/CARD14/CARD14_vep.csv') 
ASIP.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/ASIP/ASIP_vep.csv') 
MX2.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/MX2/MX2_vep.csv') 
PLA2G6.to_csv('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/PLA2G6/PLA2G6_vep.csv') 


