




import os
import io
import pandas as pd
import numpy as np
import subprocess

region = str(input("enter the region you want to regress:"))
chromosome = str(input("enter the chromosome that region is at : "))

wd = '/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/analysis/'+region
os.chdir(wd)

fisher = pd.read_csv(wd+'/'+region+'_fisher.csv', low_memory=False).drop(axis = 1, labels = { 'ID_x', 'CHR'})
fisher = fisher.sort_values('P').set_index('SNP')
fisher['F_T'] = (fisher.F_A + fisher.F_U)/2
ref_alt = fisher.A1 + '_' + fisher.A2
ref_alt = ref_alt.to_frame().rename(columns = {0:'ref_alt'})
vep = pd.read_csv(region+'_vep.csv', low_memory = False).drop(axis =1, labels = {'Unnamed: 0', 'level_2'}).set_index('SNP')
rs = vep[['Feature', 'Existing_variation', 'Consequence', 'CANONICAL', 'Protein_position','Amino_acids']]
fisher = fisher.merge(rs, right_index = True, left_index = True).sort_values('P').drop_duplicates()
fisher = fisher.reset_index().sort_values('P')
fisher['ID_x'] = fisher.SNP
top = fisher.iloc[0,:]
top_fisher = [top.ID_x]
with open('to_condition', 'w') as f:
    for var in top_fisher:
        f.write("%s\n" % var)

os.chmod('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_1.sh', 0o755)
with open('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_1.sh', 'rb') as file:
    script = file.read()
rc = subprocess.call(script, shell=True)
fisher = fisher.drop(axis = 1, labels = 'BP')

cond_1 = pd.read_csv(wd + '/plink/1ADD.csv', delimiter = ',', header = None ).drop\
(axis = 1, labels = 0).rename(columns = {1:'CHR', 2:'SNP',3:'BP', 4:'A1',5:'TEST', 6:'NMISS', 7:'OR', 8:'SATT', 9:'P'}).sort_values\
('P').set_index('SNP').fillna(1)
cond_1['-logP'] = -np.log10(cond_1.P)
cond_1 = cond_1.merge(ref_alt, right_index = True, left_index = True)
fisher = fisher.set_index('SNP')
cond__1 = cond_1.reset_index().sort_values('P')

cond_1 = cond_1.rename(columns = {'OR':'OR_cond1', 'P':'P_cond1', '-logP':'-logP_cond1'})
fisher=fisher.rename(columns = {'P': 'P_fisher', '-logP':'-logP_fisher', 'OR':'OR_fisher'})
print(cond_1) 
r_1 = cond_1.merge(fisher, right_index = True, left_index = True)

r1 = r_1[(r_1['-logP_fisher'].astype(float) > 3) & (r_1['-logP_cond1'].astype(float) < 1)]

r1['GRCh38'] = r1.ID_x + '-' + r1.BP.astype(str)
r1['Signal']  = 's1'
r1['Region'] = region

r1=r1[['Signal', 'Region', 'Existing_variation', 'GRCh38', 'ref_alt','-logP_fisher', 'Consequence', 'Protein_position', 'Amino_acids','Feature','CANONICAL',  'OR_fisher', 'OR_cond1', 'P_fisher', 'P_cond1']]
r1['CHROM'] = r1.GRCh38.str.split(':').str.get(0)
r1['start'] = r1.GRCh38.str.split(':').str.get(1).str.split('-').str.get(0)
r1['end'] = r1.GRCh38.str.split(':').str.get(1).str.split('-').str.get(1)
r1BED = r1[['CHROM','start','end',  'Existing_variation']]
r1 = r1.drop(axis =1, labels = ['CHROM', 'start', 'end'])
print(r1)
r1.to_csv('r1', index = None)

emptyr1 = r1.empty
if emptyr1:
        print('no significant variants after the first conditioning')

elif emptyr1==False:
        top_cond1 = [ top.SNP, cond__1.iloc[0,:].SNP]



with open('to_condition', 'w') as f:
    for var in top_cond1:
        f.write("%s\n" % var)

os.chmod('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_2.sh', 0o755)

with open('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_2.sh', 'rb') as file:
    script = file.read()
rc = subprocess.call(script, shell=True)

cond_2 = pd.read_csv(wd + '/plink/2ADD.csv', delimiter = ',', header = None ).drop\
(axis = 1, labels = 0).rename(columns = {1:'CHR', 2:'SNP',3:'BP', 4:'A1',5:'TEST', 6:'NMISS', 7:'OR', 8:'SATT', 9:'P'}).sort_values('P')
cond_2['-logP'] = -np.log10(cond_2.P)
cond_2 = cond_2.set_index('SNP')
cond_2 = cond_2.merge(ref_alt, right_index = True, left_index = True)
cond__2 = cond_2.reset_index().sort_values('P')
cond_2 = cond_2.rename(columns = {'OR':'OR_cond2', 'P':'P_cond2', '-logP':'-logP_cond2'})

r_2 = cond_2.merge(r_1, right_index = True, left_index = True)

r_2['ID_x'] = 'chr' + r_2.CHR_x.astype(str) + ':' + r_2.BP_x.astype(str)

r2 = r_2[(r_2['-logP_cond1'].astype(float) > 3) & (r_2['-logP_cond2'].astype(float) < 1)]

r2['GRCh38'] = r2.ID_x + '-' + r2.ID_x.str.split(':').str.get(1)
r2BED = r2[['CHR_x','BP_x','BP_x',  'Existing_variation']]

r2['Signal']  = 's2'
r2['Region'] = region

r2=r2[['Signal', 'Region', 'Existing_variation', 'GRCh38', 'ref_alt_x','-logP_fisher', 'Consequence', 'Protein_position', 'Amino_acids', 'Feature', 'CANONICAL', 'OR_fisher', 'OR_cond1', 'OR_cond2', 'P_fisher', 'P_cond1', 'P_cond2']]
r2.to_csv('r2', index = False)
emptyr2 = r2.empty

if emptyr2:
        print("There are no significant variants after the second logistic regression")
        print("Finished, candidates saved to r files, position bed file saved to candidates.bed")

        BED = (r1BED)
        BED.to_csv('candidates.bed', sep = '\t', index = None)

elif emptyr2==False:
        top_cond2 = [ top.ID_x, cond__1.iloc[0,:].SNP, cond__2.iloc[0,:].SNP]
        print(cond__2.head())

        with open('to_condition', 'w') as f:
                for var in top_cond2:
                        f.write("%s\n" % var)
        os.chmod('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_3.sh', 0o755)

        with open('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/scripts/logreg_3.sh', 'rb') as file:
                script = file.read()
        rc = subprocess.call(script, shell=True)
        cond_3 = pd.read_csv(wd + '/plink/3ADD.csv', delimiter = ',', header = None ).drop\
        (axis = 1, labels = 0).rename(columns = {1:'CHR', 2:'SNP',3:chromosome, 4:'A1',5:'TEST', 6:'NMISS', 7:'OR', 8:'SATT', 9:'P'}).sort_values('P')
        cond_3['-logP'] = -np.log10(cond_3.P)
        cond_3 = cond_3.set_index('SNP')
        cond__3 = cond_3.reset_index().sort_values('P')
        cond_3 = cond_3.rename(columns = {'OR':'OR_cond3', 'P':'P_cond3', '-logP':'-logP_cond3'})
        r3 = cond_3.merge(r_2, right_index = True, left_index = True)
        r3['ID_x'] = 'chr' + r3.CHR_x.astype(str) + ':' + r3.BP_x.astype(str)
        r3 = r3[(r3['-logP_cond2'].astype(float) > 3) & (r3['-logP_cond3'].astype(float) < 1)]
        r3['GRCh38'] = r3.ID_x + '-' + r3.ID_x.str.split(':').str.get(1)
        emptyr3 = r3.empty
        if emptyr3:
                print("There are no significant variants after the third logistic regression")
                BED = (r1BED, r2BED)
                BED1 = r1BED.join(r2BED, lsuffix='r1',rsuffix ='r2')
                BED1.to_csv('candidates.bed', sep = '\t', index = None)
        elif emptyr3==False:
                r3=r3[['ref_alt_x', '-logP_cond2','-logP_cond3', 'GRCh38', 'Existing_variation']]
                r3 = r3.merge(rs, right_index = True, left_index = True).drop_duplicates()
                r3.to_csv('r3', header = None)
                r3['CHROM'] = r3.GRCh38.str.split(':').str.get(0)
                r3['start'] = r3.GRCh38.str.split(':').str.get(1).str.split('-').str.get(0)
                r3['end'] = r3.GRCh38.str.split(':').str.get(1).str.split('-').str.get(1)
                r3BED = r3[['CHROM','start','end',  'Existing_variation']]
                print(r3)
                BED = BED1.join(r2BED, lsuffix='r1',rsuffix ='r2')
                BED.to_csv('candidates.bed', sep = '\t', index = None)
                print("Finished, candidates saved to r files, position bed file saved to candidates.bed")

  
'''
os.chmod('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/lcr_filtered/scripts/gwas-reseq2/scripts/bash_scripts/vep_annotate_candidates.sh', 0o755)

with open('/mnt/Adenina/drobles/ccastaneda/gwas-reseq2/lcr_filtered/scripts/gwas-reseq2/scripts/bash_scripts/vep_annotate_candidates.sh', 'rb') as file:
	script = file.read()
        
	rc = subprocess.call(script, shell=True)
'''

print('annotations done')
exit()


