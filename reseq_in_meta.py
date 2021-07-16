

 import pandas as pd
 meta = pd.read_csv('../archivos_fuente/metafull3.txt', delimiter = '\s+', low_memory = False)
 meta = meta.set_index('SNP')
 snp = pd.read_csv('final_SNPS_37')
 snp = snp.set_index('gr37')
 a = snp.merge(meta, right_index = True, left_index = True)

 
 a.to_csv('reseq_in_meta.csv')
