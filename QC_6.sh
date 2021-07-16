





plink --vcf filtered.recode.vcf --assoc fisher --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq2/archivos_fuente/pheno --out filteres --allow-no-sex --keep-allele-order --set-missing-var-ids @:#                                                         


cat filteres.assoc.fisher | tr -s ' ' ',' > filteres.fisher.csv


source activate vep

/mnt/Adenina/drobles/ccastaneda/home/ensembl-vep/vep -i filteres.recode.vcf  --everything --cache -o filteres.vep  --output_format vcf --no_stats --force_overwrite

python vep_to_csv.py
          
