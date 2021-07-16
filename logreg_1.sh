

plink --bfile ./plink/plink --allow-no-sex --condition-list to_condition --logistic --out ./plink/first_conditional --pheno /mnt/Adenina/drobles/ccastaneda/gwas-reseq2/archivos_fuente/pheno --keep-allele-order

cat ./plink/first_conditional.assoc.logistic  | tr -s ' ' ',' | grep 'ADD' > ./plink/1ADD.csv

 
