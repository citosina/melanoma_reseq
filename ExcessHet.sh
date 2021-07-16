source ~/anaconda3/etc/profile.d/conda.sh
module load gatk/4.1.1.0 
module load plink2

conda activate gatk2

gatk --java-options "-Xmx3g -Xms3g" VariantFiltration --variant out.recode.vcf --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O cohort_excesshet.vcf.gz
