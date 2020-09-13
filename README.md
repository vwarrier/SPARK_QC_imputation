# SPARK_QC_imputation
QC and imputation of the SPARK dataset


## Additional step: combine into 1 file for polygenic scoring
```{bash}
./plink --bfile SPARK_CEU_chr1_binary_hg19 --merge-list SPARKCEU_merge.txt --maf 0.001 --make-bed --out SPARKCEU_hg19_allchrs --threads 10
for i in {1..22}; do ./plink --bfile SPARK_CEU_chr${i}_binary_hg19 --exclude SPARKCEU_hg19_allchrs-merge.missnp --maf 0.001 --make-bed --out SPARKCEU_chr${i}_hg19_v2 --threads 10; done
./plink --bfile SPARKCEU_chr1_hg19_v2 --merge-list SPARKCEU_merge2.txt --maf 0.001 --make-bed --out SPARKCEU_hg19_allchrs --threads 10

for i in {1..22}; do rm SPARKCEU_chr${i}_hg19_v2*; done

```
### Generate Principal components
First, using the CEU files, understand relatedness

```{bash}
./king -b QC4_CEU.bed --kinship --prefix QC4_CEU_kinship
```

Next, use PCAir to create PCs

