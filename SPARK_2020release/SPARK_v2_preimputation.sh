###Scripts for new spark dataset

setwd("/mnt/b2/home4/arc/vw260/SPARK/SPARK_v2/")

ancestry = fread("SPARK.WES2.ancestry.2021_03.tsv")

head(ancestry)

EUR = subset(ancestry, population == "EUR")
fam = fread("SPARK.WES2.release.2020_06.genotype.fam")
fam2 = fam[fam$V2 %in% EUR$spid,] ##11737

write.table(fam2[,c("V1", "V2")], file = "EUR_keep.txt", row.names = F, col.names = F, quote = F)

./plink --bfile SPARK.WES2.release.2020_06.genotype --keep EUR_keep.txt --make-bed --geno 0.1 --out QC1output --threads 15
./plink --bfile QC1output --make-bed --me 0.05 0.1 --mind 0.05 --out QC2output --threads 15
./plink --bfile QC2output --make-bed --out QC3output --remove QC2output.irem --threads 15
./plink --bfile QC3output --hwe 0.000001 --make-bed --geno 0.05 --threads 15 --out QC4output 
./plink --bfile QC4output --het --out check_het --threads 15

R
library(data.table)
het = fread("check_het.het")
het$HET = (het$`N(NM)` - het$`O(HOM)`)/het$`N(NM)` #create heterozygosity stats
mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity
hetoutlier = subset(het, abs(Z) > 3)
failedsample = hetoutlier[,c(1:2)]
write.table(failedsample, file = "failedsample_het.txt", row.names = F, col.names = T, quote = F)
q()
n

./plink --bfile QC4output --make-bed --out QC5output --remove failedsample_het.txt --threads 15

R
library(data.table)
map = fread("QC5output.bim")
map$update <- gsub("GSA-", "", map$V2)
write.table(map[,c("V2", "update")], file = "updatenames.txt", row.names = F, col.names = F, quote = F)
q()
n

./plink --bfile QC5output --update-name updatenames.txt --make-bed --out QC6output 
./plink --bfile QC6output --freq --out QC6outputfreq

perl HRC-1000G-check-bim.pl -b QC6output.bim -f QC6outputfreq.frq -r 1000GP_Phase3_combined.legend -g

./plink --bfile QC6output --exclude Exclude-QC6output-1000G.txt --make-bed --out TEMP1 --threads 20
./plink --bfile TEMP1 --update-map Chromosome-QC6output-1000G.txt --update-chr --make-bed --out TEMP2 --threads 20
./plink --bfile TEMP2 --update-map Position-QC6output-1000G.txt --make-bed --out TEMP3 --threads 20
./plink --bfile TEMP3 --flip Strand-Flip-QC6output-1000G.txt --make-bed --out TEMP4 --threads 20
./plink --bfile TEMP4 --reference-allele Force-Allele1-QC6output-1000G.txt --make-bed --out QC6output-updated --threads 20
rm TEMP*
for i in {1..23}; do ./plink --bfile QC6output-updated --reference-allele Force-Allele1-QC6output-1000G.txt --chr ${i} --recode-vcf --out QC6output_file_chr${i}; done
for i in {1..23}; do vcf-sort QC6output_file_chr${i}.vcf | bgzip -c > QC6output_file_chr${i}.vcf.gz; done
for i in {1..23}; do rm QC6output_file_chr${i}.vcf | rm QC6output_file_chr${i}.log; done