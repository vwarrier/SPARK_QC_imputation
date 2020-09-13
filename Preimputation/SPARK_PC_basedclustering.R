library(umap)
library(ggplot2)
library(data.table)

## 6.2 Read additional population files and define races
load(PC.RData)

populationfile = fread("SPARK_1KG_population.txt")
setnames(PC, "Sample_name", "Sample")

populationfile = populationfile[!duplicated(populationfile$Sample),]

PC_withancestry = merge(populationfile, PC, by = "Sample")
PC_withancestry = unique(PC_withancestry)

#####We will fist identify clusters using only the 1KG dataset and then project the remaining data on it######

##6.3 UMAP with 1KG only (5 PCs, 8 PCs, and then 10 PCs)
# First subset to only 1 KG
PC_withancestry_1KG = subset(PC_withancestry, Dataset == "1000G")

# 6.3.1 Keep only the first five PCs and make plots
PC_forumap_1KG = PC_withancestry_1KG[,c(5:9)]
PC_forumap_labels_1KG = PC_withancestry_1KG[,c(1:4)]
PC_umap_1KG = umap(PC_forumap_1KG)

PC_umap_layout_1KG = PC_umap_1KG$layout
PC_umap_layout_1KG = cbind(PC_umap_layout_1KG, PC_forumap_labels_1KG)
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12))
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Population, colour = Population)) + scale_shape_manual(values=c(0,1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))

# 6.3.2 Keep only the first eight PCs and make plots
PC_forumap_1KG = PC_withancestry_1KG[,c(5:12)]
PC_forumap_labels_1KG = PC_withancestry_1KG[,c(1:4)]
PC_umap_1KG = umap(PC_forumap_1KG)

PC_umap_layout_1KG = PC_umap_1KG$layout
PC_umap_layout_1KG = cbind(PC_umap_layout_1KG, PC_forumap_labels_1KG)
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(0,1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Population, colour = Population)) + scale_shape_manual(values=c(0,1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))


# 6.3.3 Keep only the first ten PCs and make plots
PC_forumap_1KG = PC_withancestry_1KG[,c(5:14)]
PC_forumap_labels_1KG = PC_withancestry_1KG[,c(1:4)]
PC_umap_1KG = umap(PC_forumap_1KG)

PC_umap_layout_1KG = PC_umap_1KG$layout
PC_umap_layout_1KG = cbind(PC_umap_layout_1KG, PC_forumap_labels_1KG)
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(0,1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Population, colour = Population)) + scale_shape_manual(values=c(0,1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25))


###6.4 Go back to 5 
## UMAP with 1KG - 5 PCs

PC_forumap_1KG = PC_withancestry_1KG[,c(5:9)]
PC_forumap_labels_1KG = PC_withancestry_1KG[,c(1:4)]
PC_umap_1KG = umap(PC_forumap_1KG)


PC_umap_layout_1KG = PC_umap_1KG$layout
PC_umap_layout_1KG = cbind(PC_umap_layout_1KG, PC_forumap_labels_1KG)
ggplot(PC_umap_layout_1KG,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12))


##6.5 UMAP with others and project - 5 PCs
PC_withancestry_no1KG = subset(PC_withancestry, Dataset == "SPARK") ####Change ABCD to your population

PC_forumap_no1KG = PC_withancestry_no1KG[,c(5:9)]
PC_forumap_labels_no1KG = PC_withancestry_no1KG[,c(1:4)]
PC_umap_no1KG = predict(PC_umap_1KG, PC_forumap_no1KG)
PC_umap_layout_no1KG = cbind(PC_umap_no1KG, PC_forumap_labels_no1KG)


PC_all_umap = rbind(PC_umap_layout_1KG, PC_umap_layout_no1KG)

q = ggplot(PC_all_umap,aes(x=V1,y=V2,color=as.character(Superpopulation))) + geom_point(aes(shape = as.character(Superpopulation), colour = as.character(Population))) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))

q + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 1)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 1))



###6.6 Subset populations
American_chunk = subset(PC_all_umap, V1 > -1 & V1 < 8.5 & V2 < -6.5 )

p = ggplot(American_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))

TSI_chunk = subset(PC_all_umap, V1 < -4 & V1 > -6.5 & V2 < -5 & V2 > -7)
p = ggplot(TSI_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))


Finnish_chunk = subset(PC_all_umap, V1 < -12.5 & V2 < -6 & V2 > -9)
p = ggplot(Finnish_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))

CEU_chunk = subset(PC_all_umap, V1 > -12.5 & V1 < -6.5 & V2 < -8.5)
p = ggplot(CEU_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))


EastAsian_chunk = subset(PC_all_umap, V1 >14 &  V2 < -2.5)
p = ggplot(EastAsian_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))

SouthAsian_chunk = subset(PC_all_umap, V1 >-5.5 & V1 < -1 & V2 > 12)
p = ggplot(SouthAsian_chunk,aes(x=V1,y=V2,color=Population)) + geom_point(aes(shape = Superpopulation, colour = Population)) + scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p + theme(panel.grid.minor = element_line(colour="white", size=0.5)) +  scale_y_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5)) + scale_x_continuous(minor_breaks = seq(-15 , 15, 0.5), breaks = seq(-15 , 15, 0.5))



## 6.7 Keep only the samples in your data and not in 1000G and save
American_chunk_data2 = subset(American_chunk, Dataset == "SPARK") ###input your dataname
dim(American_chunk_data2)

SouthAsian_chunk_data2 = subset(SouthAsian_chunk, Dataset == "SPARK") ###input your dataname
dim(SouthAsian_chunk_data2)

EastAsian_chunk_data2 = subset(EastAsian_chunk, Dataset == "SPARK") ###input your dataname
dim(EastAsian_chunk_data2)

CEU_chunk_data2 = subset(CEU_chunk, Dataset == "SPARK") ###input your dataname
dim(CEU_chunk_data2)

TSI_chunk_data2 = subset(TSI_chunk, Dataset == "SPARK") ###input your dataname
dim(TSI_chunk_data2)

Finnish_chunk_data2 = subset(Finnish_chunk, Dataset == "SPARK") ###input your dataname
dim(Finnish_chunk_data2)


save(American_chunk_data2, file = "American_chunk_SPARK.rda") ###input your dataname
save(TSI_chunk_data2, file = "TSI_chunk_SPARK.rda") ###input your dataname
save(SouthAsian_chunk_data2, file = "SouthAsian_chunk_SPARK.rda") ###input your dataname
save(EastAsian_chunk_data2, file = "EastAsian_chunk_SPARK.rda") ###input your dataname
save(CEU_chunk_data2, file = "CEU_chunk_SPARK.rda") ###input your dataname
save(Finnish_chunk_data2, file = "Finnish_chunk_SPARK.rda") ###input your dataname

##6.8 Create keep files
fam_check = fread("QC3output.fam")
fam_check = fam_check[,c("V1", "V2")]

keep_American = fam_check[fam_check$V2 %in% American_chunk_data2$Sample,]
keep_CEU = fam_check[fam_check$V2 %in% CEU_chunk_data2$Sample,]
keep_TSI = fam_check[fam_check$V2 %in% TSI_chunk_data2$Sample,]
keep_SouthAsian = fam_check[fam_check$V2 %in% SouthAsian_chunk_data2$Sample,]
keep_EastAsian = fam_check[fam_check$V2 %in% EastAsian_chunk_data2$Sample,]
keep_Finnish = fam_check[fam_check$V2 %in% Finnish_chunk_data2$Sample,]

write.table(keep_American, file = "keep_American.txt", row.names = F, col.names = F, quote = F)
write.table(keep_CEU, file = "keep_CEU.txt", row.names = F, col.names = F, quote = F)
write.table(keep_SouthAsian, file = "keep_SouthAsian.txt", row.names = F, col.names = F, quote = F)
write.table(keep_EastAsian, file = "keep_EastAsian.txt", row.names = F, col.names = F, quote = F)
write.table(keep_TSI, file = "keep_TSI.txt", row.names = F, col.names = F, quote = F)
write.table(keep_Finnish, file = "keep_Finnish.txt", row.names = F, col.names = F, quote = F)

### 7 Run HWE###
./plink --bfile QC3output --keep keep_American.txt --hwe 0.000001  --make-bed --out QC3_american
./plink --bfile QC3output --keep keep_CEU.txt --hwe 0.000001  --make-bed --out QC3_CEU
./plink --bfile QC3output --keep keep_TSI.txt --hwe 0.000001  --make-bed --out QC3_TSI
./plink --bfile QC3output --keep keep_SouthAsian.txt --hwe 0.000001  --make-bed --out QC3_southasian
./plink --bfile QC3output --keep keep_EastAsian.txt --hwe 0.000001  --make-bed --out QC3_eastasian
./plink --bfile QC3output --keep keep_Finnish.txt --hwe 0.000001  --make-bed --out QC3_finnish

### 8 Check Het##
for i in QC3_american QC3_CEU QC3_TSI QC3_southasian QC3_eastasian QC3_finnish; do ./plink --bfile ${i} --het --out check_het${i}; done


het_check = function (input){ 
  het = read.delim(input, sep = "")
  het$HET = (het$N.NM. - het$O.HOM.)/het$N.NM. #create heterozygosity stats
  mean = mean(het$HET)
  sd = sd(het$HET)
  het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity
  hetoutlier = subset(het, abs(Z) > 3)
  failedsample = hetoutlier[,c(1:2)]
  write.table(failedsample, file = paste0(input, "het.txt"), row.names = F, col.names = T, quote = F)
}


filelist = list("check_hetQC3_southasian.het", "check_hetQC3_CEU.het", "check_hetQC3_eastasian.het",
                "check_hetQC3_TSI.het", "check_hetQC3_finnish.het",  "check_hetQC3_american.het")

lapply(filelist,het_check)

for i in american TSI CEU southasian eastasian finnish; do ./plink --bfile QC3_${i} --remove check_hetQC3_${i}.hethet.txt --make-bed --out QC4_${i}; done

### Run geno again
./plink --bfile QC4_american --merge-list QC4_mergelist.txt --geno 0.05 --make-bed --out QC5output

##Update positions in bimfiles
bimfile = fread("QC5output.bim")
bimfile$V7 = bimfile$V2
bimfile$V7 <- gsub('GSA-', '', bimfile$V7)
update_id_file = bimfile[,c("V2", "V7")]
write.table(update_id_file, file = "update_id_file.txt", row.names = F, col.names = F, quote = F)

./plink --bfile QC5output --update-name update_id_file.txt --make-bed --out QC6output


### Step 12: Make files ready for the imputation server
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zip
unzip HRC-1000G-check-bim-v4.2.9.zip

wget http://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
gunzip 1000GP_Phase3_combined.legend.gz

./plink --bfile QC6output --freq --out QC6outputfreq



perl HRC-1000G-check-bim.pl -b QC6output.bim -f QC6outputfreq.frq -r 1000GP_Phase3_combined.legend -g 


./plink --bfile QC6output --exclude Exclude-QC6output-1000G.txt --make-bed --out TEMP1
./plink --bfile TEMP1 --update-map Chromosome-QC6output-1000G.txt --update-chr --make-bed --out TEMP2
./plink --bfile TEMP2 --update-map Position-QC6output-1000G.txt --make-bed --out TEMP3
./plink --bfile TEMP3 --flip Strand-Flip-QC6output-1000G.txt --make-bed --out TEMP4
./plink --bfile TEMP4 --reference-allele Force-Allele1-QC6output-1000G.txt --make-bed --out QC6output-updated
rm TEMP*
  
  ## Seperate into CEU and non-CEU
  
  ./plink --bfile QC6output-updated --keep keep_CEU.txt --make-bed --out QC6output-updated_CEU

for i in {1..23}; do ./plink --bfile QC6output-updated_CEU --reference-allele Force-Allele1-QC6output-1000G.txt --chr ${i} --recode-vcf --out QC6output_CEU_chr${i}; done
for i in {1..23}; do vcf-sort QC6output_CEU_chr${i}.vcf | bgzip -c > QC6output_CEU_chr${i}.vcf.gz; done
for i in {1..23}; do rm QC6output_CEU_chr${i}.vcf | rm QC6output_CEU_chr${i}.log; done


###non-CEU

keep_nonCEU = do.call("rbind", list(keep_American, keep_EastAsian, keep_Finnish, keep_SouthAsian, keep_TSI))
write.table(keep_nonCEU, file = "keep_nonCEU.txt", row.names = F, col.names = T, quote = F)


./plink --bfile QC6output-updated --keep keep_nonCEU.txt --make-bed --out QC6output-updated_nonCEU

for i in {1..23}; do ./plink --bfile QC6output-updated_nonCEU --reference-allele Force-Allele1-QC6output-1000G.txt --chr ${i} --recode-vcf --out QC6output_nonCEU_chr${i}; done
for i in {1..23}; do vcf-sort QC6output_nonCEU_chr${i}.vcf | bgzip -c > QC6output_nonCEU_chr${i}.vcf.gz; done
for i in {1..23}; do rm QC6output_nonCEU_chr${i}.vcf | rm QC6output_nonCEU_chr${i}.log; done