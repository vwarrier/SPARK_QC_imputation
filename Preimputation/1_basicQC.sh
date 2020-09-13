###Basic QC

./plink --bfile SPARK.30K.array_genotype.20190829.201911update --geno 0.1  --make-bed --out QC1output

./plink --bfile QC1output --make-bed --me 0.05 0.1 --mind 0.05 --out QC2output

./plink --bfile QC2output --make-bed --out QC3output --remove failedsample.txt


