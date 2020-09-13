
./plink --bfile QC5output --freq --out QC5outputfreq

perl HRC-1000G-check-bim.pl -b QC5output.bim -f QC5outputfreq.frq -r 1000GP_Phase3_combined.legend -g