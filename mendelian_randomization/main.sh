module load gcta/1.94.1
module load R
DATA='/projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/data'
cd /projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/mendelian_randomization


##reformatting the file: SNP A1 A2 freq b se p N
##########
## endo ##
##########
cat /projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/genetic_correlation_replicate/endo_with_rsid_only.tsv | head
# add N column = 470866
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="N"} NR>1{$(NF+1)=470866} 1' \
/projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/genetic_correlation_replicate/endo_with_rsid_only.tsv \
> endo_with_rsid_only_with_N.tsv
head endo_with_rsid_only_with_N.tsv
awk 'NR==1{print "SNP A1 A2 freq b se p N"};(NR>1&&$1!="."){print $1,$4,$5,$6,$7,$8,$9,$10}' endo_with_rsid_only_with_N.tsv > endo.ss.txt
head endo.ss.txt
sort -t ' ' -k 1,1 -u endo.ss.txt > endo.ss.dedup.txt
wc endo.ss.txt endo.ss.dedup.txt
mv endo.ss.dedup.txt endo.ss.txt
#   1085840   8686720  53945829 endo.ss.txt
#   1085840   8686720  53945829 endo.ss.dedup.txt
#   2171680  17373440 107891658 total


######################
## uterine fibroids ##
######################
zcat $DATA/GCST90468154.tsv.gz | head

zcat $DATA/GCST90468154.tsv.gz | \
awk -F'\t' 'BEGIN {OFS="\t"} 
NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"} 
NR>1 {print $10, $3, $4, $7, $5, $6, $8, $11}' > uterine.ss.txt
head uterine.ss.txt
sort -t ' ' -k 1,1 -u uterine.ss.txt > uterine.ss.dedup.txt
wc uterine.ss.txt uterine.ss.dedup.txt
mv uterine.ss.dedup.txt uterine.ss.txt
# 13308323  106466584  825651794 uterine.ss.txt
# 13308048  106464384  825636540 uterine.ss.dedup.txt
# 26616371  212930968 1651288334 total



############
## asthma ## no freq???
############
zcat $DATA/ukb.childasthma.upload.final.assoc.gz | head 


############
# migraine #
############
zcat $DATA/GCST90129450_buildGRCh37.tsv.gz | head
# add N column = 211460
zcat $DATA/GCST90129450_buildGRCh37.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="N"} NR>1{$(NF+1)=211460} 1' \
| gzip > GCST90129450_with_N.tsv.gz
# check file
zcat GCST90129450_with_N.tsv.gz | head
# make gcta file
zcat GCST90129450_with_N.tsv.gz | \
awk -F'\t' 'BEGIN {OFS="\t"} 
NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"} 
NR>1 {print $3,$5,$4,$6,$7,$8,$9,$11}' > migraine.ss.txt
head migraine.ss.txt

sort -t ' ' -k 1,1 -u migraine.ss.txt > migraine.ss.dedup.txt
wc migraine.ss.txt migraine.ss.dedup.txt
  # 11323612   90588896  783061177 migraine.ss.txt
  # 11323612   90588896  783061177 migraine.ss.dedup.txt
  # 22647224  181177792 1566122354 total
mv migraine.ss.dedup.txt migraine.ss.txt
##################
# osteoarthritis #
##################
zcat $DATA/GCST90129444_buildGRCh37.tsv.gz | head
# N col = 210724
zcat $DATA/GCST90129444_buildGRCh37.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="N"} NR>1{$(NF+1)=210724} 1' \
| gzip > GCST90129444_with_N.tsv.gz
# check
zcat GCST90129444_with_N.tsv.gz | head
# make file
zcat GCST90129444_with_N.tsv.gz | \
awk -F'\t' 'BEGIN {OFS="\t"} 
NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"} 
NR>1 {print $3,$5,$4,$6,$7,$8,$9,$11}' > osteoarthritis.ss.txt
head osteoarthritis.ss.txt

sort -t ' ' -k 1,1 -u osteoarthritis.ss.txt > osteoarthritis.ss.dedup.txt
wc osteoarthritis.ss.txt osteoarthritis.ss.dedup.txt
mv osteoarthritis.ss.dedup.txt osteoarthritis.ss.txt
  # 11323612   90588896  782508355 osteoarthritis.ss.txt
  # 11323612   90588896  782508355 osteoarthritis.ss.dedup.txt
  # 22647224  181177792 1565016710 total
##################
## dysmenorrhea ##
##################
zcat $DATA/GCST90044465_buildGRCh37.tsv.gz | head
# make file
zcat $DATA/GCST90044465_buildGRCh37.tsv.gz | \
awk -F'\t' 'BEGIN {OFS="\t"}
NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"}
NR>1 {print $2, $4, $5, $7, $12, $13, $14, $6}' > dysmenorrhea.ss.txt
head dysmenorrhea.ss.txt
sort -t ' ' -k 1,1 -u dysmenorrhea.ss.txt > dysmenorrhea.ss.dedup.txt
wc dysmenorrhea.ss.txt dysmenorrhea.ss.dedup.txt
mv dysmenorrhea.ss.dedup.txt dysmenorrhea.ss.txt
  # 11832030   94656240  615687107 dysmenorrhea.ss.txt
  # 11832030   94656240  615687107 dysmenorrhea.ss.dedup.txt
  # 23664060  189312480 1231374214 total
##################
## menorrhagia ##
##################
zcat $DATA/GCST90454221.tsv.gz | head
# N = 249874
zcat $DATA/GCST90454221.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="N"} NR>1{$(NF+1)=249874} 1' \
| gzip > GCST90454221_with_N.tsv.gz
zcat GCST90454221_with_N.tsv.gz | head
# make file
zcat GCST90454221_with_N.tsv.gz | \
awk -F'\t' 'BEGIN {OFS="\t"}
NR==1 {print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"}
NR>1 {print $9, $3, $4, $7, $5, $6, $8, $12}' > menorrhagia.ss.txt
head menorrhagia.ss.txt
sort -t ' ' -k 1,1 -u menorrhagia.ss.txt > menorrhagia.ss.dedup.txt
wc menorrhagia.ss.txt menorrhagia.ss.dedup.txt
mv menorrhagia.ss.dedup.txt menorrhagia.ss.txt
  # 18299836  146398688 1076369353 menorrhagia.ss.txt
  # 18299836  146398688 1076369353 menorrhagia.ss.dedup.txt
  # 36599672  292797376 2152738706 total


#########################################################################
######################## Mendelian randomization ########################
#########################################################################
### dysmenorrhea
###run gsmr to estimate causal effect of all exposures on endometriosis 
##write the exposure and outcome file names to files for gsmr to read:
echo "dysmenorrhea dysmenorrhea.ss.txt" > exposure.txt
echo "endo endo.ss.txt" > outcome.txt
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed to eliminate SNPs that are in high LD)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-file exposure.txt outcome.txt --gwas-thresh 8e-4 --effect-plot --gsmr-direction 0 --out ENDO_dysmenorrhea-gsmr 
# p < 8.0e-04 and LD r2 < 0.05.
# $ cat ENDO_dysmenorrhea-gsmr.gsmr
# Exposure        Outcome bxy     se      p       nsnp
# dysmenorrhea    endo    0.574812        0.0428248       4.46917e-41     17


### menorrhagia
###run gsmr to estimate causal effect of all exposures on endometriosis 
##write the exposure and outcome file names to files for gsmr to read:
awk '{
  count[$1]++
  if (count[$1] > 1)
    $1 = $1"_"count[$1]-1
  print
}' menorrhagia.ss.txt > tempfile_menorrhagia.ss.txt
echo "menorrhagia tempfile_menorrhagia.ss.txt" > exposure.txt
echo "endo endo.ss.txt" > outcome.txt
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed to eliminate SNPs that are in high LD)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-file exposure.txt outcome.txt --gwas-thresh 5e-7 --effect-plot --gsmr-direction 0 --out ENDO_menorrhagia-gsmr
# p < 5.0e-07 and LD r2 < 0.05
# $ cat ENDO_menorrhagia-gsmr.gsmr
# Exposure        Outcome bxy     se      p       nsnp
# menorrhagia     endo    0.914367        0.0756137       1.15549e-33     15

##################
### migraine
###run gsmr to estimate causal effect of all exposures on endometriosis 
##write the exposure and outcome file names to files for gsmr to read:
echo "migraine migraine.ss.txt" > exposure.txt
echo "endo endo.ss.txt" > outcome.txt
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed to eliminate SNPs that are in high LD)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-file exposure.txt outcome.txt --effect-plot --gsmr-direction 0 --out ENDO_migraine-gsmr
# $ head ENDO_migraine-gsmr.gsmr 
# Exposure        Outcome bxy     se      p       nsnp
# migraine        endo    0.0685394       0.0415356       0.0989148       13


### osteoarthritis
###run gsmr to estimate causal effect of all exposures on endometriosis 
##write the exposure and outcome file names to files for gsmr to read:
echo "osteoarthritis osteoarthritis.ss.txt" > exposure.txt
echo "endo endo.ss.txt" > outcome.txt
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed to eliminate SNPs that are in high LD)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-file exposure.txt outcome.txt --effect-plot --gsmr-direction 0 --out ENDO_osteoarthritis-gsmr
# $ head ENDO_osteoarthritis-gsmr.gsmr 
# Exposure        Outcome bxy     se      p       nsnp
# osteoarthritis  endo    0.0292155       0.0362141       0.419815        17


### uterine
###run gsmr to estimate causal effect of all exposures on endometriosis 
##write the exposure and outcome file names to files for gsmr to read:
# grep -v 'NA' uterine.ss.txt > tempfile_uterine.ss.txt
echo "uterine tempfile_uterine.ss.txt" > exposure.txt
echo "endo endo.ss.txt" > outcome.txt
##Use 1000 Genomes Europeans to estimate LD among the SNPs (this is needed to eliminate SNPs that are in high LD)
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --gsmr-file exposure.txt outcome.txt --gwas-thresh 5e-4 --effect-plot --gsmr-direction 0 --out ENDO_uterine-gsmr
# 140 SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample. These SNPs have been saved in [ENDO_uterine-gsmr.freq.badsnps].
# Error: there are too many SNPs that have large difference in allele frequency. Please check the GWAS summary data.
# I manually checked the SNP rsID and the coordinates match up to the GCh37 like the rest of the data.
# Usually this error is due the inconsistent definition of allele frequency in user's own GWAS summary even the allele frequency column seems to be correct.
# I would suggest you generate the freq files from the LD reference file you provided using the PLINK function. And then compare the freq file with your own GWAS.


