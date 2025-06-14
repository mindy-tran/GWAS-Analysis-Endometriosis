

cd /projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/conditional_analysis_replicate

##load plink and gcta modules
module load gcta/1.94.1

module load python2/2.7.16
module load R/3.5.1
module load new_fugue/2010-06-02
module load plink/1.07
module load samtools/1.10
module load locuszoom/1.4


DATA='/projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/data'
zcat $DATA/GCST90205183_buildGRCh37.tsv.gz | head

### genome wide signif + add N column 470866
zcat $DATA/GCST90205183_buildGRCh37.tsv.gz |awk 'NR==1||$8<5e-8 {print $0"  470866"}'  > tophits.txt
wc -l tophits.txt # 863 tophits.txt
sort -k 8g tophits.txt | head # top hits are in chr1
# chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency beta    standard_error  p_value  470866
# 1       22462111        A       G       0.1936  0.1593  0.0148  5.118E-27  470866
# 1       22465820        T       C       0.8067  -0.1591 0.0148  5.926E-27  470866
# 1       22468215        T       C       0.1930  0.1588  0.0148  7.381E-27  470866
# 1       22450487        T       C       0.1933  0.1571  0.0147  1.17E-26  470866
# 1       22436446        C       G       0.1992  0.1519  0.0145  1.115E-25  470866
# 1       22470407        T       C       0.1871  0.1577  0.0151  1.566E-25  470866
# 1       22470451        C       G       0.1871  0.1576  0.0151  1.679E-25  470866
# 1       22422721        A       G       0.1977  0.1505  0.0145  3.08E-25  470866
# 1       22403357        A       G       0.2079  0.1412  0.0142  2.688E-23  470866




# get chr:pos --> rsID
awk 'BEGIN{OFS="\t"} {print $1 ":" $4, $2}' /projectnb/bs859/data/1000G/plinkformat/1000G_EUR.bim > chrpos_to_rsid.txt

### create reformatted file for GCTA conditional analysis: SNP A1 A2 freq b se p N 
zcat $DATA/GCST90205183_buildGRCh37.tsv.gz | \
awk 'NR==1 || $8 < 5e-8 {print $0"\t470866"}' | \
awk 'BEGIN {OFS="\t"} NR==1 {next} {print $1":"$2, $3, $4, $5, $6, $7, $8, $9}' > togcta.txt
head togcta.txt


awk 'FNR==NR {map[$1]=$2; next} {if ($1 in map) $1=map[$1]; print}' chrpos_to_rsid.txt togcta.txt > togcta_with_rsid.txt


# run the conditional analysis for chromosome 1 only.  This will take a while!
gcta64 --bfile /projectnb/bs859/data/1000G/plinkformat/1000G_EUR --cojo-file togcta_with_rsid.txt --cojo-slct --chr 1 --out chr1_EUR > chr1_EUR.log 
wc chr1_EUR.jma.cojo # 2 independent asociations identified with COJO.




#################
### Locuszoom ###
#################
# Use the SNP column to extract variants from your PLINK files (--extract) or plot LD regions (e.g., LocusZoom)
cut -f1,2,3,8 chr1_EUR.cma.cojo > locuszoom_input.txt # CHR    SNP       POS         P

cut -f1,2,3,8 chr1_EUR.jma.cojo | head # rs12037376 and rs484686

## manually make locuszoom_input_rsid.txt
head locuszoom_input.txt


## rs12037376
locuszoom \
  --metal locuszoom_input.txt \
  --markercol SNP \
  --pvalcol p \
  --refsnp rs12037376 \
  --chr 1 \
  --start 22462111 \
  --source 1000G_Nov2014 \
    --build hg19 \
    --pop EUR \
  --plotonly

## rs484686
locuszoom \
  --metal locuszoom_input.txt \
  --markercol SNP \
  --pvalcol p \
  --refsnp rs484686 \
  --chr 1 \
  --start 172152202 \
  --source 1000G_Nov2014 \
    --build hg19 \
    --pop EUR \
  --plotonly



###############
##### VEP #####
###############

##the top locus is on chromosome 1, has multiple variants with p<5e-20
# we know this from the tophits.txt and Fig. 1. Circular Manhattan plots in paper
zcat $DATA/GCST90205183_buildGRCh37.tsv.gz | awk 'NR==1||$8<5e-10 {print $0}'  > p5e-10.txt

## extract the rsids from the variants
## we are going to use VEP (Variant Effect Predictor) to get
## annotations for these 1448 variants.
cut -f2 p5e-10.txt > toVEP.txt 

### VEP results
# rs12037376 is WNT4, protein-coding, intron variant.
# rs484686 is DNM3, protein-coding, intron variant.




