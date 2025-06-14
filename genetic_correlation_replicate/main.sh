cd /projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/genetic_correlation_replicate

export LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'
DATA='/projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/data'
module load R
module load python2
module load ldsc
module load htslib/1.16
module load bcftools/1.16

# endo
zcat $DATA/GCST90205183_buildGRCh37.tsv.gz | head
# uterine fibroids
zcat $DATA/GCST90468154.tsv.gz | head
# asthma 
zcat $DATA/ukb.childasthma.upload.final.assoc.gz | head 
# migraine
zcat $DATA/GCST90129450_buildGRCh37.tsv.gz | head
# osteoarthritis
zcat $DATA/GCST90129444_buildGRCh37.tsv.gz | head
# dysmenorrhea
zcat $DATA/GCST90044465_buildGRCh37.tsv.gz | head
# Menorrhagia
zcat $DATA/GCST90454221.tsv.gz | head

ls $LDSCORES_DIR/eur_w_ld_chr
ls $LDSCORES_DIR/UKBB.ALL.ldscore
ls $LDSCORES_DIR/w_hm3.snplist


# GCh38 data chrpos_to_rsid.txt
head ../genetic_correlation/chrpos_to_rsid_standard.txt

awk -F '[:\t]' '
{
    # Remove "NC_0000" prefix and ".11" suffix from chromosome ID
    if ($1 ~ /^NC_0000[0-9]+/) {
        chr_num = substr($1, 8, 3) + 0
        print chr_num "\t" $2 "\t" $3
    } else if ($1 == "NC_000023.11") {
        print "X\t" $2 "\t" $3
    } else if ($1 == "NC_000024.11") {
        print "Y\t" $2 "\t" $3
    } else if ($1 == "NC_012920.1") {
        print "MT\t" $2 "\t" $3
    }
}' ../genetic_correlation/chrpos_to_rsid_standard.txt > chrpos_to_rsid.txt 
head chrpos_to_rsid.txt


##############
### UKBB use to map rsIDs
##############
zcat /projectnb/bs859/data/ldscore_files/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz | awk 'NR > 1 {print $1 "\t" $3 "\t" $2}' > ukbb_chr_bp_rsid.txt


############################
### ENDOMETRIOSIS ########## 
############################
awk '
    BEGIN { FS=OFS="\t" }
    FNR==NR {
        key = $1 ":" $2
        map[key] = $3
        next
    }
    NR==1 {
        print "rsID", $0
        next
    }
    {
        key = $1 ":" $2
        rsid = (key in map ? map[key] : "NA")
        print rsid, $0
    }
' ukbb_chr_bp_rsid.txt <(zcat $DATA/GCST90205183_buildGRCh37.tsv.gz) > full_endo_data_with_rsid.tsv

awk 'NR==1 || $1 != "NA"' full_endo_data_with_rsid.tsv > endo_with_rsid_only.tsv
wc -l endo_with_rsid_only.tsv # 1,085,839 SNPs


##Step 1:  format summary statistics for ldsc
munge_sumstats.py \
--sumstats endo_with_rsid_only.tsv \
--snp NA \
--a1 effect_allele \
--a2 other_allele \
--N 470866 \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out ENDO_fulldata
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  ENDO_fulldata.sumstats.gz|head -n 5
zcat  ENDO_fulldata.sumstats.gz|wc # 1,217,312 variants
## Step 2: run the LD score regression
ldsc.py \
--h2 ENDO_fulldata.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out ENDO_fulldata_h2_orig_scale_UKBB

cat ENDO_fulldata_h2_orig_scale_UKBB.log | tail -n 8
# Using two-step estimator with cutoff at 30.
# Total Observed scale h2: 0.0147 (0.0014)
# Lambda GC: 1.0405
# Mean Chi^2: 1.0917
# Intercept: 0.9464 (0.0085)
# Ratio < 0 (usually indicates GC correction).
# Analysis finished at Sun Apr 20 22:26:14 2025
# Total time elapsed: 41.12s

##################### 
###### asthma #######
##################### GCST008917

zcat $DATA/ukb.childasthma.upload.final.assoc.gz | head -n 2

# FIXME: column names
##Step 1:  format summary statistics for ldsc
munge_sumstats.py \
--sumstats $DATA/ukb.childasthma.upload.final.assoc.gz \
--snp SNP \
--a1 A1 \
--a2 A2 \
--N 357157 \
--signed-sumstats BETA,0 \
--p P \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out ATHSMA
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  ATHSMA.sumstats.gz|head -n 5
zcat  ATHSMA.sumstats.gz|wc # 1217312 variants

##################### 
###### uterine fibroid #######
#####################
# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468154/GCST90468154.tsv.gz

zcat $DATA/GCST90468154.tsv.gz | head # 13,308,323 variants, GCh37


##Step 1:  format summary statistics for ldsc
munge_sumstats.py \
--sumstats $DATA/GCST90468154.tsv.gz \
--snp rs_id \
--a1 other_allele \
--a2 effect_allele \
--N-col n \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out UTERINE_FIBROIDS
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  UTERINE_FIBROIDS.sumstats.gz|head -n 5
zcat  UTERINE_FIBROIDS.sumstats.gz|wc # 1,217,312 variants
## Step 2: run the LD score regression
ldsc.py \
--h2 UTERINE_FIBROIDS.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out UTERINE_FIBROIDS_h2_orig_scale_UKBB

cat UTERINE_FIBROIDS_h2_orig_scale_UKBB.log | tail -n 8
# Using two-step estimator with cutoff at 30.
# Total Observed scale h2: 0.0147 (0.0014)
# Lambda GC: 1.0405
# Mean Chi^2: 1.0917
# Intercept: 0.9464 (0.0085)
# Ratio < 0 (usually indicates GC correction).
# Analysis finished at Tue Apr 22 18:57:18 2025
# Total time elapsed: 43.34s



########################################## 
### migraine MIGRAINE ####################
########################################## 

# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129450/GCST90129450_buildGRCh37.tsv.gz

zcat $DATA/GCST90129450_buildGRCh37.tsv.gz | head

##Step 1:  format summary statistics for ldsc
munge_sumstats.py \
--sumstats $DATA/GCST90129450_buildGRCh37.tsv.gz \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--N 211460 \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out MIGRAINE
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  MIGRAINE.sumstats.gz|head -n 5
zcat  MIGRAINE.sumstats.gz|wc # 1,217,312 variants
## Step 2: run the LD score regression
ldsc.py \
--h2 MIGRAINE.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out MIGRAINE_h2_orig_scale_UKBB

cat MIGRAINE_h2_orig_scale_UKBB.log | tail -n 8
# Using two-step estimator with cutoff at 30.
# Total Observed scale h2: 0.0346 (0.003)
# Lambda GC: 1.1683
# Mean Chi^2: 1.2041
# Intercept: 1.0504 (0.0093)
# Ratio: 0.247 (0.0453)
# Analysis finished at Wed Apr 23 17:49:53 2025
# Total time elapsed: 40.63s



########################################## 
### osteoarthritis, hip #################### OSTEOARTHRITIS
########################################## 
# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129444/GCST90129444_buildGRCh37.tsv.gz
# osteoarthritis
zcat $DATA/GCST90129444_buildGRCh37.tsv.gz | head

munge_sumstats.py \
--sumstats $DATA/GCST90129444_buildGRCh37.tsv.gz \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--N 210724 \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out OSTEOARTHRITIS
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  OSTEOARTHRITIS.sumstats.gz|head -n 5
zcat  OSTEOARTHRITIS.sumstats.gz|wc # 1,217,312 variants


##################### 
###### dysmenorrhea/painful mentration ####### DYSMENORRHEA
#####################
# wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90084001-GCST90085000/GCST90084673/GCST90084673_buildGRCh38.tsv.gz
# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90478001-GCST90479000/GCST90478722/GCST90478722.tsv.gz # GRCh37

# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90044001-GCST90045000/GCST90044465/GCST90044465_buildGRCh37.tsv.gz
# dysmenorrhea
zcat $DATA/GCST90044465_buildGRCh37.tsv.gz | head
munge_sumstats.py \
--sumstats $DATA/GCST90044465_buildGRCh37.tsv.gz \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--N-col N \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out DYSMENORRHEA
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  DYSMENORRHEA.sumstats.gz|head -n 5
zcat  DYSMENORRHEA.sumstats.gz|wc # 1,217,312 variants


##################### 
###### Menorrhagia ####### MENORRHAGIA
#####################
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454221/GCST90454221.tsv.gz
# Menorrhagia
zcat $DATA/GCST90454221.tsv.gz | head

munge_sumstats.py \
--sumstats $DATA/GCST90454221.tsv.gz \
--snp rs_id \
--a1 effect_allele \
--a2 other_allele \
--N 249874 \
--p p_value \
--signed-sumstats beta,0 \
--merge-alleles $LDSCORES_DIR/w_hm3.snplist \
--out MENORRHAGIA
#output:  reformated, gzipped summary statistics (in your own directory)
zcat  MENORRHAGIA.sumstats.gz|head -n 5
zcat  MENORRHAGIA.sumstats.gz|wc # 1,217,312 variants



### ### ### ### ### ### ### 
### Genetic Correlation ### 
### ### ### ### ### ### ### 
## estimate the genetic correlations between uterine fibrosis, migrains, and endometriosis

## Step 3: genetic correlations
ldsc.py \
--rg ENDO_fulldata.sumstats.gz,ATHSMA.sumstats.gz,UTERINE_FIBROIDS.sumstats.gz,MIGRAINE.sumstats.gz,OSTEOARTHRITIS.sumstats.gz,DYSMENORRHEA.sumstats.gz,MENORRHAGIA.sumstats.gz \
--ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out ENDO_allTRAITS.HGT_rg
# look at results
cat ENDO_allTRAITS.HGT_rg.log | tail -n 10














##################### 
### athsma
# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007800/CHILD_ONSET_ASTHMA.20180501.allchr.assoc.GC.gz	
##################### 

# zcat $DATA/CHILD_ONSET_ASTHMA.20180501.allchr.assoc.GC.gz	 | head -n 3


# ##Step 1:  format summary statistics for ldsc
# munge_sumstats.py \
# --sumstats $DATA/ukb.childasthma.upload.final.assoc.gz \
# --snp SNP \
# --a1 A1 \
# --a2 A2 \
# --N 9676 \
# --p P \
# --signed-sumstats BETA,0 \
# --merge-alleles $LDSCORES_DIR/w_hm3.snplist \
# --out ATHSMA
# #output:  reformated, gzipped summary statistics (in your own directory)
# zcat  ATHSMA.sumstats.gz|head -n 10
# zcat  ATHSMA.sumstats.gz|wc # 1,217,312 variants
# ## Step 2: run the LD score regression
# ldsc.py \
# --h2 ATHSMA.sumstats.gz \
# --ref-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --w-ld $LDSCORES_DIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --out ATHSMA_h2_orig_scale_UKBB

# cat ATHSMA_h2_orig_scale_UKBB.log | tail -n 8
## bad results here??
# Using two-step estimator with cutoff at 30.
# Total Observed scale h2: 1.0501 (0.1507)
# Lambda GC: 1.1459
# Mean Chi^2: 1.2738
# Intercept: 1.0686 (0.0118)
# Ratio: 0.2507 (0.043)
# Analysis finished at Tue Apr 22 19:22:33 2025
# Total time elapsed: 42.1s




# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468152/GCST90468152.tsv.gz

