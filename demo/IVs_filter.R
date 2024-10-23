library(readr)
library(dplyr)
library(stringr)

set.seed(0)
trait = read.csv("./data/trait.csv") # The `trait.csv` is at the R package folder.

# data load
test = "SCZ -> MDD"
S =  stringr::str_split(test," -> ")[[1]]

risk = trait[trait$Abbreviation == tolower(S[1]),]
disease = trait[trait$Abbreviation == tolower(S[2]),]

# load data
SCZ_filepath = "./" # PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz Download from https://pgc.unc.edu/for-researchers/download-results/; PMID:35396580
MDD_filepath = "./" # PGC_UKB_depression_genome-wide.txt Download from https://pgc.unc.edu/for-researchers/download-results/; PMID:30718901
risk_data = read_delim(SCZ_filepath,delim=list(space=" ",tab="\t")[[risk$delimeter]],skip=risk$skip)
disease_data = read_delim(MDD_filepath,delim=list(space=" ",tab="\t")[[disease$delimeter]],skip=disease$skip)

# combine gwas summary by cause
risk_cols = risk[,c("rsid","beta","se","A1","A2","P.value")]
disease_cols = disease[,c("rsid","beta","se","A1","A2","P.value")]

# numeric
for(ind in c("beta","se","P.value")){
    risk_data[[risk_cols[[ind]]]] = as.numeric(risk_data[[risk_cols[[ind]]]])
    disease_data[[disease_cols[[ind]]]] = as.numeric(disease_data[[disease_cols[[ind]]]])
}

# The function of `gwas_merge` follows the merge in the CAUSE package.
X <- gwas_merge(risk_data, disease_data, snp_name_cols = c(risk_cols[[1]], disease_cols[[1]]),
                beta_hat_cols = c(risk_cols[[2]], disease_cols[[2]]),
                se_cols = c(risk_cols[[3]], disease_cols[[3]]),
                A1_cols = c(risk_cols[[4]], disease_cols[[4]]),
                A2_cols = c(risk_cols[[5]], disease_cols[[5]]),
                pval_cols = c(risk_cols[[6]], disease_cols[[6]]))

r2_thresh = 0.1
pval_thresh = 5e-8

# IVs_scz_mdd.csv
X_clump <- X[X$p1 < pval_thresh,] %>%
    rename(rsid = snp,
           pval = p1) %>%
    ieugwasr::ld_clump(dat = .,
                       clump_r2 = r2_thresh,
                       clump_p = pval_thresh,
                       plink_bin = "/home/tb/tools/GWAS//plink",
                       bfile="/home/tb/ref/1kg.v3//EUR")

# initEst_scz_mdd.csv
X1 = X[X$p1 > 0.5,]
X_clump1 <- X1[sample(seq(dim(X1)[1]),100000),] %>%
    rename(rsid = snp,
           pval = p1) %>%
    ieugwasr::ld_clump(dat = .,
                       clump_r2 = r2_thresh,
                       clump_p = 1,
                       plink_bin = "/home/tb/tools/GWAS//plink",
                       bfile="/home/tb/ref/1kg.v3/EUR")
