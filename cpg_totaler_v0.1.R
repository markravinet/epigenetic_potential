# cpg_totaler
# a lighterweight counting script to identify the number of total cpgs on all sites
# v 0.1 - 15/08/2024

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(getopt))

# specify command line options
spec <- matrix(c(
  'vcf', 'v', 1, 'character', 'vcf path',
  'flanking', 'f', 1, 'character', 'Path to flanking sites info - as bed file',
  'gff', 'g', 1, 'character', 'Path to gff annotation',
  'out', 'o', 1, 'character', 'Output file'
), ncol = 5, byrow = T)

# set command line options
opt = getopt(spec)

# show help if asked for
if (!is.null(opt$help)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# set variables for call
vcf_path <- opt$vcf
flanking_path <- opt$flanking
gff_path <- opt$gff
out <- opt$out

# testing data
# flanking_path <- "./data/22_cpg_flanking_sites.bed" # flanking sites
# vcf_path <- "./data/22_norm_filtered_gs.vcf.gz" # vcf
# gff_path <- "./Ficedula_albicollis.FicAlb1.5.109.gff3.gz" # gff
# out <- "./data/22_cpg.counts.gz"

# read in vcf
vcf <- read.vcfR(vcf_path)
# read in gff
gff <- read.table(gff_path, sep = "\t", quote = "")
# filter gff - assumes per chr vcf
my_chr <- vcf@fix[1,1] %>% as.vector()

### Dealing with flanking sites ###
# read in flanking sites
flanking <- read_tsv(flanking_path, col_names = c("pos", "sites"))
# derive sites from  vcf
sites <- vcf@fix[, c(1:2, 4, 5)] %>% as_tibble()
names(sites) <- tolower(names(sites))
sites <- sites %>% mutate(pos = as.numeric(pos))

# bind columns
flanking <- bind_cols(sites, flanking)

# split sites
s1 <- str_sub(flanking$sites, 1, 1)
s2 <- str_sub(flanking$sites, 2, 2)
s3 <- str_sub(flanking$sites, 3, 3)

# create a new tibble
cpg <- bind_cols(select(flanking, chr = chrom, pos = pos...2, ref, alt), as_tibble(data.frame(s1, s2, s3)))

# cpg1 is a site where the reference is C and the alternative is another base, because ref is always 0, these need recoding
# cpg2 is a site where the reference is G and the alternative is another base, because ref is always 0, these need recoding
# poly-cpg3 is a site where the alternative is C - this does not need recoding as ALT alleles are coded as 1
# poly-cpg4 is a site where the alternative is C - this does not need recoding as ALT alleles are coded as 1


# classify sites
status <- rep(NA, nrow(cpg))
# classify cpg - simplest case - focal site
status[which(s2 == "C" & s3 =="G")] <- "cpg1"
# classify cpg - additional case - first site
status[which(s1 == "C" & s2 =="G")] <- "cpg2"
# nb - remove gpc to prevent confusion - 26/04/23
# classify gpc
#status[which(s2 == "G" & s3 =="C")] <- "gpc"
# also need to incorporate potentially polymorphic sites - i.e. where alt allele is a C - i.e. poly_cpg
status[which(cpg$alt == "C" & cpg$s3 == "G")] <- "poly-cpg3"
# need to also extend this to sites where the first flanking allele is C (and the focal alt is G)
status[which(cpg$alt == "G" & cpg$s1 == "C")] <- "poly-cpg4"

### For CpG sites in REF
# create cpg 
cpg_vcf <- vcf[which(status == "cpg1" | status == "cpg2"), ]

# create vcf for testing
test_vcf <- cpg_vcf[, ]

# get counts
counts <- extract_gt_tidy(test_vcf, alleles = FALSE) %>% 
  separate(col = gt_GT, into = c("a1", "a2"), convert = T) %>%
  mutate(count = (a1 + a2)) 
# convert counts (cpg1 and cpg2 only)
new_count <- counts$count 
new_count[counts$count == 0] <- 2
new_count[counts$count == 2] <- 0
counts$count <- new_count

# convert back to a vcf style matrix
counts2 <- counts %>%
  select(Key, Indiv, count) %>%
  #select(-a1, -a2) %>%
  pivot_wider(names_from = Indiv, values_from = count)

# run on a normal vcf means that this produces additional columns for GT values - i.e. DP, AD, GQ
# writing these out within the same output is complicated (massive files) - so option is to extract these and 
# write them out as a separate file - can add this functionality later

### For CpG sites in ALT - i.e. polycpg
# create cpg 
cpg_vcf2 <- vcf[which(status == "poly-cpg3" | status == "poly-cpg4"), ]

# create vcf for testing
test_vcf2 <- cpg_vcf2[, ]

# get counts
counts <- extract_gt_tidy(test_vcf2, alleles = FALSE) %>% 
  separate(col = gt_GT, into = c("a1", "a2"), convert = T) %>%
  mutate(count = (a1 + a2)) 

# convert back to a vcf style matrix
counts3 <- counts %>%
  select(Key, Indiv, count) %>%
  #select(-a1, -a2) %>%
  pivot_wider(names_from = Indiv, values_from = count)


## combine count tables
final_counts <- bind_cols(rbind(test_vcf@fix[, 1:2], test_vcf2@fix[, 1:2]),
                          bind_rows(counts2, counts3)) %>% 
  mutate(status = c(rep("cpg", nrow(counts2)), rep("polycpg", nrow(counts3)))) %>%
  mutate(POS = as.numeric(POS)) %>%
  select(-Key) %>%
  arrange(POS)

# write out total per individual
# cpg total
cpg_total <- final_counts %>%
  filter(status == "cpg") %>%
  select(-CHROM, -POS, -status) %>%
  colSums(na.rm = T)
# polycpg total
polycpg_total <- final_counts %>%
  filter(status == "polycpg") %>%
  select(-CHROM, -POS, -status) %>%
  colSums(na.rm = T)

# final table
cpg_totals <- bind_rows(cpg_total, polycpg_total)
# add columns
cpg_totals$chr <- rep(my_chr, nrow(cpg_totals))
cpg_totals$status <- c("cpg", "polycpg")

# write out
# counts
write_csv(final_counts, out)
# totals
path <- out %>% sub("\\.gz$", "", .)
out2 <- paste0(path, ".total")
# final totals
write_csv(cpg_totals, out2)
