# cpg_counter
# script to identify and count cpg sites from a vcf
# v 0.1 - 11/07/2023
# issues fixed - 17/07/2023
# v 0.2 - 19/07/2023 - updated to count sites which are not CpGs
# v 0.3 - 18/03/2025 - updated to incorporate windows
# v 0.4 - 19/11/2025 - updated to incorporate windows
# v 0.5 - 05/01/2026 - updated to incorporate windows

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(getopt))

# specify command line options
spec <- matrix(c(
  'vcf', 'v', 1, 'character', 'vcf path',
  'start', 's', 1, 'character', 'vcf start position',
  'stop', 'p', 1, 'character', 'vcf stop position',
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
start <- as.numeric(opt$start)
stop <- as.numeric(opt$stop)
flanking_path <- opt$flanking
gff_path <- opt$gff
out <- opt$out

# testing data
# flanking_path <- "./data/1_1_5000000_cpg_flanking_sites.bed" # flanking sites
# vcf_path <- "./data/1_norm_filtered_gs_depth_filtered_1_1_5000000.vcf.gz" # vcf
# gff_path <- "./Ficedula_albicollis.FicAlb1.5.109.gff3.gz" # gff
# out <- "./data/1_cpg.counts.gz"
# start <- 1
# stop <- 5000000

# read in vcf
vcf <- read.vcfR(vcf_path)
# read in gff
gff <- read.table(gff_path, sep = "\t", quote = "")
# filter gff - assumes per chr vcf
my_chr <- vcf@fix[1,1] %>% as.vector()

sprintf("Working on chromosome %s, windows:%s-%s", my_chr, start, stop)

# Filter the vcf
# split into specified window
# pull site positions from vcf and make a numeric vector
positions <- as.numeric(vcf@fix[,2])
# find positions which fulfil window requirements
idx <- which(positions >=start & positions <=stop)

# filter vcf
vcf <- vcf[idx, ]

### Dealing with flanking sites ###
# read in flanking sites
flanking <- read_tsv(flanking_path, col_names = c("pos", "sites"))
# derive sites from  vcf
sites <- vcf@fix[, c(1:2, 4, 5)] %>% as_tibble()
names(sites) <- tolower(names(sites))
sites <- sites %>% mutate(pos = as.numeric(pos))

# bind columns
flanking <- bind_cols(sites, flanking[idx,])

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

# output some status stats
x <- status
x[is.na(x)] <- "variant"
x <- as.factor(x)
z <- x %>% table()
sprintf("There are %s CpGs, %s polyCpGs and %s variants in this window and %s total variable sites", 
        z[1]+z[2], z[3] + z[4], z[5], sum(z))

### Sorting out duplicates
cpg <- as_tibble(data.frame(cpg,status))
dist_diff <- c(NA, diff(cpg$pos))
cpg$dist_diff <- dist_diff == 1

# If the previous site is a G or polymorphic for a G (i.e. polycpg) you can ignore it
# set an index for where sites are consecutive (just to make code below clearer)
x <- which(cpg$dist_diff == TRUE)
# run through dataframe and identify consecutive sites where previous alt or ref is a G (i.e. we can ignore)
y <- sapply(2:nrow(cpg), function(x) {
  cpg$dist_diff[x] == TRUE & (cpg$ref[x-1] == "G" | cpg$alt[x-1] == "G")
}, simplify = TRUE) 
# nb - you have to append an NA here - i.e. skipping first row
isG <- c(NA, y)
# add to data.frame
cpg$isG <- isG

# run through dataframe and identify consecutive sites where previous alt or ref is a G (i.e. we can ignore)
y <- sapply(2:nrow(cpg), function(x) {
  cpg$dist_diff[x] == TRUE & (cpg$ref[x-1] == "C" | cpg$alt[x-1] == "C")
}, simplify = TRUE) 
# nb - you have to append an NA here - i.e. skipping first row
isC <- c(NA, y)
# add to data.frame
cpg$isC <- isC

# do we need to adjust?
cpg$filter <- cpg$dist_diff == TRUE & cpg$isC == TRUE & !is.na(cpg$status)

### For CpG sites in REF
# create cpg 
cpg_vcf <- vcf[which(status[idx] == "cpg1" | status[idx] == "cpg2"), ]

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
cpg_vcf2 <- vcf[which(status[idx] == "poly-cpg3" | status[idx] == "poly-cpg4"), ]

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

### now also create an output that is for sites which are not cpgs
non_vcf <- vcf[which(is.na(status[idx])), ]

# get counts
counts <- extract_gt_tidy(non_vcf, alleles = FALSE) %>% 
  separate(col = gt_GT, into = c("a1", "a2"), convert = T) %>%
  mutate(count = (a1 + a2)) 

# convert back to a vcf style matrix
counts4 <- counts %>%
  select(Key, Indiv, count) %>%
  #select(-a1, -a2) %>%
  pivot_wider(names_from = Indiv, values_from = count)

# note that for this count table, 0 is homozygous REF, 1 is heterozygote, 2 is homozygous ALT

## combine count tables
final_counts <- bind_cols(rbind(test_vcf@fix[, 1:2], test_vcf2@fix[, 1:2], non_vcf@fix[,1:2]),
                          bind_rows(counts2, counts3, counts4)) %>% 
  mutate(status = c(rep("cpg", nrow(counts2)), rep("polycpg", nrow(counts3)),rep("non", nrow(counts4)))) %>%
  mutate(POS = as.numeric(POS)) %>%
  select(-Key) %>%
  arrange(POS)

# sort gff - also filter based on gene
my_gff <- as_tibble(gff) %>%
  filter(V1 == my_chr & V3 == "gene") %>%
  select(start = V4, stop = V5, gene = V9) %>%
  arrange(start)

# get variant positions
variants <- final_counts %>% select(CHROM, POS)

# go through per row and extract snps
my_genes_list <- apply(my_gff, 1, function(x){
  z <- filter(variants, POS >= x[1] & POS <= x[2])
  z %>% mutate(gene = rep(x[3], nrow(z)))
})

# make a data.frame
my_genes <- bind_rows(my_genes_list) %>% select(-CHROM)

# merge with counts 
cpg_final <- left_join(final_counts, my_genes, by = "POS")

# remove any duplicate rows (created by genes with multiple hits)
cpg_final <- distinct(cpg_final, POS, .keep_all = T)

### correct for consecutive sites ###
# Add the filter column - ISSUE IS HERE
cpg_final$filter <- cpg$filter
# set index
idx <- which(cpg_final$filter == TRUE)
n_cons_sites <- length(idx) 
# how many consecutive sites?
sprintf("There are %s consecutive sites - %s percent of the total %s sites", n_cons_sites, signif(n_cons_sites/nrow(cpg_final)*100, 3), nrow(cpg_final))
# replacement counts
replace <- sapply(idx, function(y) {
  z <- cpg_final[c(y-1, y), ]
  z %>% summarise_all(min)
}, simplify = F) %>% bind_rows()
replace$filter <- rep(1, nrow(replace))

# remove original counts
cpg_final2 <- cpg_final[-(c(idx, idx-1)), ]
# append replacements
cpg_final <- bind_rows(cpg_final2, replace) %>% arrange(POS) 

sprintf("Corrected %s consecutive sites", n_cons_sites)

# write out
write_csv(cpg_final, out)

