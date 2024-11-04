## Introduction

This repository contains a set of bash commands and R scripts that are necessary to quantify the CpG positions across a population genomic dataset. The pipeline is relatively simple and we provide instructions here on how to use it.

The pipeline depends on a number of other tools including

- [bcftools](https://github.com/samtools/bcftools)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) 

To use everything here, you will need a reference genome, a set of genomic data - i.e. variant calls - from a number of individuals and ideally, information on the genome annotation, stored as a [gff](https://en.wikipedia.org/wiki/General_feature_format).

**NB: all examples here use variables as placeholders - i.e. $VCF - remember to switch them for your own paths/files**
 
### Step 1: Extract flanking sites

In order to identify CpG positions, we need flanking sites from the reference genome. First, we extract the positions of your variant calls:

```
bcftools query -f '%CHROM\t%POS\n' $VCF > init.bed
```

We next calculate the positions of the flanking sites - this is very easy with `awk`:

```
awk '{print $1"\t"$2-2"\t"$2+1}' init.bed > variant.bed
```

Finally we use bedtools to extract these regions from the reference genome:

```
bedtools getfasta -fi $REF -bed variant.bed -tab > cpg_flanking_sites.bed
```

We now have a file of flanking sites for the downstream analysis.

### Step 2: Count CpG positions

In order to count the CpG positions, we use the R script `cpg_counter_v0.2.R`. It can be used as a commandline tool like so:

```
Rscript cpg_counter_v0.2.R -v ${VCF} -f cpg_flanking_sites.bed -g ${GFF} -o $cpg.count.gz
```

Where the options are as follows:

- `-v` - path to vcf
- `-f` - path to flanking sites bedfile
- `-g` - path to gff 
- `-o` - output file path

The script is designed to run on a single chromosome or block of a chromosome. It cannot run across the entire genome at once. Therefore it is advised you split the analyses across the chromosomes or genome windows in order to speed it up.

### Step 3: Count CpG totals

If you want to calculate the total CpG site numbers, you can use the `cpg_totaler_v0.1.R` script. This is more lightweight than the CpG counter facility and is designed to run on very large datasets (i.e. all variant and invariant sites in a genome). It's usuage is almost identical to the counting script.

```
Rscript cpg_totaler_v0.1.R -v ${VCF} -f ${CHR}_cpg_flanking_sites.bed -g ${GFF} -o ${CHR}_cpg.count.gz
```

Where the options are as follows:

- `-v` - path to vcf
- `-f` - path to flanking sites bedfile
- `-g` - path to gff 
- `-o` - output file path

This script should also be run across chromosomes.

### Output explained