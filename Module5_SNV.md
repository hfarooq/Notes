---
layout: tutorial_page
permalink: /BiCG_2017
title: bioinformaticsdotca
header1: Bioinformatics for Cancer Genomics 2017
header2: Lab Module 5 - Somatic Mutations and Annotation
image: /site_images/CBW-CSHL-graphic-square.png
home:
https://bioinformaticsdotca/BiCG_2017/tree/c1a499370d05dccf91ab61207b711295e0f31c08
---

# Lab Module 5 - Somatic Mutations and Annotation

This lab is designed to provide an introduction into Somatic Nucleotide Variation detection using two common programs: [Strelka](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts271) and [MuTect](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html)

## Setup

First login into the server, and then enter the workspace directory:

```
cd ~/workspace
```

In order to keep our files in one location, we're going to create a new directory for this module and enter it:

```
mkdir Module5
cd Module5
```

Now to ease the use of commands that will follow, we're going to set some environment variables. These variables will stand in for the directories where our softwares of interest are location.

>**Environment variables are only temporary. This means that once you log out of the server, these variables will be erased and must be reset if you'd like to use them again.**

Now we'll set out environment variables for MuTect and SnpEff (which will be used for processing the outputs from the callers):
```
MUTECT_DIR=/usr/local/mutect
SNPEFF_DIR=/usr/local/snpEff
```

## Linking the Sequencing and Referencing Data

For this lab module, we'll be using exome data on the HCC1395 breast cancer cell line. The tumour and normal bams have already been processed and placed on the server. So we'll create a soft link to the directory that it's stored:
```
ln -s /home/ubuntu/CourseData/CG_data/Module5/HCC1395
```
For mutation calling we also need a reference genome, and so we're going to make a link to that folder as well:
```
ln -s /home/ubuntu/CourseData/CG_data/Module5/ref_data
```
For this lab we're going to limit our analysis to just the 7MB and 8MB region of chromosome 17 to ensure processing occurs quickly. The files we'll be focusing on can be viewed using the following command:
```
ls HCC1395/HCC1395_exome*ordered.17*
```
You should see the files:
* HCC1395_exome_normal_ordered.17.7MB-8MB.bam
* HCC1395_exome_normal_ordered.17.7MB-8MB.bam.bai
* HCC1395_exome_tumour_ordered.17.7MB-8MB.bam
* HCC1395_exome_tumour_ordered.17.7MB-8MB.bam.bai

Let's take a bit of time to look into the bam files we have, starting with the header. The header contains information about the reference genome, commands used to generate the file, and on read groups (information about sample/flow cell/lane etc.)
```
samtools view -H HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam | less -S
```

The main contents of the file contain the one aligned read (read name, chromosome, position etc.) per line:
```
samtools view HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam | less -S
```

In order to determine statistics about the bam file, we can use samtools' built in tool:
```
samtools flagstat HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam
```

# Predicting SNVs

## MuTect

* Insert Mutect Introduction *

In order to run MuTect, we simply run the java application of the program in one go. For MuTect to run, we need to have common SNV's present from dbsnp and cosmic to help filter out noise and focus on damaging mutations. These files are stored in the same path as our MUTECT_DIR variable.
Looking at the dbsnp file:
```
less $MUTECT_DIR/dbsnp_132_b37.leftAligned.vcf
```
And now the cosmic file:
```
less $MUTECT_DIR/b37_cosmic_v54_120711.vcf
```

In order to run MuTect, we'll be running the following command:
```
mkdir results; mkdir results/mutect;
/usr/local/java6/bin/java -Xmx4g -jar $MUTECT_DIR/muTect-1.1.4.jar \
--analysis_type MuTect --reference_sequence ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.fa \
--dbsnp $MUTECT_DIR/dbsnp_132_b37.leftAligned.vcf --cosmic $MUTECT_DIR/b37_cosmic_v54_120711.vcf \
--intervals 17:7000000-8000000 --input_file:normal HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam \
--input_file:tumor HCC1395/HCC1395_exome_tumour_ordered.17.7MB-8MB.bam --vcf results/mutect/mutect.call_stats.vcf
```

The `-Xmx4g` option allocates 4 GB of RAM for MuTect to run it's analysis. MuTect will output a vcf file with the calls.

The calls file contains information about the each SNV including information such as the reference allele, alternate allele and whether the SNV should be kept based on MuTect's quality filtering.

```
less -S results/mutect.call_stats.vcf
```


## Strelka

We will first call mutations using Strelka. Create a local copy of the Strelka config file.  Strelka provides aligner specific config files for bwa, eland, and isaac.  Each file contains default configuration parameters that work well with the aligner.  The bam files we are working with were created using bwa, so we select that config file and make a local copy to make changes.

```
mkdir config
cp /usr/local/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini
```

Since our data is exome and so the coverage of the file is different, we need to change the `isSkipDepthFilters` parameter in the strelka_config_bwa.ini file. We'll first create a new config file for exome analysis:

```
cp config/strelka_config_bwa.ini config/strelka_config_bwa_exome.ini
```

Now let's edit the `config/strelka_config_bwa_exome.ini` and change the `isSkipDepthFilters = 0` to `isSkipDepthFilters = 1`.  We will use this using the vim editor:

```
vim config/strelka_config_bwa_exome.ini
```

In order to quit the editor while saving changes, press **ESC**, following by **:x!** and pressing enter. The reason why we do this is described on the [Strelka FAQ page](https://sites.google.com/site/strelkasomaticvariantcaller/home/faq):

> The depth filter is designed to filter out all variants which are called above a multiple of the mean chromosome depth, the default configuration is set to filter variants with a depth greater than 3x the chromosomal mean. If you are using exome/targeted sequencing data, the depth filter should be turned off...
> 
> However in whole exome sequencing data, the mean chromosomal depth will be extremely small, resulting in nearly all variants being (improperly) filtered out.

If you were doing this for whole genome sequencing data, then you should leave this parameter set to 0 as the depth of coverage won't be as high. 

A Strelka analysis is performed in 2 steps.  In the first step we provide Strelka with all the information it requires to run the analysis, including the tumour and normal bam filenames, the config, and the reference genome.  Strelka will create an output directory with the setup required to run the analysis.

```
configureStrelkaWorkflow.pl \
    --tumor HCC1395/HCC1395_exome_tumour_ordered.17.7MB-8MB.bam \
    --normal HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam \
    --ref ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.fa \
    --config config/strelka_config_bwa_exome.ini \
    --output-dir results/strelka/
```

The output directory will contain a _makefile_ that can be used with the tool _make_.  One benefit of the makefile style workflow is that it can be easily parallelized using _qmake_ on a grid engine cluster.  

To run the Strelka analysis, use make and specify the directory constructed by `configureStrelkaWorkflow.pl` with make's `'-C'` option.

```
make -C results/strelka/ -j 2
```

The `-j 2` parameter specifies that we want to use 2 cores to run Strelka. Change this number to increase or decrease the parallelization of the job. The more cores the faster the job will be, but the higher the load on the machine that is running Strelka. 

> If you have access to a grid engine cluster, you can replace the command `make` with `qmake` to launch Strelka on the cluster.

Strelka has the benefit of calling SNVs and small indels.  Additionally, Strelka calculates variant quality and filters the data in two tiers.  The filenames starting with `passed` contain high quality candidates, and filenames starting with `all` contain high quality and marginal quality candidates.

The Strelka results are in VCF format, with additional fields explained on the [strelka website](https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output).

```
less -S results/strelka/results/passed.somatic.snvs.vcf
less -S results/strelka/results/passed.somatic.indels.vcf
```

## Converting the VCF format into a tabular format

The VCF format is sometimes not useful for visualization and data exploration purposes which often requires the data to be in tabular format. We can convert from VCF format to tabular format using the extractField() function from SnpSift/SnpEff. Since each mutation caller has a different set of output values in the VCF file, the command needs be adjusted for the mutation caller. 

For example, to convert the Strelka VCF file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "."  results/strelka/results/passed.somatic.snvs.vcf CHROM POS ID REF ALT QUAL FILTER QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > results/strelka/results/passed.somatic.snvs.txt 
```

To convert the MuTect VCF file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "." results/mutect/mutect.call_stats.vcf CHROM POS ID REF ALT QUAL FILTER  > results/mutect/mutect.call_stats.txt
```

The -e parameter specifies how to represent empty fields. In this case, the "." character is placed for any empty fields. This facilitates loading and completeness of data. For more details on the extractField() function see the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#Extract).

## Data Exploration

### Integrative Genomics Viewer (IGV) 

A common step after prediction of SNVs is to visualize these mutations in IGV. Let's load these bam into IGV. Open IGV, then:

1. Change the genome to hg19 (if it isn't already)
2. File -> Load from URL ...
    * http://cbwxx.dyndns.info//Module5/HCC1395/HCC1395_exome_tumour_ordered.17.7MB-8MB.bam
    * http://cbwxx.dyndns.info//Module5/HCC1395/HCC1395_exome_normal_ordered.17.7MB-8MB.bam

Where the xx is your student number. Once the tumour and normal bam have been loaded into IGV, we can investigate a few predicted positions in IGV:

* 17:7491818
* 17:7578406
* 17:7482929
* 17:7710987

Manually inspecting these predicted SNVs in IGV is a good way to verify the predictions and also identify potential false positives:

> When possible, you should always inspect SNVs

### Exploration in R

While IGV is good for visualizing individual mutations, looking at more global characteristics would require loading the data into an analysis language like R:

We will use exome-wide SNV predictions for Strelka for these analyses. These processed tabular text files along with the `analyzeSNVResults.Rmd` RMarkdown file that contains the R code for the analysis can downloaded as a package. 

Open the `analyzeSNVResults.Rmd` in RStudio now. 
