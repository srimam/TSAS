# Tn-seq analysis software (TSAS 2.0)

Transposon mutagenesis followed by high-throughput sequencing (Tn-seq) is a robust approach for genome-wide identification of putative essential genes in any organism of interest for which the appropriate genetic tools are available (or can be developed). However, easily accessible computational tools for the streamlined analysis of the data from Tn-seq experiments are currently limited.
This document provides instructions for how to use the Tn-seq analysis software (TSAS) described in Burger et al. 2017, which provides tools for various kinds of analysis of Tn-seq data sets. 

TSAS 2.0 is a rewritten of TSAS in Python to simplify use, while also significantly speed up the analysis.

## Required Input Data
The following input data files are required to run TSAS:
a.	Aligned reads file obtained from mapping short reads to a reference genome. TSAS currently handles Bowtie (version 1 and 2), SOAP and Eland result formats. This could be extended in future updates.
b.	Genome sequence file in Fasta format. The name(s) of the chromosome(s) in this file must match that used for mapping aligned reads to the reference genome.
c.	GFF v3 file containing the coordinates of genomic elements in the reference genome. Again, the chromosome names in this file must match those in the genome sequence and aligned reads files.

## Usage

### Usage options

```python
usage: TSAS.py [-h] -a {1,2} -gs GENOME_SEQUENCE -gff GFF -f
               {Bowtie,SOAP,Eland} -tr TREATMENT [-c CONTROL]
               [-seq TARGET_SEQUENCE] [-mh MIN_HITS] [-cl CLIPPING]
               [-cap {0,1,2,3}] [-w {0,1}] [-r {long,short}] [-t]

TSAS codebase for TnSeq data analysis

optional arguments:
  -h, --help            show this help message and exit
  -a {1,2}, --analysis_type {1,2}
                        Type of Tn-seq analysis to perform. Either 1 or 2
                        sample analysis (Default = 1)
  -gs GENOME_SEQUENCE, --genome_sequence GENOME_SEQUENCE
                        Genome sequence file i.e. .fna file
  -gff GFF, --gff GFF   GFF v3 file
  -f {Bowtie,SOAP,Eland}, --mapped_read_format {Bowtie,SOAP,Eland}
                        Format of mapped reads file
  -tr TREATMENT, --treatment TREATMENT
                        Mapped read file(s) for treatment sample(s).
                        Replicates provided as comma separated string with no
                        spaces. If replicate samples are provided, analysis
                        will average over the replicates. For 1-sample
                        analysis, samples should be listed here alone
  -c CONTROL, --control CONTROL
                        Mapped read file(s) for control sample(s). Replicates
                        provided as comma separated string with no spaces.
                        This parameter is only used with two sample analysis
                        and will be ignored by one-sample analysis. For a one-
                        sample analysis, sample(s) should be list under
                        Treatment. If replicate samples are provided, analysis
                        will average over the replicates
  -seq TARGET_SEQUENCE, --target_sequence TARGET_SEQUENCE
                        Transposon sequence specificity (i.e., target sequence
                        of transpson e.g., TA for Mariner transpons). Leave
                        blank from random transposons like Tn5.
  -mh MIN_HITS, --min_hits MIN_HITS
                        Threshold for minimum number of hits at a unique
                        insertion site for it to be considered a true site.
                        Default = 1
  -cl CLIPPING, --clipping CLIPPING
                        Percentage of start and end of gene to be discarded
                        when determining essentiality (e.g., a value of 5
                        equals 5 percent to be ignored @ start and @ end).
                        Default = 0
  -cap {0,1,2,3}, --capping {0,1,2,3}
                        Capping of reads at unique sites to minimize
                        PCR/sequencing bias (0 - no capping; or 1 - capping
                        with average hits per insertion site + 2 st. dev.
                        (default); or 2 - capping with average hits per
                        insertion site; or 3 - capping with median hits per
                        insertion site). This parameter is only relevant to 2
                        sample analysis
  -w {0,1}, --weight {0,1}
                        Weight reads per gene based on number of unique
                        insertions in a gene (total hits*(unique hits per
                        gene/average unique hits per gene)) 0 - No weights 1 -
                        Weights used (- only relevant to 2 sample analysis)
  -r {long,short}, --result_format {long,short}
                        Result format. If Long it provides additional columns
                        for capped and weighted reads, as well as unadjusted
                        pvalues)
  -t , --threads        Number of threads to use for analysis (Default = no of
                        available cpus)
```

### Example runs

### One sample analysis 
```python
python3 TSAS.py -a 1 \
-gs example_data/NC_002516.fna \
-gff example_data/NC_002516.gff \
-f Bowtie \
-tr example_data/sample1BC2.fastq_mapped_M1,example_data/sample2BC2.fastq_mapped_M1 \
-mh 1 \
-cl 10 \
-seq TA
```

### Two sample analysis
```python
python3 TSAS.py -a 2 \
-gs example_data/NC_002516.fna \
-gff example_data/NC_002516.gff \
-f Bowtie \
-c example_data/sample1BC2.fastq_mapped_M1,example_data/sample2BC2.fastq_mapped_M1 \
-tr example_data/sample3BC2.fastq_mapped_M1,example_data/sample4BC2.fastq_mapped_M1 \
-mh 1 \
-cl 10 \
-cap 3 \
-w 1
```

Output files formats are identical to TSAS 1. See TSAS 1 user guide for description of results columns.

Example data folder will need to be uncompressed after cloning repo, before running test commands.

## Dependencies

Python >=3.7

Numpy >=1.21

Scipy >=1.7


TSAS will probably work with older versions of Numpy and Scipy but hasn't been tested with these.

