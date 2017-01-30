Tn-seq analysis software (TSAS)

Transposon mutagenesis followed by high-throughput sequencing (Tn-seq) is a robust approach for genome-wide identification of putative essential genes in any organism of interest for which the appropriate genetic tools are available (or can be developed). However, easily accessible computational tools for the streamlined analysis of the data from Tn-seq experiments are currently limited.
This document provides detailed instructions for how to use the Tn-seq analysis software (TSAS) described in Burger et al. 2017, which provides tools for various kinds of analysis of Tn-seq data sets. See TSAS User Guide.docx for detailed description of usage and results

System requirements
The only system requirement for this tool is to have Java runtime installed and configured on your machine. Due to occasional backward compatibility issues between different versions of Java, it may be beneficial to also have Java software development kit (SDK) installed and configured to enable recompilation of the code if necessary. 
Since the TSAS code is all written in Java it should be platform independent and run on any machine with Java installed. It has been tested in Linux (CentOS 6), Windows, and Mac environments with Java 1.8.

Required Data
The following input data files are required to run TSAS:
a.	Aligned reads file obtained from mapping short reads to a reference genome. TSAS currently handles Bowtie (version 1 and 2), SOAP and Eland result formats. This could be extended in future updates.
b.	Genome sequence file in Fasta format. The name(s) of the chromosome(s) in this file must match that used for mapping aligned reads to the reference genome.
c.	GFF v3 file containing the coordinates of genomic elements in the reference genome. Again, the chromosome names in this file must match those in the genome sequence and aligned reads files.
