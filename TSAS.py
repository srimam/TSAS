import os
import time
import logging
import sys
import argparse
from analysis.tools import Tools
from analysis.one_sample_analysis import OneSampleAnalysis
from analysis.two_sample_analysis import TwoSampleAnalysis

__version__ = '2.0.0'

LOG_FILE = 'TSAS.log'

def main(user_options, rundir):
    
    start = time.time()
    tools = Tools(LOG_FILE)
    logger = Tools.get_logger(LOG_FILE,__name__, logging.DEBUG)
    logger.info(f"Command: {' '.join(sys.argv)}")

    #Read in genome .fna file
    replicon_seqs, replicon_names = tools.parse_genome_file(user_options.genome_sequence)

    genome_size = 0
    replicon_info = []

    for id, seq in replicon_seqs.items():
        genome_size+=len(seq)
        replicon_info.append(f"{id}\t{len(seq)} bp\n")

    logger.info(f"Input Parameters\n\nAnalysis type : {user_options.analysis_type} sample analysis\nName of genome FASTA file : {user_options.genome_sequence}\nName of gff file : {user_options.gff}\nMapped reads file format : {user_options.mapped_read_format}\nName of mapped reads file (Treatment) : {user_options.treatment}\nName of mapped reads file (Control) : {user_options.control}\nThreshold for unique insertions : {user_options.min_hits}\n% of start and end to ignore : {user_options.clipping}\nCapping (2-sample only) : {user_options.capping}\nWeights (2-sample only) : {user_options.weight}\nTransposon sequence specificity: {user_options.target_sequence}\nResult format : {user_options.result_format}")
    logger.info(f"Number of replicons : {len(replicon_seqs.keys())}\n\nReplicon stats:\n{''.join(replicon_info)}\nGenome size : {genome_size} bp")

    #Read in GFF file
    logger.info("Parsing GFF and creating gene index...")
    coordinates, gene_index = tools.parse_gff(user_options.gff)
    if coordinates=={}:
        logger.error("Unable to parse GFF file to obtain coordinates. Please check the right file has been provided. Exiting...")
        sys.exit()
    else:
        logger.info("Finished parsing GFF file...")

    if user_options.analysis_type == '2':#ie two-sample analysis
        TwoSampleAnalysis(user_options.treatment.split(","), user_options.control.split(","), coordinates, gene_index, replicon_names, replicon_seqs, genome_size, user_options.min_hits, 
                        user_options.clipping, user_options.capping, user_options.weight, user_options.mapped_read_format, user_options.target_sequence, user_options.result_format, user_options.analysis_type, user_options.threads, LOG_FILE).run_analysis()
    else:#ie one-sample analysis
        OneSampleAnalysis(user_options.treatment.split(","), coordinates, gene_index, replicon_names, replicon_seqs, genome_size, user_options.min_hits, 
                        user_options.clipping, user_options.capping, user_options.mapped_read_format, user_options.target_sequence, user_options.result_format, user_options.analysis_type, user_options.threads, LOG_FILE).run_analysis()

    #logger.info(f"Run took {round((time.time()-start)/60,2)} mins to complete...")
    logger.info(f"Run took {time.time()-start} secs to complete...")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='TSAS codebase for TnSeq data analysis')
    parser.add_argument('-a', '--analysis_type', choices = ['1','2'], default = '1', help='Type of Tn-seq analysis to perform. Either 1 or 2 sample analysis (Default = 1)', required=True)
    parser.add_argument('-gs', '--genome_sequence', help='Genome sequence file i.e. .fna file', required=True)
    parser.add_argument('-gff', '--gff', help='GFF v3 file', required=True)
    parser.add_argument('-f', '--mapped_read_format', choices = ['Bowtie','SOAP', 'Eland'], default = 'Bowtie', help='Format of mapped reads file', required=True)
    parser.add_argument('-tr', '--treatment', help='Mapped read file(s) for treatment sample(s). Replicates provided as comma separated string with no spaces. If replicate samples are provided, analysis will average over the replicates. For 1-sample analysis, samples should be listed here alone', required=True)
    parser.add_argument('-c', '--control', help='Mapped read file(s) for control sample(s). Replicates provided as comma separated string with no spaces. This parameter is only used with two sample analysis and will be ignored by one-sample analysis. For a one-sample analysis, sample(s) should be list under Treatment. If replicate samples are provided, analysis will average over the replicates')
    parser.add_argument('-seq', '--target_sequence', default = '', help='Transposon sequence specificity (i.e., target sequence of transpson e.g., TA for Mariner transpons). Leave blank from random transposons like Tn5.')
    parser.add_argument('-mh', '--min_hits', default = 1, type = int, help='Threshold for minimum number of hits at a unique insertion site for it to be considered a true site. Default = 1')
    parser.add_argument('-cl', '--clipping', default = 0, type = float, help='Percentage of start and end of gene to be discarded when determining essentiality (e.g., a value of 5 equals 5 percent to be ignored @ start and @ end). Default = 0')
    parser.add_argument('-cap', '--capping', default = '1', choices = ['0','1','2','3'], help='Capping of reads at unique sites to minimize PCR/sequencing bias (0 - no capping; or 1 - capping with average hits per insertion site + 2 st. dev. (default); or 2 - capping with average hits per insertion site; or 3 - capping with median hits per insertion site). This parameter is only relevant to 2 sample analysis')
    parser.add_argument('-w', '--weight', default = '0', choices = ['0','1'], help='Weight reads per gene based on number of unique insertions in a gene (total hits*(unique hits per gene/average unique hits per gene)) 0 - No weights 1 - Weights used (- only relevant to 2 sample analysis)')
    parser.add_argument('-r', '--result_format', default = 'short', choices = ['long','short'], help='Result format. If Long it provides additional columns for capped and weighted reads, as well as unadjusted pvalues)')
    parser.add_argument('-t', "--threads", default = os.cpu_count(), type = int, help="Number of threads to use for analysis (Default = no of available cpus)", metavar="")

    user_options=parser.parse_args()
    main(user_options, os.getcwd())