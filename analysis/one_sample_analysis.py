import logging
from .tools import Tools
from scipy import stats


LOG_FILE = 'TSAS.log'
tools = Tools(LOG_FILE)

class OneSampleAnalysis:
    '''
        Class to conduct TSAS one sample analysis
    '''

    def __init__(self, mapped_reads_exp: list, coordinates: dict, gene_index: dict, replicon_names: dict, replicon_seqs: dict, genome_size: int, 
                min_hits: int, clipping: float, capping: str, mapped_read_format: str, target_seq: str, result_format: str, analysis_type: str, threads: int, filename: str) -> None:
        
        self.logger = Tools.get_logger(filename,__name__, logging.DEBUG)
        self.mapped_reads_exp = mapped_reads_exp
        self.coordinates = coordinates
        self.gene_index = gene_index
        self.replicon_names = replicon_names
        self.replicon_seqs = replicon_seqs
        self.genome_size = genome_size
        self.min_hits = min_hits
        self.clipping = clipping
        self.capping = capping
        self.mapped_read_format = mapped_read_format
        self.target_seq = target_seq
        self.result_format = result_format
        self.analysis_type = analysis_type
        self.threads = threads

    def run_analysis(self) -> None:
        
        genome_wide_insertion_sites = 0
        if(self.target_seq != ""):
            for id, seq in self.replicon_seqs.items():
                genome_wide_insertion_sites+=tools.count_insertion_sites(seq,self.target_seq)

            if genome_wide_insertion_sites == 0:
                self.logger.info("Number of insertion sites for provided sequence couldn't be determined... Assuming random insertion sequence!")   
            else:
                self.logger.info(f"Total potential {self.target_seq} insertion sites in genome: {genome_wide_insertion_sites}")

        #Process mapped reads file(s)
        mapped_reads = tools.mp_process_mapped_reads(self.mapped_reads_exp, self.replicon_names, self.min_hits, self.mapped_read_format, self.logger, self.threads)   
        
        #Summarize mapped reads to gene level data
        experiment_stats = tools.gene_level_analysis(self.mapped_reads_exp, mapped_reads, self.coordinates, self.gene_index, self.target_seq, self.replicon_seqs, self.min_hits, self.clipping, self.capping, self.analysis_type, self.logger)
        
        #Average replicates
        mean_stats = tools.average_replicates(self.mapped_reads_exp, experiment_stats)

        #Compute significance using binomial distribution
        pvalues = {}
        corrected_pvalues = {}
        bonferroni_pvalues = {}
        corrected_rev_pvalues = {}
        bonferroni_rev_pvalues = {}

        rev_pvalues = {}
        for gid in self.coordinates.keys():
            pvalue = 1
            if genome_wide_insertion_sites == 0:#if no sequence specific transposon insertion site or it could not be evalulated for some reason
                pvalue = stats.binom.sf(round(mean_stats[1][gid],0), #mean_unique_hits_per_gene (k)
                experiment_stats[self.mapped_reads_exp[0]][0][gid][4], #available_gene_insertion_sites (n)
                mean_stats[0]/self.genome_size)#mean_unique_hits_genome/genome_wide_insertion_sites (p)
            else:#for sequence specific insertions
                if mean_stats[1][gid] > experiment_stats[self.mapped_reads_exp[0]][0][gid][4]:#ie more insertion sites observed than expected based on provided sequence specificity
                    pvalue = stats.binom.sf(experiment_stats[self.mapped_reads_exp[0]][0][gid][4], #mean_unique_hits_per_gene (k)
                    experiment_stats[self.mapped_reads_exp[0]][0][gid][4], #available_gene_insertion_sites (n)
                    mean_stats[0]/genome_wide_insertion_sites)#mean_unique_hits_genome/genome_wide_insertion_sites (p)
                    self.logger.warning(f"More insertion sites observed than expected for {gid} based on provided specificity. Capping at max expected. Please verify target sequence specified...")
                else:
                    pvalue = stats.binom.sf(round(mean_stats[1][gid],0), #mean_unique_hits_per_gene (k)
                    experiment_stats[self.mapped_reads_exp[0]][0][gid][4], #available_gene_insertion_sites (n)
                    mean_stats[0]/genome_wide_insertion_sites)#mean_unique_hits_genome/genome_wide_insertion_sites (p)

            rev_pvalues[gid] = pvalue
            pvalues[gid] = 1 - pvalue

        #BH and Bonferroni correction
        corrected_pvalues, bonferroni_pvalues = Tools.correct_pvalues(pvalues, 0)
        corrected_rev_pvalues, bonferroni_rev_pvalues = Tools.correct_pvalues(rev_pvalues, 0)

        #Output results
        self.logger.info("Printing out results...")
        with open("Essential_genes.txt", "w") as outfile:
            if self.result_format == "long":
                outfile.write("Gene ID\tAnnotation\tGene length (bp)\tNo. of Unique hits\tNormalize Unique hits (hits/bp)\tMean Total number of reads\tPvalue (Essential)\tAdj. Pvalue (Essential)\tFWER (Essential)\tPvalue (Improved fitness)\tAdj. Pvalue (Improved fitness)\tFWER (Improved fitness)\n")
                for gid, pval in dict(sorted(corrected_pvalues.items(), key=lambda x:(x[1]))).items():
                    annotation = ""
                    if len(self.coordinates[gid]) > 4:
                        annotation = self.coordinates[gid][4]
                    outfile.write(f"{gid}\t{annotation}\t{self.coordinates[gid][2]-self.coordinates[gid][1]+1}\t{mean_stats[1][gid]}\t{mean_stats[1][gid]/experiment_stats[self.mapped_reads_exp[0]][0][gid][2]}\t{mean_stats[2][gid]}\t{pvalues[gid]}\t{corrected_pvalues[gid]}\t{bonferroni_pvalues[gid]}\t{rev_pvalues[gid]}\t{corrected_rev_pvalues[gid]}\t{bonferroni_rev_pvalues[gid]}\n")
            else:
                outfile.write("Gene ID\tAnnotation\tGene length (bp)\tNo. of Unique hits\tNormalize Unique hits (hits/bp)\tMean Total number of reads\tAdj. Pvalue (Essential)\tFWER (Essential)\tAdj. Pvalue (Improved fitness)\tFWER (Improved fitness)\n")
                for gid, pval in dict(sorted(corrected_pvalues.items(), key=lambda x:(x[1]))).items():
                    annotation = ""
                    if len(self.coordinates[gid]) > 4:
                        annotation = self.coordinates[gid][4]
                    outfile.write(f"{gid}\t{annotation}\t{self.coordinates[gid][2]-self.coordinates[gid][1]+1}\t{mean_stats[1][gid]}\t{mean_stats[1][gid]/experiment_stats[self.mapped_reads_exp[0]][0][gid][2]}\t{mean_stats[2][gid]}\t{corrected_pvalues[gid]}\t{bonferroni_pvalues[gid]}\t{corrected_rev_pvalues[gid]}\t{bonferroni_rev_pvalues[gid]}\n")
