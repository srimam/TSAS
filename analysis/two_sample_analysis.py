import logging
import numpy as np
from .tools import Tools

LOG_FILE = 'TSAS.log'
tools = Tools(LOG_FILE)

class TwoSampleAnalysis:
    '''
        Class to conduct TSAS two sample analysis
    '''

    def __init__(self, mapped_reads_treatment: list, mapped_reads_control: list, coordinates: dict, gene_index: dict, replicon_names: dict, replicon_seqs: dict, genome_size: int, 
                min_hits: int, clipping: float, capping: str, weight: str, mapped_read_format: str, target_seq: str, result_format: str, analysis_type: str, threads: int, filename: str) -> None:
        
        self.logger = Tools.get_logger(filename,__name__, logging.DEBUG)
        self.mapped_reads_treatment = mapped_reads_treatment
        self.mapped_reads_control = mapped_reads_control
        self.coordinates = coordinates
        self.gene_index = gene_index
        self.replicon_names = replicon_names
        self.replicon_seqs = replicon_seqs
        self.genome_size = genome_size
        self.min_hits = min_hits
        self.clipping = clipping
        self.capping = capping
        self.weight = weight
        self.mapped_read_format = mapped_read_format
        self.target_seq = target_seq
        self.result_format = result_format
        self.analysis_type = analysis_type
        self.threads = threads

    def run_analysis(self) -> None:
        #Process mapped reads file(s)
        mapped_reads_c = tools.mp_process_mapped_reads(self.mapped_reads_control, self.replicon_names, self.min_hits, self.mapped_read_format, self.logger, self.threads)   
        mapped_reads_t = tools.mp_process_mapped_reads(self.mapped_reads_treatment, self.replicon_names, self.min_hits, self.mapped_read_format, self.logger, self.threads)   
        
        #Summarize mapped reads to gene level data
        experiment_stats_c = tools.gene_level_analysis(self.mapped_reads_control, mapped_reads_c, self.coordinates, self.gene_index, self.target_seq, self.replicon_seqs, self.min_hits, self.clipping, self.capping, self.analysis_type, self.logger)
        experiment_stats_t = tools.gene_level_analysis(self.mapped_reads_treatment, mapped_reads_t, self.coordinates, self.gene_index, self.target_seq, self.replicon_seqs, self.min_hits, self.clipping, self.capping, self.analysis_type, self.logger)
        
        #combine data
        mapped_reads = {}
        mapped_reads.update(mapped_reads_c)
        mapped_reads.update(mapped_reads_t)

        experiment_stats = {}
        experiment_stats.update(experiment_stats_c)
        experiment_stats.update(experiment_stats_t)

        total_reads_per_gene_weighted = {}
        if self.weight=='1':
            total_reads_per_gene_weighted = tools.weight_reads(experiment_stats, mapped_reads, self.logger)
        else:#if no weight just use the precomputed total_reads_per_gene
            for reads_file in mapped_reads:
                reads_per_gene = {}
                for gid in experiment_stats[reads_file][0].keys():#gene level stats keys (gids)
                    reads_per_gene[gid] = experiment_stats[reads_file][0][gid][1]#Gene level stats, gid,total_reads_per_gene
                total_reads_per_gene_weighted[reads_file]=reads_per_gene

        self.logger.info("Calculating gene level fold enrichment and significance between treatment and control")

        #compute summary stats for combined replicates
        summary_stats_c = tools.average_replicates2(mapped_reads_c, experiment_stats_c, total_reads_per_gene_weighted)
        summary_stats_t = tools.average_replicates2(mapped_reads_t, experiment_stats_t, total_reads_per_gene_weighted)

        #compute ratios
        ratios = tools.compute_ratios(summary_stats_c, summary_stats_t)

        #compuate tstat pvalues
        tstat_pvals = {}
        tstat_adj_pvals = {}
        tstat_pvals_reads = {}
        tstat_adj_pvals_reads = {}
        if len(mapped_reads_c) >=2 and len(mapped_reads_t) >=2: #if there are at least 2 control and 2 test samples compute tstat
            tstat_pvals, tstat_adj_pvals = tools.compute_tsat(len(mapped_reads_c), len(mapped_reads_t), summary_stats_c[5], summary_stats_t[5], summary_stats_c[4], summary_stats_t[4], ratios[2])
            tstat_pvals_reads, tstat_adj_pvals_reads = tools.compute_tsat(len(mapped_reads_c), len(mapped_reads_t), summary_stats_c[1], summary_stats_t[1], summary_stats_c[0], summary_stats_t[0], ratios[2])


        #compute proportions pvalues
        prop_pvalues_ins, prop_adj_pvalues_ins = tools.compute_proportions(summary_stats_t[4],summary_stats_c[4],summary_stats_t[6],summary_stats_c[6], ratios[2])
        prop_pvalues_reads, prop_adj_pvalues_reads = tools.compute_proportions(summary_stats_t[0],summary_stats_c[0],summary_stats_t[8],summary_stats_c[8], ratios[2])

        #compute fisher/chi2 pvalues
        fisher_pvalues_ins, fisher_adj_pvalues_ins = tools.compute_fisher(summary_stats_t[4],summary_stats_c[4],summary_stats_t[6],summary_stats_c[6], ratios[2])
        fisher_pvalues_reads, fisher_adj_pvalues_reads = tools.compute_fisher(summary_stats_t[0],summary_stats_c[0],summary_stats_t[8],summary_stats_c[8], ratios[2])
        

        #Output results
        self.logger.info("Printing out results...")
        with open("Conditional_essentiality.txt", "w") as outfile:
            if self.result_format == "long":
                if len(tstat_pvals) > 0:
                    outfile.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Raw Reads(treatment)\tAve. Raw Read (control)\tAve. Capped_reads(treatment)\tAve. Capped reads (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (Reads)\tpvalue (proportions_insertions)\tAdj. pvalue (proportions_insertions)\tpvalue (proportions_reads)\tAdj. pvalue (proportions_reads)\tpvalue (Fisher_insertions)\tAdj. pvalue (Fisher_insertions)\tpvalue (Fisher_reads)\tAdj. pvalue (Fisher_reads)\tpvalue (t-test insertions)\tAdj. pvalue (t-test insertions)\tpvalue (t-test reads)\tAdj. pvalue (t-test reads)\n")
                else:
                    outfile.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Raw Reads(treatment)\tAve. Raw Read (control)\tAve. Capped_reads(treatment)\tAve. Capped reads (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (Reads)\tpvalue (proportions_insertions)\tAdj. pvalue (proportions_insertions)\tpvalue (proportions_reads)\tAdj. pvalue (proportions_reads)\tpvalue (Fisher_insertions)\tAdj. pvalue (Fisher_insertions)\tpvalue (Fisher_reads)\tAdj. pvalue (Fisher_reads)\n")
                for gid, pval in dict(sorted(prop_adj_pvalues_reads.items(), key=lambda x:(x[1]))).items():
                    annotation = ""
                    if len(self.coordinates[gid]) > 4:
                        annotation = self.coordinates[gid][4]
                    if len(tstat_pvals) > 0:
                        outfile.write(f"{gid}\t{annotation}\t{summary_stats_t[4][gid]}\t{summary_stats_c[4][gid]}\t{summary_stats_t[3][gid]}\t{summary_stats_c[3][gid]}\t{summary_stats_t[2][gid]}\t{summary_stats_c[2][gid]}\t{summary_stats_t[0][gid]}\t{summary_stats_c[0][gid]}\t{ratios[1][gid]}\t{np.log2(ratios[1][gid])}\t{ratios[0][gid]}\t{np.log2(ratios[0][gid])}\t{prop_pvalues_ins[gid]}\t{prop_adj_pvalues_ins[gid]}\t{prop_pvalues_reads[gid]}\t{prop_adj_pvalues_reads[gid]}\t{fisher_pvalues_ins[gid]}\t{fisher_adj_pvalues_ins[gid]}\t{fisher_pvalues_reads[gid]}\t{fisher_adj_pvalues_reads[gid]}\t{tstat_pvals[gid]}\t{tstat_adj_pvals[gid]}\t{tstat_pvals_reads[gid]}\t{tstat_adj_pvals_reads[gid]}\n")
                    else:
                        outfile.write(f"{gid}\t{annotation}\t{summary_stats_t[4][gid]}\t{summary_stats_c[4][gid]}\t{summary_stats_t[3][gid]}\t{summary_stats_c[3][gid]}\t{summary_stats_t[2][gid]}\t{summary_stats_c[2][gid]}\t{summary_stats_t[0][gid]}\t{summary_stats_c[0][gid]}\t{ratios[1][gid]}\t{np.log2(ratios[1][gid])}\t{ratios[0][gid]}\t{np.log2(ratios[0][gid])}\t{prop_pvalues_ins[gid]}\t{prop_adj_pvalues_ins[gid]}\t{prop_pvalues_reads[gid]}\t{prop_adj_pvalues_reads[gid]}\t{fisher_pvalues_ins[gid]}\t{fisher_adj_pvalues_ins[gid]}\t{fisher_pvalues_reads[gid]}\t{fisher_adj_pvalues_reads[gid]}\n")

            else:
                if len(tstat_pvals) > 0:
                    outfile.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (reads)\tAdj. pvalue (proportions_insertions)\tAdj. pvalue (proportions_reads)\tAdj. pvalue (Fisher_insertions)\tAdj. pvalue (Fisher_reads)\tAdj. pvalue (t-test insertions)\tAdj. pvalue (t-test reads)\n")
                else:
                    outfile.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (reads)\tAdj. pvalue (proportions_insertions)\tAdj. pvalue (proportions_reads)\tAdj. pvalue (Fisher_insertions)\tAdj. pvalue (Fisher_reads)\n")
                for gid, pval in dict(sorted(prop_adj_pvalues_reads.items(), key=lambda x:(x[1]))).items():
                    annotation = ""
                    if len(self.coordinates[gid]) > 4:
                        annotation = self.coordinates[gid][4]
                    if len(tstat_pvals) > 0:
                        outfile.write(f"{gid}\t{annotation}\t{summary_stats_t[4][gid]}\t{summary_stats_c[4][gid]}\t{summary_stats_t[0][gid]}\t{summary_stats_c[0][gid]}\t{ratios[1][gid]}\t{np.log2(ratios[1][gid])}\t{ratios[0][gid]}\t{np.log2(ratios[0][gid])}\t{prop_adj_pvalues_ins[gid]}\t{prop_adj_pvalues_reads[gid]}\t{fisher_adj_pvalues_ins[gid]}\t{fisher_adj_pvalues_reads[gid]}\t{tstat_adj_pvals[gid]}\t{tstat_adj_pvals_reads[gid]}\n")
                    else:
                        outfile.write(f"{gid}\t{annotation}\t{summary_stats_t[4][gid]}\t{summary_stats_c[4][gid]}\t{summary_stats_t[0][gid]}\t{summary_stats_c[0][gid]}\t{ratios[1][gid]}\t{np.log2(ratios[1][gid])}\t{ratios[0][gid]}\t{np.log2(ratios[0][gid])}\t{prop_adj_pvalues_ins[gid]}\t{prop_adj_pvalues_reads[gid]}\t{fisher_adj_pvalues_ins[gid]}\t{fisher_adj_pvalues_reads[gid]}\n")