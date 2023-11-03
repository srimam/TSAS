import logging
import sys
import re
import warnings
import numpy as np
import multiprocessing as mp
from scipy import stats
from itertools import repeat

class Tools:

    def __init__(self, filename) -> None:
        self.logger = Tools.get_logger(filename,__name__, logging.DEBUG)

    @staticmethod
    def get_logger(filename: str, name: str, level=logging.DEBUG):
        """
        Method to create logger for each run
        Input:
            filename: name of output logfile
            name: location where logging output occurs e.g., class name
            level: level at which logging occurs
        Output:
            logger
        """
        formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s\n%(message)s\n', datefmt='%Y-%m-%d %H:%M:%S')
        sh = logging.StreamHandler(sys.stdout)
        fh = logging.FileHandler(filename)
        sh.setFormatter(formatter)
        fh.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.addHandler(sh)
        logger.addHandler(fh)
        logger.setLevel(level)

        return logger


    def parse_genome_file(self, input_file: str) -> list:
        """
        Method to parse genome fna file
        Input:
            input_file: genome fna file
        Output:
            list of dictionaries holding replicon/chromosome sequences and replicon names
        """

        #read in genome sequence with associated chromosome IDs into dictionary
        replicon_seqs = {}
        replicon_names={}
        chr_ID = ""
        seq = []
        count=0
        
        with open(input_file) as infile:
            for line in infile:
                if line.startswith(">") and count==0:
                    chr_ID = line.replace(">","").strip().split(" ")[0]
                    replicon_names[chr_ID] = line.replace(">","").strip()
                    count=1
                elif line.startswith(">") and count==1:
                    replicon_seqs[chr_ID]="".join(seq)
                    chr_ID = line.replace(">","").strip().split(" ")[0]
                    replicon_names[chr_ID] = line.replace(">","").strip()
                    seq = []
                else:
                    seq.append(line.strip())
            
        replicon_seqs[chr_ID]="".join(seq)

        return [replicon_seqs, replicon_names]


    def run_sanity_checks(self):
        pass


    def parse_gff(self, gff_file: str) -> list:
        """Method to read in and parse gff files returning a dict with coordinates {chr, start, end, orientation} and create a gene index mapping each nucleotide to a gene
            Input:
                gff_file: gffv3 file
            
            Output:
                List with a dict having coordinates for each gene and its functional annotation if any and the gene index
        """

        #dictionary to stored retrieved coordinates
        coordinates = {}
        gene_index = {}
        count = 0
        current_chr = ""
        is_gene = 0
        gid = ""
        with open(gff_file) as infile:
            for line in infile:
                cols = line.split("\t")
                if len(cols)>6 and cols[2]=="gene":
                    if count==0:
                        current_chr = cols[0]
                    if current_chr!=cols[0]:
                        current_chr = cols[0]
                    
                    for attr in cols[8].split(";"):
                        if attr.startswith("locus_tag="):
                            gid = attr.replace("locus_tag=","").strip()
                        elif attr.startswith("Name="):
                            gid = attr.replace("Name=","").strip()
                        elif attr.startswith("ID="):
                            gid = cols[8].split(";")[0].replace("ID=","").strip()
                        
                    chr = cols[0]
                    start = int(cols[3])
                    end = int(cols[4])
                    ori = cols[6]
                    coordinates[gid]=[chr,start, end, ori]
                    count=1
                    is_gene=1

                    #map each nucleotide to a gene id to speed up downstream processing
                    for i in range(start,end+1):
                        gene_index[f"{cols[0]}\t{i}"]=gid

                if is_gene==1 and len(cols) > 2 and cols[2] in ["rRNA", "tRNA", "CDS", "mRNA", "transcript"]:
                    annotation = ""
                    for attr in cols[8].split(";"):
                        if attr.startswith("product="):
                            annotation = attr.replace("product=","").replace("%28","(").replace("%29",")").replace("%2B","+").replace("%2F","/").replace("%2C",",").replace("%27","'").strip()

                    coordinates[gid].append(annotation)
                    is_gene = 0

        if count==0:
            self.logger.error("Please provide valid GFF3 file for gene coordinates...")
            return []
        else:
            return [coordinates, gene_index]


    def count_insertion_sites(self, target_string: str, search_string: str) -> int:
        """
        Method to count number of insertion sites in a sequence given a search string
        Input:
            target_string: sequence to search (eg gene sequence or genome sequence)
            search_string: string to search for (eg TA for mariner transposon)
        Output:
            Number of sites in target sequence
        """

        searchStrings = []
        complement = search_string.upper().replace("A","t").replace("T","a").replace("C","g").replace("G","c").upper()#complement
        reverse = search_string.upper()[::-1]#reverse
        reverse_complement = search_string.upper().replace("A","t").replace("T","a").replace("C","g").replace("G","c").upper()[::-1]#reverse complement
        
        searchStrings.append(search_string.upper())#original search string
        if complement not in searchStrings:
            searchStrings.append(complement)
        if reverse not in searchStrings:
            searchStrings.append(reverse)
        if reverse_complement not in searchStrings:
            searchStrings.append(reverse_complement)

        sites = []
        for seq in searchStrings:
            sites.extend([m.start() for m in re.finditer('(?='+seq+')', target_string.upper())])
        
        return len(sites)

    @staticmethod
    def split(l, n):
        """
        Method to split a list up evenly into a specified number of smaller lists
        Input:
            l: list
            n: number of chunks to split into
        Output:
            Nested list of input list in N chunks
        """
        k, m = divmod(len(l), n)
        return (l[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
    

    def process_mapped_reads(self, mapped_reads_exp: list, replicon_names: dict, min_hits: int, format: str, logger) -> dict:
        """
        Method to process mapped reads file
        Input:
            mapped_reads_exp: list of names of mapped reads files
            replicon_names: dictionary of replicon/chromosome names
            min_hits: Threshold for minimum number of hits at a unique insertion site for it to be considered a true site
            format: mapped reads format (Bowtie, SOAP or Eland)
        Output:
            dictionary holding a list with dictionary of reads mapped to unique chromosomal locations per sample and the total number of reads per sample
        """
        mapped_reads = {}
        replicon_index = 0
        hit_index = 0
        
        if format == "Bowtie":
            replicon_index = 2
            hit_index = 3
        elif format == "SOAP":
            replicon_index = 7
            hit_index = 8
        elif format == "Eland":
            replicon_index = 6
            hit_index = 7
        else:
            logger.error(f"Unrecognized mapped file format : {format}. Exiting run...")
            sys.exit()

        for reads_file in mapped_reads_exp:
            totalReads = 0
            unique_locations = {}
            with open(reads_file) as infile:
                for line in infile:
                    cols = line.split("\t")
                    if len(line) > 0 and not line.startswith("@") and len(cols) > 2 and cols[replicon_index].strip() != "*":
                        totalReads += 1
                        if cols[replicon_index].strip() not in replicon_names.keys():
                            logger.error(f"Names of replicons (either chromosomes or plasmids) in the sequence file and mapped reads file to do no match!!\n\n {cols[replicon_index].strip()} in mapped reads file but not in sequence file.\n\nExiting run...")
                            sys.exit()

                        unique_locations[f"{cols[replicon_index].strip()}\t{cols[hit_index].strip()}"] = unique_locations.get(f"{cols[replicon_index].strip()}\t{cols[hit_index].strip()}",0)+1

            logger.info(f"Finished processing aligned reads {reads_file}...\n\n Total number of mapped reads : {totalReads}")
                        
            #writing wig file with unnormalized reads
            wig_file = f"{reads_file.split('/')[-1].split('.')[0]}_Reads_per_uniqueLocation.wig"
            logger.info(f"Writing results: {wig_file}")
            with open(wig_file,"w") as outfile:
                outfile.write("track type=wiggle_0 autoScale=on name=\"TnSeq track\" description=\"TnSeq\"\n")
                replicon = ""
                for key in dict(sorted(unique_locations.items(), key=lambda x:(x[0].split("\t")[0],int(x[0].split("\t")[1])))):
                    if replicon != key.split("\t")[0]:
                        replicon = key.split("\t")[0]
                        outfile.write(f"variableStep  chrom={replicon}  span=1\n")

                    if unique_locations[key] >= min_hits:
                        loc = key.split('\t')[1]
                        outfile.write(f"{loc}\t{unique_locations[key]}\n")

            mapped_reads[reads_file]=[unique_locations, totalReads]

        return mapped_reads


    def mp_process_mapped_reads(self, mapped_reads_exp: list, replicon_names: dict, min_hits: int, format: str, logger, threads: int) -> dict:
        """
        Method to multiprocess process_mapped_reads Method
        """
        mapped_reads = {}
        chunks = list(Tools.split(mapped_reads_exp, min(threads,len(mapped_reads_exp))))
        pool = mp.Pool(processes=threads)
        results = pool.starmap(self.process_mapped_reads, zip(chunks,repeat(replicon_names),repeat(min_hits),repeat(format),repeat(logger)))
        pool.close()
        pool.join()

        #combine the individual chunks completed for each thread
        for chunk in results:
            mapped_reads.update(chunk)

        return mapped_reads


    def gene_level_analysis(self, mapped_reads_exp: list, mapped_reads: dict, coordinates: dict, gene_index: dict, target_seq: str, replicon_seqs: dict, min_hits: int, clippings: float, capping: str, analysis_type: str, logger) -> dict:
        """
        Method to process mapped reads file
        Input:
            mapped_reads_exp: list of names of mapped reads files
            mapped_reads: dictionary of mapped reads per samples derived from process_mapped_reads Method
            coordinates: dictionary of gene coordinates for parsing gff file
            gene_index: index mapping each nucleotide to a gene id
            target_seq: transposon specificity sequence (eg TA for mariner transpon)
            replicon_seqs: dictionary of replicon/chromosome sequences
            min_hits: Threshold for minimum number of hits at a unique insertion site for it to be considered a true site
            clippings: Percentage of start and end of gene to be discarded when determining essentiality (e.g., a value of 5 equals 5% to be ignored @ start and @ end)
            capping: Capping of reads at unique sites to minimize PCR/sequencing bias
            analysis_type: one or two sample analysis
            logger: logger

        Output:
            dictionary holding a list with dictionary of gene-level statistics per sample and a list of global (genome-level) per sample
        """
        logger.info("Commencing gene-based analysis...")
        experiment_stats = {}
        
        for reads_file in mapped_reads_exp:
            #Global stats
            total_unique_hits_genome=0#total number of unique hits across entire genome
            total_reads_genome=0#total number of reads across entire genome
            total_unique_hits_genes=0#total number of unique hits occurring within genes/coding regions
            total_reads_within_genes=0#total number of reads occurring within genes/coding regions
            total_norm_unique_hits_per_gene = 0#total normalized unique hits per gene
            max = 0#maximum number of reads at any given insertion site across the genome falling within a gene
            cap = 0

            if analysis_type == "2":
                cap = Tools.set_cap(capping, mapped_reads[reads_file][0], reads_file, logger)

            #Gene level stats
            gene_level_stats = {}
            gene_hits = {}#total number of hits (above threshold) per gene
            total_reads_per_gene = {}#number of reads per gene
            total_reads_per_gene_uncapped = {}#number of reads per gene without any capping

            for location, insertions in dict(sorted(mapped_reads[reads_file][0].items(), key=lambda x:(x[0].split("\t")[0],int(x[0].split("\t")[1])))).items():
                gid = gene_index.get(location,"")
                if gid!="":
                    genesize = coordinates[gid][2]-coordinates[gid][1]+1#end - start + 1
                    clipping = round(genesize*clippings/100)
                    
                    if int(location.split("\t")[1]) >= coordinates[gid][1]+clipping and int(location.split("\t")[1]) <= coordinates[gid][2]-clipping and insertions >= min_hits:
                        total_unique_hits_genes+=1
                        total_reads_within_genes+=insertions
                        gene_hits[gid] = gene_hits.get(gid,0)+1
                        total_reads_per_gene_uncapped[gid] = total_reads_per_gene_uncapped.get(gid,0)+insertions

                        if analysis_type == "2":
                            if capping!='0':
                                if insertions > cap:
                                    total_reads_per_gene[gid] = total_reads_per_gene.get(gid,0)+cap
                                else:
                                    total_reads_per_gene[gid] = total_reads_per_gene.get(gid,0)+insertions
                            else:
                                total_reads_per_gene[gid] = total_reads_per_gene.get(gid,0)+insertions
                        else:
                            total_reads_per_gene[gid] = total_reads_per_gene.get(gid,0)+insertions

                        if insertions > max:
                            max = insertions
                    
                if insertions >= min_hits:
                    total_unique_hits_genome+=1
                    total_reads_genome+=insertions

            for gid in coordinates.keys():#for each gene collated gene level stats
                hits = gene_hits.get(gid, 0)
                genesize = coordinates[gid][2]-coordinates[gid][1]+1#end - start + 1
                clipping = round(genesize*clippings/100)
                available_gene_insertion_sites = 0
                if target_seq !="":
                    available_gene_insertion_sites = self.count_insertion_sites(replicon_seqs[coordinates[gid][0]][coordinates[gid][1]-1+clipping:coordinates[gid][2]-clipping],target_seq)
                else:
                    available_gene_insertion_sites = genesize-(clipping*2)

                norm_unique_hits_per_gene = 0
                if available_gene_insertion_sites > 0:
                    norm_unique_hits_per_gene = hits/available_gene_insertion_sites#number of hits within a gene normalized to total number of available insertion sites
                    total_norm_unique_hits_per_gene+=hits/available_gene_insertion_sites#cumulative normalized value
                

                gene_level_stats[gid] = [hits, total_reads_per_gene.get(gid, 0), genesize-(clipping*2), norm_unique_hits_per_gene, available_gene_insertion_sites, total_reads_per_gene_uncapped.get(gid, 0)]

            global_stats = [total_unique_hits_genome, total_reads_genome, total_unique_hits_genes, total_reads_within_genes, total_norm_unique_hits_per_gene, max]

            logger.info(f"Total number of unique hits in {reads_file} : {total_unique_hits_genome}\nTotal number of reads in {reads_file} : {total_reads_genome}\nTotal number of unique hits within genes in {reads_file} : {total_unique_hits_genes}\nTotal number of reads within genes in {reads_file} : {total_reads_within_genes}\n")

            experiment_stats[reads_file] = [gene_level_stats, global_stats]

        return experiment_stats


    def mp_gene_level_analysis(self, mapped_reads_exp: list, mapped_reads: dict, coordinates: dict, gene_index: dict, target_seq: str, replicon_seqs: dict, min_hits: int, clippings: float, capping: str, analysis_type: str, logger, threads: int) -> list:
        """
        Method to multiprocess gene_level_analysis
        """
        experiment_stats = {}
        chunks = list(Tools.split(mapped_reads_exp, min(threads,len(mapped_reads_exp))))
        pool = mp.Pool(processes=threads)
        results = pool.starmap(self.gene_level_analysis, zip(chunks,repeat(mapped_reads),repeat(coordinates),repeat(gene_index),repeat(target_seq),repeat(replicon_seqs),repeat(min_hits),repeat(clippings),repeat(capping),repeat(analysis_type),repeat(logger)))
        pool.close()
        pool.join()

        #combine the individual chunks completed for each thread
        for chunk in results:
            experiment_stats.update(chunk)

        return experiment_stats


    @staticmethod
    def average_replicates(mapped_reads_exp: list, experiment_stats: dict) -> list:
        """
        Method to compute average statistics over replicates (for one sample analysis)
        Input:
            mapped_reads_exp: list of names of mapped reads files
            experiment_stats: dictionary of stats produced by gene_level_analysis method
        Output:
            List of summary statistics
        """
        mean_unique_hits_genome = 0
        mean_unique_hits_per_gene = {}
        mean_total_reads_per_gene = {}

        for reads_file in mapped_reads_exp:
            mean_unique_hits_genome+=experiment_stats[reads_file][1][0]/len(mapped_reads_exp) #ie avearage of total_unique_hits_genome

            for gid, stats in experiment_stats[reads_file][0].items():
                mean_unique_hits_per_gene[gid] = mean_unique_hits_per_gene.get(gid,0) + (stats[0]/len(mapped_reads_exp)) #ie avearage of hits
                mean_total_reads_per_gene[gid] = mean_total_reads_per_gene.get(gid,0) + (stats[1]/len(mapped_reads_exp)) #ie avearage of total_reads_per_gene

        return [mean_unique_hits_genome, mean_unique_hits_per_gene, mean_total_reads_per_gene]
    
    @staticmethod
    def average_replicates2(mapped_reads_exp: dict, experiment_stats: dict, total_reads_per_gene_weighted: dict) -> list:
        """
        Method to compute average statistics over replicates (for two sample analysis)
        Input:
            mapped_reads_exp: list of names of mapped reads files
            experiment_stats: dictionary of stats produced by gene_level_analysis method
            total_reads_per_gene_weighted: dictionary of dictionary of weighted reads per gene (per sample)
        Output:
            List of summary statistics
        """
        grouped_total_reads_per_gene_weighted = {}
        grouped_total_reads_per_gene = {}
        grouped_total_reads_per_gene_uncapped = {}
        grouped_unique_hits_per_gene = {}
        grouped_total_unique_hits_genes = []
        grouped_totalReads = []

        mean_total_reads_per_gene_weighted = {}
        sd_total_reads_per_gene_weighted = {}
        mean_total_reads_per_gene = {}
        mean_total_reads_per_gene_uncapped = {}
        mean_unique_hits_per_gene = {}
        sd_unique_hits_per_gene = {}
        mean_total_weighted_reads = 0

        #Gather replicate values
        for reads_file in mapped_reads_exp:
            grouped_totalReads.append(mapped_reads_exp[reads_file][1])
            grouped_total_unique_hits_genes.append(experiment_stats[reads_file][1][2])#global stats, total_unique_hits_genes
            for gid, stats in experiment_stats[reads_file][0].items():#gene level stats
                if gid not in grouped_total_reads_per_gene_weighted:
                    grouped_total_reads_per_gene_weighted[gid] = [total_reads_per_gene_weighted[reads_file][gid]]
                    grouped_unique_hits_per_gene[gid] = [stats[0]]#hits
                    grouped_total_reads_per_gene[gid] = [stats[1]]#total_reads_per_gene
                    grouped_total_reads_per_gene_uncapped[gid] = [stats[5]]#total_reads_per_gene_uncapped
                else:
                    grouped_total_reads_per_gene_weighted[gid].append(total_reads_per_gene_weighted[reads_file][gid])
                    grouped_unique_hits_per_gene[gid].append(stats[0])
                    grouped_total_reads_per_gene[gid].append(stats[1])
                    grouped_total_reads_per_gene_uncapped[gid].append(stats[5])
        
        #compute stats
        for gid in grouped_total_reads_per_gene_weighted:
            mean_total_reads_per_gene_weighted[gid] = np.mean(grouped_total_reads_per_gene_weighted[gid])
            sd_total_reads_per_gene_weighted[gid] = np.std(grouped_total_reads_per_gene_weighted[gid], ddof = 1)
            mean_total_weighted_reads+=np.mean(grouped_total_reads_per_gene_weighted[gid])
            mean_total_reads_per_gene[gid] = np.mean(grouped_total_reads_per_gene[gid])
            mean_total_reads_per_gene_uncapped[gid] = np.mean(grouped_total_reads_per_gene_uncapped[gid])
            mean_unique_hits_per_gene[gid] = np.mean(grouped_unique_hits_per_gene[gid])
            sd_unique_hits_per_gene[gid] = np.std(grouped_unique_hits_per_gene[gid], ddof = 1)
        
        mean_total_unique_hits_genes = np.mean(grouped_total_unique_hits_genes)
        mean_totalReads = np.mean(grouped_totalReads)

        return [mean_total_reads_per_gene_weighted, sd_total_reads_per_gene_weighted, mean_total_reads_per_gene, mean_total_reads_per_gene_uncapped, mean_unique_hits_per_gene, sd_unique_hits_per_gene, mean_total_unique_hits_genes, mean_totalReads, mean_total_weighted_reads]

    @staticmethod
    def set_cap(capping: str, unique_locations: dict, reads_file: str, logger) -> list:
        """
        Method to capped reads on unique insertion sites
        Input:
            capping: capping option to use (0 - no capping; or 1 - capping with average hits per insertion site + 2 st. dev. (default); or 2 - capping with average hits per insertion site; or 3 - capping with median hits per insertion site)
            unique_locations: unique location dictionary produced by process_mapped_reads function
            reads_file: name of mapped reads file
        Output:
            Value to cap maximum reads at unique location at
        """
        
        logger.info("Capping reads...\n")
        cap = 0
        mean = np.mean(list(unique_locations.values()))#average number of hits per location
        sd = np.std(list(unique_locations.values()), ddof = 1)#standard deviation of hits per location
        median = np.median(list(unique_locations.values()))#median hits per location
        mean2SD = mean + 2*sd

        if capping == "1":
            cap = round(mean2SD,0)
        elif capping == "2":
            cap = round(mean,0)
        else:
            cap = round(median,0)

        logger.info(f"{reads_file}\nMean: {mean}\tMedian: {median}\tMean2SD: {mean2SD}\nRead cap set at: {cap} reads per unique location")
    
        return cap
    
    @staticmethod
    def weight_reads(experiment_stats: dict, mapped_reads: dict, logger) -> list:
        """
        Method to weight reads based on number of unique insertions (for two sample analysis)
        Input:
            experiment_stats: dictionary of stats produced by gene_level_analysis method
            mapped_reads_exp: list of names of mapped reads files
            logger: logger
        Output:
            Dictionary of dictionary of weighted reads per gene (per sample)
        """

        logger.info("Weighting reads...\n")
        norm_hits = 0
        total_reads_per_gene_weighted = {}
        for reads_file in mapped_reads:
            norm_hits+=experiment_stats[reads_file][1][4]#Global stats,total_norm_unique_hits_per_gene

        ave_total_norm_hits = norm_hits/(len(mapped_reads))

        for reads_file in mapped_reads:
            weighted_reads_per_gene = {}
            for gid in experiment_stats[reads_file][0].keys():#gene level stats keys (gids)
                total_reads_per_gene = experiment_stats[reads_file][0][gid][1]#Gene level stats, gid,total_reads_per_gene
                norm_unique_hits_per_gene = experiment_stats[reads_file][0][gid][3]#Gene level stats, gid,norm_unique_hits_per_gene
                weighted = 0
                try:
                    weighted = total_reads_per_gene*(norm_unique_hits_per_gene/(ave_total_norm_hits/len(experiment_stats[reads_file][0])))# this weight total reads per gene to based on the number of insertions in that gene relative to the rate of insertion ger gene across all experiments
                except Exception as e:
                    pass #catching and ignoring divide by zero exceptions. These values set to zero
                if weighted > 0 and weighted < 1.0:
                    weighted = 1.0

                weighted_reads_per_gene[gid] = weighted

            total_reads_per_gene_weighted[reads_file]=weighted_reads_per_gene

        return total_reads_per_gene_weighted

    @staticmethod
    def compute_tsat(num_of_samples_c: int, num_of_samples_t: int, sd_unique_hits_c: dict, sd_unique_hits_t: dict, mean_unique_hits_c: dict, mean_unique_hits_t: dict, exclude: int) -> list:
        """
        Method to conduct t-test (for two sample analysis)
        Input:
            num_of_samples_c: number of samples in control group
            num_of_samples_t: number of samples in treatment group
            sd_unique_hits_c: standard deviation of unique hits (or reads) in control group
            sd_unique_hits_t: standard deviation of unique hits (or reads) in treatment group
            mean_unique_hits_c: mean of unique hits (or reads) in control group
            mean_unique_hits_t: mean of unique hits (or reads) in treatment group
            exclude: number of genes to exclude from correction for multiple testing (ie gene with no insertions)
        Output:
            List of dictionaries of pvalues and adjusted pvalues
        """
        
        grand_sd = {}
        tstat = {}
        pvalues = {}
        adj_pvalues = {}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            #compute grand SD for each gene
            for gid in sd_unique_hits_c:
                #(((number of treatment samples*std of number unique hits treat)+(number of control samples*std of number unique hits ctl))/total number of all samples - 2)**5
                grand_sd[gid]=((((num_of_samples_t-1)*(sd_unique_hits_t[gid]**2))+((num_of_samples_c-1)*(sd_unique_hits_c[gid]**2)))/(num_of_samples_t+num_of_samples_c-2))**0.5

            #compute t-statistic for each gene
            for gid in sd_unique_hits_c:
                #(mean_unique_hits_per_gene_treat - mean_unique_hits_per_gene_ctl)/(grand_sd*sqrt((1/number of treatment samples)+(1/number of ctl samples)))
                tstat[gid] = (mean_unique_hits_t[gid]-mean_unique_hits_c[gid])/(grand_sd[gid]*((1/num_of_samples_t)+(1/num_of_samples_c))**0.5)

            #compute pvalues
            for gid in sd_unique_hits_c:
                if np.isnan(tstat[gid]):
                    pvalues[gid] = 1
                if mean_unique_hits_t[gid] <=2 and mean_unique_hits_c[gid] <=2:
                    pvalues[gid] = 1
                else:
                    pvalues[gid] = stats.t.sf(abs(tstat[gid]),(num_of_samples_t+num_of_samples_c-2))

            #Modified BH correction
            adj_pvalues = Tools.correct_pvalues(pvalues, exclude)[0]

        return [pvalues, adj_pvalues]


    @staticmethod
    def compute_ratios(summary_stats_c: list, summary_stats_t: list) -> list:
        """
        Method to compute ratios (for two sample analysis)
        Input:
            summary_stats_c: list of summary stats per gene for control group obtain from average_replicates2
            summary_stats_t: list of summary stats per gene for treatment group obtain from average_replicates2
        Output:
            List of dictionaries of ratios of reads and insertions, along with the number of genes to exclude from correction for multiple testing (ie gene with no insertions)
        """
        ratio_reads = {}
        ratio_ins = {}
        exclude = 0

        for gid in summary_stats_c[0]:#mean_total_reads_per_gene_weighted
            #ratio for reads
            if (summary_stats_c[0][gid] == 0 and summary_stats_t[0][gid] == 0) or (summary_stats_c[0][gid] == 1 and summary_stats_t[0][gid] == 1):
                ratio_reads[gid] = 1
            elif summary_stats_c[0][gid] == 0:
                #(mean_total_reads_per_gene_weighted treat/mean_totalReads treat)/(1.0/mean_totalReads ctl);
                ratio_reads[gid] = (summary_stats_t[0][gid]/summary_stats_t[7])/(0.5/summary_stats_c[7])
            elif summary_stats_t[0][gid] == 0:
                ratio_reads[gid] = (0.5/summary_stats_t[7])/(summary_stats_c[0][gid]/summary_stats_c[7])
            else:
                ratio_reads[gid] = (summary_stats_t[0][gid]/summary_stats_t[7])/(summary_stats_c[0][gid]/summary_stats_c[7])

            #ratio for insertions
            if (summary_stats_c[4][gid] == 0 and summary_stats_t[4][gid] == 0):
                ratio_ins[gid] = 1
                exclude+=1
            elif (summary_stats_c[4][gid] == 1 and summary_stats_t[4][gid] == 1):
                ratio_ins[gid] = 1
            elif summary_stats_c[0][gid] == 0:
                #(mean_unique_hits_per_gene treat/mean_total_unique_hits_genes treat)/(1.0/mean_total_unique_hits_genes ctl);
                ratio_ins[gid] = (summary_stats_t[4][gid]/summary_stats_t[6])/(0.1/summary_stats_c[6])
            elif summary_stats_t[0][gid] == 0:
                ratio_ins[gid] = (0.1/summary_stats_t[6])/(summary_stats_c[4][gid]/summary_stats_c[6])
            else:
                ratio_ins[gid] = (summary_stats_t[4][gid]/summary_stats_t[6])/(summary_stats_c[4][gid]/summary_stats_c[6])

        return [ratio_reads, ratio_ins, exclude]


    @staticmethod
    def compute_proportions(mean_stat_per_gene_t: dict, mean_stat_per_gene_c: dict, mean_total_stat_t: float, mean_total_stat_c: float, exclude: int) -> list:
        """
        Method to compute proportions pvalues (for two sample analysis)
        Input:
            mean_stat_per_gene_c: mean of unique hits (or reads) per gene in control group
            mean_stat_per_gene_t: mean of unique hits (or reads) per gene in treatment group
            mean_total_stat_c: mean of total unique hits (or reads) in control group
            mean_total_stat_t: mean of total unique hits (or reads) in treatment group
            exclude: number of genes to exclude from correction for multiple testing (ie gene with no insertions)
        Output:
            List of dictionaries of pvalues and adjusted pvalues
        """
        pvalues = {}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            for gid in mean_stat_per_gene_c:
                #insertions
                try:
                    diff = (mean_stat_per_gene_t[gid]/mean_total_stat_t) - (mean_stat_per_gene_c[gid]/mean_total_stat_c)
                    cum_prop = (mean_stat_per_gene_t[gid]+mean_stat_per_gene_c[gid])/(mean_total_stat_t+mean_total_stat_c)
                    prop_test_stat = diff/((cum_prop*(1-cum_prop))*((1/mean_total_stat_t)+(1/mean_total_stat_c)))**0.5
                except Exception as e:
                    prop_test_stat = np.nan

                if np.isnan(prop_test_stat):
                    pvalues[gid] = 1
                else:
                    pvalues[gid] = stats.norm.sf(abs(prop_test_stat))*2 #two-tailed

        #Modified BH correction
        adj_pvalues = Tools.correct_pvalues(pvalues, exclude)[0]
    
        return [pvalues,adj_pvalues]
    
    @staticmethod
    def compute_fisher(mean_stat_per_gene_t: dict, mean_stat_per_gene_c: dict, mean_total_stat_t: float, mean_total_stat_c: float, exclude: int) -> list:
        """
        Method to compute fisher test pvalues (for two sample analysis)
        Input:
            mean_stat_per_gene_c: mean of unique hits (or reads) per gene in control group
            mean_stat_per_gene_t: mean of unique hits (or reads) per gene in treatment group
            mean_total_stat_c: mean of total unique hits (or reads) in control group
            mean_total_stat_t: mean of total unique hits (or reads) in treatment group
            exclude: number of genes to exclude from correction for multiple testing (ie gene with no insertions)
        Output:
            List of dictionaries of pvalues and adjusted pvalues
        """
        pvalues = {}

        for gid in mean_stat_per_gene_c:
            if (mean_total_stat_t - mean_stat_per_gene_t[gid] < 0) or (mean_total_stat_c - mean_stat_per_gene_c[gid] < 0):
                pvalues[gid] = 1
            else:
                a = round(mean_stat_per_gene_t[gid],0)
                b = round(mean_stat_per_gene_c[gid],0)
                c = round(mean_total_stat_t - mean_stat_per_gene_t[gid],0)
                d = round(mean_total_stat_c - mean_stat_per_gene_c[gid],0)
                table = [[a,b],[c,d]]
                if a == 0 and b !=0:#set min count to 1
                    table = [[1,b],[c,d]]
                elif a != 0 and b ==0:#set min count to 1
                    table = [[a,1],[c,d]]

                if a==0 and b==0:
                    pvalues[gid] = 1
                else:
                    #pvalues[gid] = stats.fisher_exact(table,alternative='two-sided')[1]
                    pvalues[gid] = stats.chi2_contingency(table)[1]#switched from fisher to chi2 for speed (due to large numbers). Results very similar
        
        #Modified BH correction
        adj_pvalues = Tools.correct_pvalues(pvalues, exclude)[0]

        return [pvalues,adj_pvalues]
    
    @staticmethod
    def correct_pvalues(pvalues: dict, exclude: int) -> list:
        """
        Method to compute fisher test pvalues (for two sample analysis)
        Input:
            pvalues: dictionary with pvalues to be corrected
            exclude: number of genes to exclude from correction for multiple testing (ie gene with no insertions)
        Output:
            List of dictionaries of pvalues and adjusted pvalues
        """
        adj_pvalues = {}
        bonferroni_pvalues = {}
        denominator = 1

        for gid, pval in dict(sorted(pvalues.items(), key=lambda x:(x[1]))).items():
            adj_pvalues[gid] = min(1,((len(pvalues)-exclude)/denominator)*pval)
            bonferroni_pvalues[gid] = min(1,len(pvalues)*pval)

            if denominator < (len(pvalues)-exclude):
                denominator+=1

        return [adj_pvalues,bonferroni_pvalues]