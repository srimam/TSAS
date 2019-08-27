package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import jsc.distributions.*;
import jsc.contingencytables.*;
import jsc.tests.H1;
/*
Code to analyze TnSeq data, comparing two samples (control and treatment)

usage: java TnSeq_analysis_twoSample
*/

public class TnSeq_analysis_twoSample
{
	@SuppressWarnings("unchecked")
	public static void execute() throws IOException
	{
	try
	{
		Calendar s = Calendar.getInstance();
		System.out.println(s.getTime());
		
		Pattern p = Pattern.compile(">");
		Pattern p1 = Pattern.compile("\t");
		Pattern p2 = Pattern.compile(" ");
		Pattern p3 = Pattern.compile(";");
		Pattern p4 = Pattern.compile("locus_tag=");
		Pattern p5 = Pattern.compile("Name=");
		Pattern p6 = Pattern.compile("=");
		Pattern p7 = Pattern.compile(",");
		Pattern p8 = Pattern.compile("product=");
		Pattern p9 = Pattern.compile("ID=");
		Pattern p10 = Pattern.compile("[|]");


//Reading parameters file
		FileInputStream f3 = new FileInputStream(new File("Parameters.txt"));
		BufferedReader b3 = new BufferedReader(new InputStreamReader(f3));


//Processing Parameters.txt file
		String lines="",genome_fasta="",gff="",min_hits="",clippings="",capping="",weights="",format="",result_format="";
		ArrayList<String> mapped_reads_reps_ctl=new ArrayList<String>();//**
		ArrayList<String> mapped_reads_reps_exp=new ArrayList<String>();//**
		System.out.println("\nInput Parameters\n");
		while ((lines = b3.readLine())!=null)
		{	
			if(lines.length()>0 && lines.charAt(0)!='#')
			{
				String[] result = p6.split(lines);
				if(result[0].trim().equalsIgnoreCase("Genome_sequence"))
				{ 
					if(result.length>1) genome_fasta = result[1].trim();
					else genome_fasta = "";
					System.out.println("Name of genome FASTA file : "+genome_fasta);
				}
				else if(result[0].trim().equalsIgnoreCase("GFF")) 
				{ 
					if(result.length>1) gff = result[1].trim();
					else gff = "";
					System.out.println("Name of gff file : "+gff);
				}
				else if(result[0].trim().equalsIgnoreCase("Mapped_format")) 
				{ 
					if(result.length>1)format = result[1].trim();
					else format = "";
					System.out.println("Mapped reads file format : "+format);
				}
				else if(result[0].trim().equalsIgnoreCase("Control")) 
				{ 
					String[] result1;
					if(result.length>1) result1 = p7.split(result[1].trim());//**
					else result1=new String[]{};
					int counting=0;
					for(int i=0;i<result1.length;i++)
					{
						mapped_reads_reps_ctl.add(result1[i].trim());
						if(counting==0)
							System.out.println("Name of mapped reads file (Control) : "+mapped_reads_reps_ctl.get(i));
						else
							System.out.print(", "+mapped_reads_reps_ctl.get(i));
						counting++;
					}
				}
				else if(result[0].trim().equalsIgnoreCase("Treatment")) 
				{ 
					String[] result1;
					if(result.length>1) result1 = p7.split(result[1].trim());//**
					else result1=new String[]{};
					int counting=0;
					for(int i=0;i<result1.length;i++)
					{
						mapped_reads_reps_exp.add(result1[i].trim());
						if(counting==0)
							System.out.println("\nName of mapped reads file (Treatment) : "+mapped_reads_reps_exp.get(i));
						else
							System.out.print(", "+mapped_reads_reps_exp.get(i));
						counting++;
					}
				}
				else if(result[0].trim().equalsIgnoreCase("Min. hits")) 
				{ 
					if(result.length>1) min_hits = result[1].trim();
					else min_hits = "";
					try{Double.parseDouble(min_hits);}
					catch(Exception e){min_hits="5";}
					System.out.println("\nThreshold for unique insertions : "+min_hits +" hits");
				}
				else if(result[0].trim().equalsIgnoreCase("Clipping")) 
				{ 
					if(result.length>1) clippings = result[1].trim();
					else clippings = "";
					try{Double.parseDouble(clippings);}
					catch(Exception e){clippings="5";}
					System.out.println("% of start and end to ignore : "+clippings+" %");
				}
				else if(result[0].trim().equalsIgnoreCase("Capping")) 
				{ 
					if(result.length>1) capping = result[1].trim();
					else capping = "";
					if(capping.equals("0")) System.out.println("Capping : No");
					else if(capping.equals("1")) System.out.println("Capping : Yes (mean + 2SD)");
					else if(capping.equals("2")) System.out.println("Capping : Yes (mean)");
					else if(capping.equals("3")) System.out.println("Capping : Yes (median)");
					else
					{
						capping = "1";
						System.out.println("Capping : Yes (Capping set to default value of 1 (mean + 2SD))");
					}
				}
				else if(result[0].trim().equalsIgnoreCase("Weights")) 
				{ 
					if(result.length>1) weights = result[1].trim();
					else weights = "";
					if(weights.equals("0")) System.out.println("Weights : No");
					else
					{
						weights = "1";
						System.out.println("Weights : Yes");
					}
				}
				else if(result[0].trim().equalsIgnoreCase("Result")) 
				{ 
					if(result.length>1) result_format = result[1].trim();
					else result_format = "";
					if(result_format.equalsIgnoreCase("long")) System.out.println("Result format : Long");
					else
					{
						result_format = "short";
						System.out.println("Result format : short");
					}
				}
			}
		}
		

//Checking file validity
		
		File newFile;
		//File containing fasta format of genome
		if(genome_fasta.equals(""))
		{
			System.out.println("Please provide a valid genome sequence FASTA file... Exiting run!!");
			return;
		}
		newFile = new File(genome_fasta);
		if(!newFile.exists())
		{
			System.out.println("\n\n"+genome_fasta+ " (fasta file) not found... Exiting run!!");
			return;
		}

		//GFF file
		if(gff.equals(""))
		{
			System.out.println("Please provide a valid GFF file... Exiting run!!");
			return;
		}
		newFile = new File(gff);
		if(!newFile.exists())
		{
			System.out.println("\n\n"+gff+ " (gff file) not found... Exiting run!!");
			return;
		}

		//Mapped file format
		if(!(format.equalsIgnoreCase("Bowtie") || format.equalsIgnoreCase("SOAP") || format.equalsIgnoreCase("Eland")))
		{
			System.out.println("Unrecognized mapped file format ("+format+")... Exiting run!!");
			return;
		}

		//Aligned reads control
		if(mapped_reads_reps_ctl.size()==0)
		{
			System.out.println("Please provide a valid mapped reads file (control)... Exiting run!!");
			return;
		}

		//Aligned reads Treatment
		if(mapped_reads_reps_exp.size()==0)
		{
			System.out.println("Please provide a valid mapped reads file (treatment)... Exiting run!!");
			return;
		}
		
		//Make sure different file names are specified for treatment and control and that the specified files actually exist...
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			if(mapped_reads_reps_ctl.contains(mapped_reads_reps_exp.get(i)))
			{
				System.out.println("Some control and treatment files have identical names. Please provide a different control or treatment files... Exiting run!!");
				return;
			}
			newFile = new File(mapped_reads_reps_exp.get(i));
			if(!newFile.exists())
			{
				System.out.println("\n\n"+mapped_reads_reps_exp.get(i)+ " (mapped reads treatment file) not found... Exiting run!!");
				return;
			}
		}

		for(int i=0;i<mapped_reads_reps_ctl.size();i++)
		{
			newFile = new File(mapped_reads_reps_ctl.get(i));
			if(!newFile.exists())
			{
				System.out.println("\n\n"+mapped_reads_reps_ctl.get(i)+ " (mapped reads control file) not found... Exiting run!!");
				return;
			}
		}

		
//Parsing fasta file for genome...
		FileInputStream f = new FileInputStream(new File(genome_fasta));
		BufferedReader b = new BufferedReader(new InputStreamReader(f));
		ArrayList<String> replicon_names = new ArrayList<String>();
		ArrayList<String> replicon_seqs = new ArrayList<String>();
		ArrayList<String> replicon_desc = new ArrayList<String>();
		StringBuffer allSeq= new StringBuffer();
		lines="";
		while ((lines = b.readLine())!=null)
		{	
			if(lines.charAt(0)=='>')
			{
				String[] result = p2.split(lines);
				replicon_names.add(result[0].trim().replace(">",""));
				replicon_desc.add(lines.trim().replace(">",""));
			}
		}
		f.close();b.close();
		System.out.println("\nNumber of replicons : "+replicon_names.size());
		
		f = new FileInputStream(new File(genome_fasta));
		b = new BufferedReader(new InputStreamReader(f));
		lines="";
		while ((lines = b.readLine())!=null)
		{	
			allSeq.append(lines);
		}
		f.close();b.close();

		String[] seqResult = p.split(allSeq);
		for(int i = 0;i<replicon_names.size();i++)
		{
			replicon_seqs.add(seqResult[i+1].replace(replicon_desc.get(i),""));
		}
		int genomeSize=0;
		for(int i = 0;i<replicon_names.size();i++)
		{
			System.out.println(replicon_names.get(i)+"\t"+replicon_seqs.get(i).length() + " bp");
			genomeSize+=replicon_seqs.get(i).length();
		}
		System.out.println("Genome size : "+genomeSize+ " bp\n");


//Parsing gff file for gene coordinates and annotation...

		FileInputStream f1 = new FileInputStream(new File(gff));
		BufferedReader b1 = new BufferedReader(new InputStreamReader(f1));
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> annotation = new ArrayList<String>();
		ArrayList<String> replicon = new ArrayList<String>();
		ArrayList<Integer> end = new ArrayList<Integer>();
		ArrayList<Integer> start = new ArrayList<Integer>();
		lines = "";
		int tracker=0;
		while ((lines = b1.readLine())!=null)
		{
			String[] result = p1.split(lines);
			if(result.length>2 && result[2].equals("gene"))
			{
				if(!replicon_names.contains(result[0]))
				{
					System.out.println("Names of replicons (either chromosomes or plasmids) in the sequence file and gff file to do no match!!\n\n"+ result[0] + " in gff file but not in sequence file.\n\n"+"Exiting run...");
					return;
				}
				replicon.add(result[0]);
				String[] result2 = p3.split(result[8]);
				Matcher m = p4.matcher(result[8]);
				Matcher m2 = p5.matcher(result[8]);
				if(m.find())
				{
					for(int i = 0;i<result2.length;i++)
					{
						m=p4.matcher(result2[i]);
						if(m.find())
						{
							geneID.add(result2[i].trim().replace("locus_tag=","").replace("\"",""));
							i=result2.length;
							tracker++;
						}
					}
				}
				else if(m2.find())
				{
					for(int i = 0;i<result2.length;i++)
					{
						m=p5.matcher(result2[i]);
						if(m.find())
						{
							geneID.add(result2[i].trim().replace("Name=",""));
							i=result2.length;
							tracker++;
						}
					}
				}
				else
				{
					for(int i = 0;i<result2.length;i++)
					{
						m=p9.matcher(result2[i]);
						if(m.find())
						{
							geneID.add(result2[i].trim().replace("ID=",""));
							i=result2.length;
							tracker++;
						}
					}
				}
				if(result[6].trim().equals("-"))
				{
					end.add(Integer.parseInt(result[4]));
					start.add(Integer.parseInt(result[3]));
				}
				else
				{
					end.add(Integer.parseInt(result[4]));
					start.add(Integer.parseInt(result[3]));
				}
				if(tracker==2)
				{
					annotation.add("Pseudo gene");
					tracker=1;
				}
				
			}
			if(tracker==1)
			{
				if(result.length>2 && (result[2].equals("rRNA") || result[2].equals("tRNA") || result[2].equals("CDS") || result[2].equals("mRNA") || result[2].equals("transcript")))
				{
					String[] result2 = p3.split(result[8]);
					int find=0;
					for(int i = 0;i<result2.length;i++)
					{
						Matcher m=p8.matcher(result2[i]);
						if(m.find())
						{
							annotation.add(result2[i].trim().replace("product=","").replace("%28","(").replace("%29",")").replace("%2B","+").replace("%2F","/").replace("%2C",",").replace("%27","'"));
							i=result2.length;
							find++;
						}
					}
					if(find==0)
						annotation.add("");
					tracker=0;
				}
			}
		}
		System.out.println("Finished parsing GFF...\n");
		f1.close();b1.close();
		
		
//Processing mapped reads files
	
		//Parsing and processing mapped reads (Treatment)
		System.out.println("Processing aligned reads (Treatment)...\n");
		ArrayList<ArrayList> unique_locations=new ArrayList<ArrayList>();
		ArrayList<ArrayList> hits_per_unique_location=new ArrayList<ArrayList>();
		int[] total_reads=new int[mapped_reads_reps_exp.size()];
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			FileInputStream f4 = new FileInputStream(new File(mapped_reads_reps_exp.get(i)));
			BufferedReader b4 = new BufferedReader(new InputStreamReader(f4));
			ArrayList<Object> treatment_results =  TnSeq_analysis_twoSample.process_mapped_reads(b4,replicon_names,replicon_seqs,mapped_reads_reps_exp.get(i),min_hits,format);
			if(treatment_results.size()==0) return;
			else
			{
				unique_locations.add((ArrayList<String>)treatment_results.get(0));
				hits_per_unique_location.add((ArrayList<Integer>)treatment_results.get(1));
				total_reads[i] = (Integer)treatment_results.get(3);
			}
			f4.close();b4.close();
		}

		//Parsing and processing mapped reads (Control)
		System.out.println("Processing aligned reads (Control)...\n");
		ArrayList<ArrayList> unique_locations_ctl=new ArrayList<ArrayList>(); 
		ArrayList<ArrayList> hits_per_unique_location_ctl=new ArrayList<ArrayList>();
		int[] total_reads_ctl=new int[mapped_reads_reps_ctl.size()];
		for(int i=0;i<mapped_reads_reps_ctl.size();i++)
		{
			FileInputStream f2 = new FileInputStream(new File(mapped_reads_reps_ctl.get(i)));
			BufferedReader b2 = new BufferedReader(new InputStreamReader(f2));
			ArrayList<Object> ctl_results =  TnSeq_analysis_twoSample.process_mapped_reads(b2,replicon_names,replicon_seqs,mapped_reads_reps_ctl.get(i),min_hits,format);
			if(ctl_results.size()==0) return;
			else
			{
				unique_locations_ctl.add((ArrayList<String>)ctl_results.get(0));
				hits_per_unique_location_ctl.add((ArrayList<Integer>)ctl_results.get(1));
				total_reads_ctl[i] = (Integer)ctl_results.get(3);
			}
			f2.close();b2.close();
		}

	
//Identify hits per gene

		//Analyzing Treatment reads
		System.out.println("Analyzing Treatment reads...\n");
		ArrayList<ArrayList> TotalReads_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> uniqueHits_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> norm_uniqueHits_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> GeneSize=new ArrayList<ArrayList>();
		int[] max =new int[mapped_reads_reps_exp.size()];
		int[] total_uniqueHits_genes=new int[mapped_reads_reps_exp.size()];
		int[] total_reads_genes=new int[mapped_reads_reps_exp.size()];
        	double[] total_norm_uniqueHits_per_gene=new double[mapped_reads_reps_exp.size()];
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			ArrayList<Object> treatment_results2 = TnSeq_analysis_twoSample.gene_level_analysis(min_hits,clippings,mapped_reads_reps_exp.get(i),end,start,unique_locations.get(i),hits_per_unique_location.get(i),replicon,geneID,annotation);

			if(treatment_results2.size()==0) 
			{
				System.out.println("Error occurred while performing gene-level analysis of treatment data... Exiting run");
				return;
			}
			else
			{
				uniqueHits_per_gene.add((ArrayList<Double>) treatment_results2.get(0));
				norm_uniqueHits_per_gene.add((ArrayList<Double>) treatment_results2.get(1));
				TotalReads_per_gene.add((ArrayList<Double>) treatment_results2.get(2));
				GeneSize.add((ArrayList<Integer>) treatment_results2.get(3));
				max[i] = (Integer)treatment_results2.get(4);
				total_reads_genes[i] = (Integer)treatment_results2.get(5);
				total_uniqueHits_genes[i] = (Integer)treatment_results2.get(6);
            			total_norm_uniqueHits_per_gene[i] =(Double)treatment_results2.get(7);
			}
		}

		//Analyzing Control reads
		System.out.println("Analyzing Control reads...\n");
		ArrayList<ArrayList> TotalReads_per_gene_ctl=new ArrayList<ArrayList>();
		ArrayList<ArrayList> uniqueHits_per_gene_ctl=new ArrayList<ArrayList>();
		ArrayList<ArrayList> norm_uniqueHits_per_gene_ctl=new ArrayList<ArrayList>();
		ArrayList<ArrayList> GeneSize_ctl=new ArrayList<ArrayList>();
		int[] max_ctl =new int[mapped_reads_reps_ctl.size()];
		int[] total_uniqueHits_genes_ctl=new int[mapped_reads_reps_ctl.size()];
		int[] total_reads_genes_ctl=new int[mapped_reads_reps_ctl.size()];
        	double[] total_norm_uniqueHits_per_gene_ctl=new double[mapped_reads_reps_ctl.size()];
		for(int i=0;i<mapped_reads_reps_ctl.size();i++)
		{
			ArrayList<Object> ctl_results2 = TnSeq_analysis_twoSample.gene_level_analysis(min_hits,clippings,mapped_reads_reps_ctl.get(i),end,start,unique_locations_ctl.get(i),hits_per_unique_location_ctl.get(i),replicon,geneID,annotation);

			if(ctl_results2.size()==0) 
			{
				System.out.println("Error occurred while performing gene-level analysis of control data... Exiting run");
				return;
			}
			else
			{
				uniqueHits_per_gene_ctl.add((ArrayList<Double>) ctl_results2.get(0));
				norm_uniqueHits_per_gene_ctl.add((ArrayList<Double>) ctl_results2.get(1));
				TotalReads_per_gene_ctl.add((ArrayList<Double>) ctl_results2.get(2));
				GeneSize_ctl.add((ArrayList<Integer>) ctl_results2.get(3));
				max_ctl[i] = (Integer)ctl_results2.get(4);
				total_reads_genes_ctl[i] = (Integer)ctl_results2.get(5);
				total_uniqueHits_genes_ctl[i] = (Integer)ctl_results2.get(6);
            			total_norm_uniqueHits_per_gene_ctl[i] =(Double)ctl_results2.get(7);
			}
		}


//Capping analysis...

		ArrayList<ArrayList> rawReads=new ArrayList<ArrayList>();
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			rawReads.add(TotalReads_per_gene.get(i));
			if(!capping.equals("0"))
			{
				TotalReads_per_gene.set(i,TnSeq_analysis_twoSample.capping(capping,weights,clippings,min_hits,max[i],total_reads[i],unique_locations.get(i), hits_per_unique_location.get(i),geneID,replicon,end,start));
			}		
		}
		ArrayList<ArrayList> rawReads_ctl=new ArrayList<ArrayList>();
		for(int i=0;i<mapped_reads_reps_ctl.size();i++)
		{
			rawReads_ctl.add(TotalReads_per_gene_ctl.get(i));
			if(!capping.equals("0"))
			{
				TotalReads_per_gene_ctl.set(i,TnSeq_analysis_twoSample.capping(capping,weights,clippings,min_hits,max_ctl[i],total_reads_ctl[i],unique_locations_ctl.get(i),hits_per_unique_location_ctl.get(i),geneID,replicon,end,start));
			}		
		}


//Weighting reads...

		ArrayList[] TotalReads_per_gene_weighted = new ArrayList[mapped_reads_reps_exp.size()];
		ArrayList[] TotalReads_per_gene_weighted_ctl = new ArrayList[mapped_reads_reps_ctl.size()];
		double counting = 0.0;double norm_hits=0.0;
		for(int i=0;i<mapped_reads_reps_ctl.size();i++)
		{
			norm_hits+=total_norm_uniqueHits_per_gene_ctl[i];
			counting++;
		}
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			norm_hits+=total_norm_uniqueHits_per_gene[i];
			counting++;
		}
		double ave_norm_hits=norm_hits/counting;
		if(!weights.equals("0"))
		{
			for(int j=0;j<TotalReads_per_gene.size();j++)
			{
				TotalReads_per_gene_weighted[j]=new ArrayList<Double>();
				for(int i =0;i<TotalReads_per_gene.get(j).size();i++)
				{	
                			double weighted = ((Double)TotalReads_per_gene.get(j).get(i)*((Double)norm_uniqueHits_per_gene.get(j).get(i)/(ave_norm_hits/(double)geneID.size())));
               				 if(weighted<1.0)
                   		 		TotalReads_per_gene_weighted[j].add(1.0);//treatment
                			else
                   				TotalReads_per_gene_weighted[j].add(weighted);
				}
			}
			for(int j=0;j<TotalReads_per_gene_ctl.size();j++)
			{
				TotalReads_per_gene_weighted_ctl[j]=new ArrayList<Double>();
				for(int i =0;i<TotalReads_per_gene_ctl.get(j).size();i++)
				{	
                			double weighted = ((Double)TotalReads_per_gene_ctl.get(j).get(i)*((Double)norm_uniqueHits_per_gene_ctl.get(j).get(i)/(ave_norm_hits/(double)geneID.size())));
                			if(weighted<1.0)
                    				TotalReads_per_gene_weighted_ctl[j].add(1.0);
                			else
                   		 		TotalReads_per_gene_weighted_ctl[j].add(weighted);//control
				}
			}
		}
		else
		{
			for(int i =0;i<TotalReads_per_gene.size();i++)
			{
				TotalReads_per_gene_weighted[i]=TotalReads_per_gene.get(i);//treatment
			}
			for(int i =0;i<TotalReads_per_gene_ctl.size();i++)
			{
				TotalReads_per_gene_weighted_ctl[i]=TotalReads_per_gene_ctl.get(i);//control
			}
		}
	
	
//Calculating difference in insertions between samples and statistical significance...

		System.out.println("fCalculating gene level fold enrichment and significance between treatment and control...\n");
		ArrayList<Double>  average_exp = new ArrayList<Double>();
		ArrayList<Double>  average_ctl = new ArrayList<Double>();
		ArrayList<Double>  average_tot_reads_exp = new ArrayList<Double>();
		ArrayList<Double>  average_tot_reads_ctl = new ArrayList<Double>();
		ArrayList<Double>  sd_exp = new ArrayList<Double>();
		ArrayList<Double>  sd_ctl = new ArrayList<Double>();
		ArrayList<Double>  sd_ins_exp = new ArrayList<Double>();
		ArrayList<Double>  sd_ins_ctl = new ArrayList<Double>();
		ArrayList<Double>  average_uniqueHits_exp = new ArrayList<Double>();
		ArrayList<Double>  average_uniqueHits_ctl = new ArrayList<Double>();
		ArrayList<Double>  average_rawReads_exp = new ArrayList<Double>();
		ArrayList<Double>  average_rawReads_ctl = new ArrayList<Double>();
		ArrayList<Double>  average_capReads_exp = new ArrayList<Double>();
		ArrayList<Double>  average_capReads_ctl = new ArrayList<Double>();
		ArrayList<Double>  average_unique_hits_genes_exp = new ArrayList<Double>();
		ArrayList<Double>  average_unique_hits_genes_ctl = new ArrayList<Double>();
		for(int i = 0;i<geneID.size();i++)
		{
			double count1=0.0,temp=0.0,temp1=0.0,temp2=0.0,temp3=0.0,temp4=0.0,temp5=0.0,X=0.0,Y=0.0;
			for(int j=0;j<mapped_reads_reps_exp.size();j++)
			{
				temp+=(Double)TotalReads_per_gene_weighted[j].get(i);
				temp1+=(double)total_reads[j];
				temp2+=(Double)uniqueHits_per_gene.get(j).get(i);
				temp3+=(Double)rawReads.get(j).get(i);
				temp4+=(Double)TotalReads_per_gene.get(j).get(i);
				temp5+=(double)total_uniqueHits_genes[j];
				count1++;
			}
			average_exp.add(temp/count1);
			average_tot_reads_exp.add(temp1/count1);
			average_uniqueHits_exp.add(temp2/count1);
			average_rawReads_exp.add(temp3/count1);
			average_capReads_exp.add(temp4/count1);
			average_unique_hits_genes_exp.add(temp5/count1);
			for(int j=0;j<mapped_reads_reps_exp.size();j++)
			{
				X=X+Math.pow(((Double)TotalReads_per_gene_weighted[j].get(i)-average_exp.get(i)),2);//reads
				Y=Y+Math.pow(((Double)uniqueHits_per_gene.get(j).get(i)-average_uniqueHits_exp.get(i)),2);//insertions
			}
			sd_exp.add(Math.sqrt((X/((double)mapped_reads_reps_exp.size()-1.0))));
			sd_ins_exp.add(Math.sqrt((Y/((double)mapped_reads_reps_exp.size()-1.0))));
			count1=0.0;X=0.0;temp=0.0;temp1=0.0;temp3=0.0;temp4=0.0;temp2=0.0;temp5=0.0;
			for(int j=0;j<mapped_reads_reps_ctl.size();j++)
			{
				temp+=(Double)TotalReads_per_gene_weighted_ctl[j].get(i);
				temp1+=(double)total_reads_ctl[j];
				temp2+=(Double)uniqueHits_per_gene_ctl.get(j).get(i);
				temp3+=(Double)rawReads_ctl.get(j).get(i);
				temp4+=(Double)TotalReads_per_gene_ctl.get(j).get(i);
				temp5+=(double)total_uniqueHits_genes_ctl[j];
				count1++;
			}
			average_ctl.add(temp/count1);
			average_tot_reads_ctl.add(temp1/count1);
			average_uniqueHits_ctl.add(temp2/count1);
			average_rawReads_ctl.add(temp3/count1);
			average_capReads_ctl.add(temp4/count1);
			average_unique_hits_genes_ctl.add(temp5/count1);
			for(int j=0;j<mapped_reads_reps_ctl.size();j++)
			{
				X=X+Math.pow(((Double)TotalReads_per_gene_weighted_ctl[j].get(i)-average_ctl.get(i)),2);//reads
				Y=Y+Math.pow(((Double)uniqueHits_per_gene_ctl.get(j).get(i)-average_uniqueHits_ctl.get(i)),2);//insertions
			}
			sd_ctl.add(Math.sqrt((X/((double)mapped_reads_reps_ctl.size()-1.0))));
			sd_ins_ctl.add(Math.sqrt((Y/((double)mapped_reads_reps_ctl.size()-1.0))));
		}

//Code to generate t-statistic if sufficient number of samples provided using insertions
		//Calculated t-statistics for Unequal sample sizes, equal variance
		//first calculate grand sd...
		ArrayList<Double>  t_stat = new ArrayList<Double>();
		ArrayList<Double>  grand_sd = new ArrayList<Double>();
		ArrayList<Double>  pvalues = new ArrayList<Double>();
		ArrayList<Double> pvalues2 = new ArrayList<Double>();
		ArrayList<Double> adj_pvalues = new ArrayList<Double>();
		if(mapped_reads_reps_exp.size()>=2 && mapped_reads_reps_ctl.size()>=2)//ie if not at least 2 controls and treatments ignore significance testing
		{
			
			for(int i = 0;i<geneID.size();i++)//calculate grand standard deviation
			{
				grand_sd.add(Math.pow((((mapped_reads_reps_exp.size()-1)*Math.pow(sd_ins_exp.get(i),2)+(mapped_reads_reps_ctl.size()-1)*Math.pow(sd_ins_ctl.get(i),2))/(mapped_reads_reps_exp.size()+mapped_reads_reps_ctl.size()-2)),0.5));
			}
			for(int i = 0;i<geneID.size();i++)//calculate T stat
			{
          			Double temp=(average_uniqueHits_exp.get(i)-average_uniqueHits_ctl.get(i))/(grand_sd.get(i)*Math.pow(((1.0/(double)mapped_reads_reps_exp.size())+(1.0/(double)mapped_reads_reps_ctl.size())),0.5));
				t_stat.add(temp);
			}
			//determine significance
			StudentsT t = new StudentsT((double)(mapped_reads_reps_exp.size()+mapped_reads_reps_ctl.size()-2));
			for(int i = 0;i<geneID.size();i++)
			{
				if(new Double(t_stat.get(i)).isNaN())
				{
					pvalues.add(Double.NaN);
					pvalues2.add(Double.NaN);
				}
				else if(t_stat.get(i)<0.0)
				{
					pvalues.add(t.cdf(t_stat.get(i)));
					pvalues2.add(t.cdf(t_stat.get(i)));
				}
				else
				{
					pvalues.add(1.0-t.cdf(t_stat.get(i)));
					pvalues2.add(1.0-t.cdf(t_stat.get(i)));
				}
			}
			//correct for multiple testing (BH)
			Collections.sort(pvalues);
			double counter1=1.0;
			for(int i = 0;i<pvalues.size();i++)
			{
				adj_pvalues.add(Math.min(1.0,(((double)pvalues.size()/counter1)*(double)pvalues.get(i))));
				counter1++;
			}
		}


//Calculating ratio...
		int test_ins_count=0,test_read_count=0;
		ArrayList<Double>  ratio = new ArrayList<Double>();//read ratio
		ArrayList<Double>  ratio_ins = new ArrayList<Double>();//insertion ratio
		for(int i = 0;i<geneID.size();i++)
		{
			if(average_ctl.get(i)==0.0&&average_exp.get(i)==0.0)
			{
				ratio.add(1.0);
				test_read_count++;
			}
			else if(average_ctl.get(i)==1.0&&average_exp.get(i)==1.0)
			{
				ratio.add(1.0);
				test_read_count++;
			}
			else if(average_ctl.get(i)==0.0)
				ratio.add((average_exp.get(i)/average_tot_reads_exp.get(i))/(1.0/average_tot_reads_ctl.get(i)));
			else if(average_exp.get(i)==0.0)
				ratio.add((1.0/average_tot_reads_exp.get(i))/(average_ctl.get(i)/average_tot_reads_ctl.get(i)));
			else
				ratio.add((average_exp.get(i)/average_tot_reads_exp.get(i))/(average_ctl.get(i)/average_tot_reads_ctl.get(i)));
		}

		for(int i = 0;i<geneID.size();i++)
		{
			if(average_uniqueHits_exp.get(i)==0.0&&average_uniqueHits_ctl.get(i)==0.0)
			{
				ratio_ins.add(1.0);
				test_ins_count++;
			}
			else if(average_uniqueHits_exp.get(i)==1.0&&average_uniqueHits_ctl.get(i)==1.0)
				ratio_ins.add(1.0);
			else if(average_uniqueHits_ctl.get(i)==0.0)
				ratio_ins.add((average_uniqueHits_exp.get(i)/average_unique_hits_genes_exp.get(i))/(0.5/average_unique_hits_genes_ctl.get(i)));
			else if(average_uniqueHits_exp.get(i)==0.0)
				ratio_ins.add((0.5/average_unique_hits_genes_exp.get(i))/(average_uniqueHits_ctl.get(i)/average_unique_hits_genes_ctl.get(i)));
			else
				ratio_ins.add((average_uniqueHits_exp.get(i)/average_unique_hits_genes_exp.get(i))/(average_uniqueHits_ctl.get(i)/average_unique_hits_genes_ctl.get(i)));
		}


//Calculating significance using difference in proportions of insertions
		ArrayList<Object> prop_results=TnSeq_analysis_twoSample.proportions(average_uniqueHits_ctl,average_uniqueHits_exp,average_unique_hits_genes_exp, average_unique_hits_genes_ctl,geneID,test_ins_count);
		
		ArrayList<Double>  prop_insertion_pvalue = (ArrayList<Double>)prop_results.get(0);
		ArrayList<Double>  prop_insertion_pvalue2 = (ArrayList<Double>)prop_results.get(1);
		ArrayList<Double>  adj_prop_insertion_pvalue = (ArrayList<Double>)prop_results.get(2);

//Calculating significance using difference in proportions of reads
		prop_results=TnSeq_analysis_twoSample.proportions(average_ctl,average_exp,average_tot_reads_exp, average_tot_reads_ctl,geneID,test_read_count);

		ArrayList<Double>  prop_reads_pvalue = (ArrayList<Double>)prop_results.get(0);
		ArrayList<Double>  prop_reads_pvalue2 = (ArrayList<Double>)prop_results.get(1);
		ArrayList<Double>  adj_prop_reads_pvalue = (ArrayList<Double>)prop_results.get(2);


//Calculating significance for insertions using Fisher's exact test
		prop_results=TnSeq_analysis_twoSample.fisher(average_uniqueHits_ctl,average_uniqueHits_exp,average_unique_hits_genes_exp, average_unique_hits_genes_ctl,geneID,test_ins_count);
		
		ArrayList<Double>  fisher_insertion_pvalue = (ArrayList<Double>)prop_results.get(0);
		ArrayList<Double>  fisher_insertion_pvalue2 = (ArrayList<Double>)prop_results.get(1);
		ArrayList<Double>  adj_fisher_insertion_pvalue = (ArrayList<Double>)prop_results.get(2);

//Calculating significance for reads using Fisher's exact test
		prop_results=TnSeq_analysis_twoSample.fisher(average_ctl,average_exp,average_tot_reads_exp, average_tot_reads_ctl,geneID,test_read_count);

		ArrayList<Double>  fisher_reads_pvalue = (ArrayList<Double>)prop_results.get(0);
		ArrayList<Double>  fisher_reads_pvalue2 = (ArrayList<Double>)prop_results.get(1);
		ArrayList<Double>  adj_fisher_reads_pvalue = (ArrayList<Double>)prop_results.get(2);


//Print out results
		FileOutputStream out = new FileOutputStream("Conditional_essentiality.txt");
		BufferedWriter buf = new BufferedWriter(new OutputStreamWriter(out));

		int count=0;
		for(int i=0;i<geneID.size();i++)
		{
			if(result_format.equalsIgnoreCase("Long"))
			{
				count++;
				if(mapped_reads_reps_exp.size()>=2 && mapped_reads_reps_ctl.size()>=2)
				{
					if(count==1)
                    				buf.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Raw Reads(treatment)\tAve. Raw Read (control)\tAve. Capped_reads(treatment)\tAve. Capped reads (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (Reads)\tpvalue (proportions_insertions)\tAdj. pvalue (proportions_insertions)\tpvalue (proportions_reads)\tAdj. pvalue (proportions_reads)\tpvalue (Fisher_insertions)\tAdj. pvalue (Fisher_insertions)\tpvalue (Fisher_reads)\tAdj. pvalue (Fisher_reads)\tpvalue (t-test)\tAdj. pvalue (t-test)\n");
					buf.write(geneID.get(i)+"\t"+annotation.get(i) +"\t"+average_uniqueHits_exp.get(i)+"\t"+average_uniqueHits_ctl.get(i) +"\t"+average_rawReads_exp.get(i)+"\t"+average_rawReads_ctl.get(i) +"\t"+average_capReads_exp.get(i)+"\t"+average_capReads_ctl.get(i)+"\t"+average_exp.get(i)+"\t"+average_ctl.get(i)+"\t"+ratio_ins.get(i)+"\t"+(Math.log(ratio_ins.get(i))/Math.log(2))+"\t"+ratio.get(i)+"\t"+(Math.log(ratio.get(i))/Math.log(2))+"\t"+prop_insertion_pvalue2.get(i)+"\t"+ adj_prop_insertion_pvalue.get(prop_insertion_pvalue.indexOf(prop_insertion_pvalue2.get(i)))+"\t"+prop_reads_pvalue2.get(i)+"\t"+ adj_prop_reads_pvalue.get(prop_reads_pvalue.indexOf(prop_reads_pvalue2.get(i)))+"\t"+fisher_insertion_pvalue2.get(i)+"\t"+ adj_fisher_insertion_pvalue.get(fisher_insertion_pvalue.indexOf(fisher_insertion_pvalue2.get(i)))+"\t"+fisher_reads_pvalue2.get(i)+"\t"+ adj_fisher_reads_pvalue.get(fisher_reads_pvalue.indexOf(fisher_reads_pvalue2.get(i)))+"\t"+pvalues2.get(i)+"\t"+ adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\n");
				}
				else
				{
					if(count==1)
		            			buf.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Raw Reads(treatment)\tAve. Raw Read (control)\tAve. Capped reads(treatment)\tAve. Capped reads (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (reads)\tpvalue (proportions_insertions)\tAdj. pvalue (proportions_insertions)\tpvalue (proportions_reads)\tAdj. pvalue (proportions_reads)\tpvalue (Fisher_insertions)\tAdj. pvalue (Fisher_insertions)\tpvalue (Fisher_reads)\tAdj. pvalue (Fisher_reads)\n");
					buf.write(geneID.get(i)+"\t"+annotation.get(i) +"\t"+average_uniqueHits_exp.get(i)+"\t"+average_uniqueHits_ctl.get(i) +"\t"+average_rawReads_exp.get(i)+"\t"+average_rawReads_ctl.get(i) +"\t"+average_capReads_exp.get(i)+"\t"+average_capReads_ctl.get(i)+"\t"+average_exp.get(i)+"\t"+average_ctl.get(i)+"\t"+ratio_ins.get(i)+"\t"+(Math.log(ratio_ins.get(i))/Math.log(2))+"\t"+ratio.get(i)+"\t"+(Math.log(ratio.get(i))/Math.log(2))+"\t"+prop_insertion_pvalue2.get(i)+"\t"+ adj_prop_insertion_pvalue.get(prop_insertion_pvalue.indexOf(prop_insertion_pvalue2.get(i)))+"\t"+prop_reads_pvalue2.get(i)+"\t"+ adj_prop_reads_pvalue.get(prop_reads_pvalue.indexOf(prop_reads_pvalue2.get(i)))+"\t"+fisher_insertion_pvalue2.get(i)+"\t"+ adj_fisher_insertion_pvalue.get(fisher_insertion_pvalue.indexOf(fisher_insertion_pvalue2.get(i)))+"\t"+fisher_reads_pvalue2.get(i)+"\t"+ adj_fisher_reads_pvalue.get(fisher_reads_pvalue.indexOf(fisher_reads_pvalue2.get(i)))+"\n");
				}
			}
			else
			{
				count++;
				if(mapped_reads_reps_exp.size()>=2 && mapped_reads_reps_ctl.size()>=2)
				{
					if(count==1)
                    				buf.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (reads)\tAdj. pvalue (proportions_insertions)\tAdj. pvalue (proportions_reads)\tAdj. pvalue (Fisher_insertions)\tAdj. pvalue (Fisher_reads)\tAdj. pvalue (t-test)\n");
					buf.write(geneID.get(i)+"\t"+annotation.get(i) +"\t"+average_uniqueHits_exp.get(i)+"\t"+average_uniqueHits_ctl.get(i) +"\t"+average_exp.get(i)+"\t"+average_ctl.get(i)+"\t"+ratio_ins.get(i)+"\t"+(Math.log(ratio_ins.get(i))/Math.log(2))+"\t"+ratio.get(i)+"\t"+(Math.log(ratio.get(i))/Math.log(2))+"\t"+ adj_prop_insertion_pvalue.get(prop_insertion_pvalue.indexOf(prop_insertion_pvalue2.get(i)))+"\t"+ adj_prop_reads_pvalue.get(prop_reads_pvalue.indexOf(prop_reads_pvalue2.get(i)))+"\t"+ adj_fisher_insertion_pvalue.get(fisher_insertion_pvalue.indexOf(fisher_insertion_pvalue2.get(i)))+"\t"+ adj_fisher_reads_pvalue.get(fisher_reads_pvalue.indexOf(fisher_reads_pvalue2.get(i)))+"\t"+ adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\n");
				}
				else
				{
					if(count==1)
		            			buf.write("Gene ID\tAnnotation\tAveUnique hits (treatment)\tAve. Unique hits (control)\tAve. Weighted reads(treatment)\tAve. Weighted reads (control)\tRatio_Insertions (Treatment/control)\tLog-fold Change (Insertions)\tRatio_reads (Treatment/control)\tLog-fold Change (reads)\tAdj. pvalue (proportions_insertions)\tAdj. pvalue (proportions_reads)\tAdj. pvalue (Fisher_insertions)\tAdj. pvalue (Fisher_reads)\n");
					buf.write(geneID.get(i)+"\t"+annotation.get(i) +"\t"+average_uniqueHits_exp.get(i)+"\t"+average_uniqueHits_ctl.get(i) +"\t"+average_exp.get(i)+"\t"+average_ctl.get(i)+"\t"+ratio_ins.get(i)+"\t"+(Math.log(ratio_ins.get(i))/Math.log(2))+"\t"+ratio.get(i)+"\t"+(Math.log(ratio.get(i))/Math.log(2))+"\t"+ adj_prop_insertion_pvalue.get(prop_insertion_pvalue.indexOf(prop_insertion_pvalue2.get(i)))+"\t"+ adj_prop_reads_pvalue.get(prop_reads_pvalue.indexOf(prop_reads_pvalue2.get(i)))+"\t"+ adj_fisher_insertion_pvalue.get(fisher_insertion_pvalue.indexOf(fisher_insertion_pvalue2.get(i)))+"\t"+ adj_fisher_reads_pvalue.get(fisher_reads_pvalue.indexOf(fisher_reads_pvalue2.get(i)))+"\n");
				}
			}
		}
		buf.flush();buf.close();out.close();
		System.out.println("Done...");
		s = Calendar.getInstance();
		System.out.println(s.getTime());
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}

///////////////////METHODS START///////////////////////////////////////////////

//method to process aligned reads file for uniques hits and reads
private static ArrayList<Object> process_mapped_reads(BufferedReader b2,ArrayList<String> replicon_names,ArrayList<String> replicon_seqs,String mapped_reads_exp,String min_hits,String format) throws Exception
{
	int totalReads=0;
	String lines="";
	Pattern p1 = Pattern.compile("\t");
	ArrayList<Object> returnValue = new ArrayList<Object>();
	ArrayList<String> unique_locations = new ArrayList<String>();
	ArrayList<Integer> hits_per_unique_location = new ArrayList<Integer>();
	Hashtable<String, String> unique_locations2 = new Hashtable<String, String>(100,0.7f);
	if(format.equalsIgnoreCase("Bowtie"))
	{
		while ((lines = b2.readLine())!=null)
		{	
			String[] result = p1.split(lines);
			if(lines.charAt(0)!='@' && !result[2].trim().equals("*"))
			{
			totalReads++;//total number of reads
			if(!replicon_names.contains(result[2].trim()))
			{
				System.out.println("Names of replicons (either chromosomes or plasmids) in the sequence file and mapped reads file to do no match!!\n\n"+ result[2].trim() + " in mapped reads file but not in sequence file.\n\n"+"Exiting run...");
				return returnValue;
			}
			if(!unique_locations2.containsKey(result[2].trim()+"\t"+result[3].trim()))//bowtie
			{
				unique_locations.add(result[2].trim()+"\t"+result[3].trim());//total number of unique locations hit
				unique_locations2.put(result[2].trim()+"\t"+result[3].trim(),result[2].trim()+"\t"+result[3].trim());//total number of unique locations hit
				hits_per_unique_location.add(1);//total number of hits per unique location
			}
			else
			{
				hits_per_unique_location.set(unique_locations.indexOf(result[2].trim()+"\t"+result[3].trim()),(hits_per_unique_location.get(unique_locations.indexOf(result[2].trim()+"\t"+result[3].trim()))+1));
			}
			}
		}
	}
	else if(format.equalsIgnoreCase("SOAP"))
	{
		while ((lines = b2.readLine())!=null)
		{
			if(lines.charAt(0)!='@')
			{
			totalReads++;//total number of reads
			String[] result = p1.split(lines);
			if(!replicon_names.contains(result[7].trim()))
			{
				System.out.println("Names of replicons (either chromosomes or plasmids) in the sequence file and mapped reads file to do no match!!\n\n"+ result[7].trim() + " in mapped reads file but not in sequence file.\n\n"+"Exiting run...");
				return returnValue;
			}
			if(!unique_locations2.containsKey(result[7].trim()+"\t"+result[8].trim()))//bowtie
			{
				unique_locations.add(result[7].trim()+"\t"+result[8].trim());//total number of unique locations hit
				unique_locations2.put(result[7].trim()+"\t"+result[8].trim(),result[7].trim()+"\t"+result[8].trim());//total number of unique locations hit
				hits_per_unique_location.add(1);//total number of hits per unique location
			}
			else
			{
				hits_per_unique_location.set(unique_locations.indexOf(result[7].trim()+"\t"+result[8].trim()),(hits_per_unique_location.get(unique_locations.indexOf(result[7].trim()+"\t"+result[8].trim()))+1));
			}
			}
		}
	}
	else if(format.equalsIgnoreCase("Eland"))
	{
		while ((lines = b2.readLine())!=null)
		{
			if(lines.charAt(0)!='@')
			{
			totalReads++;//total number of reads
			String[] result = p1.split(lines);
			if(!replicon_names.contains(result[6].trim()))
			{
				System.out.println("Names of replicons (either chromosomes or plasmids) in the sequence file and mapped reads file to do no match!!\n\n"+ result[6].trim() + " in mapped reads file but not in sequence file.\n\n"+"Exiting run...");
				return returnValue;
			}
			if(!unique_locations2.containsKey(result[6].trim()+"\t"+result[7].trim()))//bowtie
			{
				unique_locations.add(result[6].trim()+"\t"+result[7].trim());//total number of unique locations hit
				unique_locations2.put(result[6].trim()+"\t"+result[7].trim(),result[6].trim()+"\t"+result[7].trim());//total number of unique locations hit
				hits_per_unique_location.add(1);//total number of hits per unique location
			}
			else
			{
				hits_per_unique_location.set(unique_locations.indexOf(result[6].trim()+"\t"+result[7].trim()),(hits_per_unique_location.get(unique_locations.indexOf(result[6].trim()+"\t"+result[7].trim()))+1));
			}
			}
		}
	}
	else
	{
		System.out.println("Unrecognized mapped file format : "+format);
		return returnValue;
	}
	
	System.out.println("Finished processing aligned reads ("+mapped_reads_exp+")...\n");
	System.out.println("Total number of mapped reads : "+totalReads+"\n");
	returnValue.add(unique_locations);//0
	returnValue.add(hits_per_unique_location);//1
	returnValue.add(unique_locations2);//2
	returnValue.add(totalReads);//3			
		
	//writing wig file with unnormalized reads
	System.out.println("Writing results: "+mapped_reads_exp.trim().replace(".txt","").replace(".mapped","")+"_Reads_per_uniqueLocation.wig\n");
	FileOutputStream out = new FileOutputStream(mapped_reads_exp.trim().replace(".txt","").replace(".mapped","")+"_Reads_per_uniqueLocation.wig");
	BufferedWriter buf = new BufferedWriter(new OutputStreamWriter(out));
	buf.write("track type=wiggle_0 autoScale=on name=\"TnSeq track\" description=\"TnSeq\"\n");
	for(int i = 0;i<replicon_names.size();i++)//for each replicon
	{
		buf.write("variableStep  chrom="+replicon_names.get(i)+"  span=1\n");
		for(int j = 0;j<replicon_seqs.get(i).length();j++)//for each base
		{
			if(unique_locations2.containsKey(replicon_names.get(i)+"\t"+j))
			{
				if(hits_per_unique_location.get(unique_locations.indexOf(replicon_names.get(i)+"\t"+j))>=Integer.parseInt(min_hits))
					buf.write(j+"\t"+hits_per_unique_location.get(unique_locations.indexOf(replicon_names.get(i)+"\t"+j))+"\n");
			}
		}
	}
	buf.flush();out.close();buf.close();
	return returnValue;
}

//method to identify hits on a per gene bases
private static ArrayList<Object> gene_level_analysis(String min_hits, String clippings,String mapped_reads_exp,ArrayList<Integer> end,ArrayList<Integer> start,ArrayList<String> unique_locations,ArrayList<Integer> hits_per_unique_location,ArrayList<String> replicon,ArrayList<String> geneID,ArrayList<String> annotation) throws Exception
{
	System.out.println("Commencing gene-based analysis...\n");
	Pattern p1 = Pattern.compile("\t");
	ArrayList<Object> returnValue = new ArrayList<Object>();
	ArrayList<Double> uniqueHits_per_gene = new ArrayList<Double>();
	ArrayList<Double> norm_uniqueHits_per_gene = new ArrayList<Double>();
	ArrayList<Double> TotalReads_per_gene = new ArrayList<Double>();
	ArrayList<Integer> GeneSize = new ArrayList<Integer>();
	int uniqueLocation_hits_threshold=Integer.parseInt(min_hits);
	int total_uniqueHits_genes=0,total_reads_genes=0,max=0,total_uniqueHits_genome=0;
    	double total_norm_uniqueHits_per_gene=0.0;
	for(int i = 0;i<geneID.size();i++)//for each gene
	{	
		int geneSize = end.get(i)-start.get(i)+1;
		int clipping=Math.round(geneSize*Float.parseFloat(clippings)/100);
		int hits=0;
		int totalReads_per_gene=0;
		total_uniqueHits_genome=0;//reset each time only need final value
		for(int j = 0;j<unique_locations.size();j++)//for each unique hit
		{
			String[] result = p1.split(unique_locations.get(j));
			if(replicon.get(i).equals(result[0])&& Integer.parseInt(result[1])>=(start.get(i)+clipping) && Integer.parseInt(result[1])<=(end.get(i)-clipping) && hits_per_unique_location.get(j)>=uniqueLocation_hits_threshold)//check if unique hit within gene i
			{

				hits++;
				totalReads_per_gene+=hits_per_unique_location.get(j);
				total_uniqueHits_genes++;
				total_reads_genes+=hits_per_unique_location.get(j);
				if(hits_per_unique_location.get(j)>max) max = hits_per_unique_location.get(j);
			}
			if(hits_per_unique_location.get(j)>=uniqueLocation_hits_threshold)//check if hits in a given location is upto threshold
			{
				total_uniqueHits_genome++;
			}
		}
		uniqueHits_per_gene.add((double)hits);
		TotalReads_per_gene.add((double)totalReads_per_gene);
		GeneSize.add((geneSize-(clipping*2)));

	}
	System.out.println("Total number of unique hits ("+mapped_reads_exp+") : "+ total_uniqueHits_genome+"\n");
	System.out.println("Total number of unique hits within genes ("+mapped_reads_exp+") : "+total_uniqueHits_genes+"\n");
	System.out.println("Total number of reads within genes ("+mapped_reads_exp+") : "+total_reads_genes+"\n");
	for(int i = 0;i<geneID.size();i++)
	{
        	total_norm_uniqueHits_per_gene+=(double)uniqueHits_per_gene.get(i)/(double)GeneSize.get(i);
		norm_uniqueHits_per_gene.add((double)uniqueHits_per_gene.get(i)/(double)GeneSize.get(i));
	}
	
	returnValue.add(uniqueHits_per_gene);//0
	returnValue.add(norm_uniqueHits_per_gene);//1
	returnValue.add(TotalReads_per_gene);//2
	returnValue.add(GeneSize);//3
	returnValue.add(max);//4
	returnValue.add(total_reads_genes);//5
	returnValue.add(total_uniqueHits_genes);//6
    	returnValue.add(total_norm_uniqueHits_per_gene);//7
	return returnValue;
}
//Capping function
private static ArrayList<Double> capping(String capping,String weights,String clippings,String min_hits,int max,int total_reads,ArrayList<String> unique_locations,ArrayList<Integer> hits_per_unique_location,ArrayList<String> geneID,ArrayList<String> replicon,ArrayList<Integer> end,ArrayList<Integer> start)
{
	System.out.println("Capping reads...\n");
	int cap=0;
	double mean=0.0;
	double median=0.0;
	double mean2SD=0.0;
	Pattern p1 = Pattern.compile("\t");
	if(capping.equals("0")) cap = max+1;
	else
	{
		mean = (double)total_reads/(double)unique_locations.size();//average number of reads per unique location
		double X=0.0;
		for(int i = 0;i<hits_per_unique_location.size();i++)//for unique location
		{	
			X=X+Math.pow(((double)hits_per_unique_location.get(i)-mean),2);
		}
		double sdX = Math.sqrt((X/((double)unique_locations.size()-1.0)));//std. dev.
		ArrayList<Integer> hits_per_unique_location_2 = new ArrayList<Integer>(hits_per_unique_location);
		Collections.sort(hits_per_unique_location_2);
		median = hits_per_unique_location_2.get((int)hits_per_unique_location_2.size()/2);//median
		mean2SD = Math.round((mean + (sdX*2)));//mean + (sdX*2);//
		if(capping.equals("1")) cap = Math.round((float)mean2SD);
		else if(capping.equals("2")) cap = Math.round((float)mean);
		else if(capping.equals("3")) cap = Math.round((float)median);
	}
	System.out.println("Mean :"+mean);
	System.out.println("Median :"+median);
	System.out.println("Mean2SD :"+mean2SD);
	System.out.println("Max :"+max);
	System.out.println("Read cap set at : "+cap + " reads per unique location");

	ArrayList<Double> TotalReads_per_gene_capped = new ArrayList<Double>();
	for(int i = 0;i<geneID.size();i++)//for each gene
	{	
		int geneSize = end.get(i)-start.get(i)+1;
		int clipping=Math.round(geneSize*Float.parseFloat(clippings)/100);
		int uniqueLocation_hits_threshold=Integer.parseInt(min_hits);
		int totalReads_per_gene=0;

		for(int j = 0;j<unique_locations.size();j++)//for each unique hit
		{
			String[] result = p1.split(unique_locations.get(j));
			if(replicon.get(i).equals(result[0])&& Integer.parseInt(result[1])>=(start.get(i)+clipping) && Integer.parseInt(result[1])<=(end.get(i)-clipping) && hits_per_unique_location.get(j)>=uniqueLocation_hits_threshold)//check if unique hit within gene i
			{	
				if(hits_per_unique_location.get(j)>cap) totalReads_per_gene+=cap;
				else totalReads_per_gene+=hits_per_unique_location.get(j);			
			}
		}
		TotalReads_per_gene_capped.add((double)totalReads_per_gene);
	}
	return TotalReads_per_gene_capped;
}

//Calculating significance using difference in proportions
private static ArrayList<Object> proportions(ArrayList<Double>  average_uniqueHits_ctl, ArrayList<Double>  average_uniqueHits_exp, ArrayList<Double>  average_unique_hits_genes_exp, ArrayList<Double>  average_unique_hits_genes_ctl,ArrayList<String> geneID,int test_ins_count)
{
		ArrayList<Object> returnValue = new ArrayList<Object>();
		ArrayList<Double>  prop_insertion_pvalue = new ArrayList<Double>();
		ArrayList<Double>  prop_insertion_pvalue2 = new ArrayList<Double>();
		ArrayList<Double>  adj_prop_insertion_pvalue = new ArrayList<Double>();
		ArrayList<Double>  prop_test_stat_temp = new ArrayList<Double>();
		for(int i = 0;i<geneID.size();i++)//for each gene
		{
			double diff = (average_uniqueHits_exp.get(i)/average_unique_hits_genes_exp.get(i))-(average_uniqueHits_ctl.get(i)/average_unique_hits_genes_ctl.get(i));//difference in count proportions
			double cum_prop = (average_uniqueHits_exp.get(i)+average_uniqueHits_ctl.get(i))/(average_unique_hits_genes_exp.get(i)+average_unique_hits_genes_ctl.get(i));//cumulative proportion calculation...
			double prop_test_stat=diff/(Math.sqrt((cum_prop*(1-cum_prop))*((1/average_unique_hits_genes_exp.get(i))+(1/average_unique_hits_genes_ctl.get(i)))));
			prop_test_stat_temp.add(prop_test_stat);
			if(new Double(prop_test_stat).isNaN())
			{
				prop_insertion_pvalue.add(1.0);//two tailed...
				prop_insertion_pvalue2.add(1.0);//two tailed...
			}
			else if(prop_test_stat<0)
			{
				prop_insertion_pvalue.add(new Normal().cdf(prop_test_stat)*2);//two tailed...
				prop_insertion_pvalue2.add(new Normal().cdf(prop_test_stat)*2);//two tailed...
			}
			else
			{

				prop_insertion_pvalue.add((1.0-new Normal().cdf(prop_test_stat))*2);//two tailed...
				prop_insertion_pvalue2.add((1.0-new Normal().cdf(prop_test_stat))*2);//two tailed...
			}
		}
		//correct for multiple testing (BH)
		Collections.sort(prop_insertion_pvalue);
		double counter1=1.0;
		for(int i = 0;i<prop_insertion_pvalue.size();i++)
		{
			//adj_prop_insertion_pvalue.add(Math.min(1.0,(((double)prop_insertion_pvalue.size()/counter1)*(double)prop_insertion_pvalue.get(i))));
			adj_prop_insertion_pvalue.add(Math.min(1.0,(((double)(prop_insertion_pvalue.size()-test_ins_count)/counter1)*(double)prop_insertion_pvalue.get(i))));
			if(counter1<(prop_insertion_pvalue.size()-test_ins_count))			
				counter1++;
		}
	returnValue.add(prop_insertion_pvalue);//0
	returnValue.add(prop_insertion_pvalue2);//1
	returnValue.add(adj_prop_insertion_pvalue);//2
	return returnValue;
}

//Calculating significance using Fisher's exact test
private static ArrayList<Object> fisher(ArrayList<Double>  average_uniqueHits_ctl, ArrayList<Double>  average_uniqueHits_exp, ArrayList<Double>  average_unique_hits_genes_exp, ArrayList<Double>  average_unique_hits_genes_ctl,ArrayList<String> geneID,int test_ins_count)
{
		ArrayList<Object> returnValue = new ArrayList<Object>();
		ArrayList<Double>  fisher_insertion_pvalue = new ArrayList<Double>();
		ArrayList<Double>  fisher_insertion_pvalue2 = new ArrayList<Double>();
		ArrayList<Double>  adj_fisher_insertion_pvalue = new ArrayList<Double>();
		for(int i = 0;i<geneID.size();i++)//for each gene
		{
			if((int)(average_unique_hits_genes_exp.get(i)-average_uniqueHits_exp.get(i))<0 || (int)(average_unique_hits_genes_ctl.get(i)-average_uniqueHits_ctl.get(i))<0)
			{
				fisher_insertion_pvalue.add(Double.NaN);
				fisher_insertion_pvalue2.add(Double.NaN);
			}
			else
			{
				ContingencyTable2x2 table = new ContingencyTable2x2((int)average_uniqueHits_exp.get(i).doubleValue(),(int)average_uniqueHits_ctl.get(i).doubleValue(),(int)(average_unique_hits_genes_exp.get(i)-average_uniqueHits_exp.get(i)),(int)(average_unique_hits_genes_ctl.get(i)-average_uniqueHits_ctl.get(i)));
				if((int)average_uniqueHits_exp.get(i).doubleValue()==0 && (int)average_uniqueHits_ctl.get(i).doubleValue()!=0)//set min count to 1
					table = new ContingencyTable2x2(1,(int)average_uniqueHits_ctl.get(i).doubleValue(),(int)(average_unique_hits_genes_exp.get(i)-average_uniqueHits_exp.get(i)),(int)(average_unique_hits_genes_ctl.get(i)-average_uniqueHits_ctl.get(i)));
				else if((int)average_uniqueHits_exp.get(i).doubleValue()!=0 && (int)average_uniqueHits_ctl.get(i).doubleValue()==0)//set min count to 1
					table = new ContingencyTable2x2((int)average_uniqueHits_exp.get(i).doubleValue(),1,(int)(average_unique_hits_genes_exp.get(i)-average_uniqueHits_exp.get(i)),(int)(average_unique_hits_genes_ctl.get(i)-average_uniqueHits_ctl.get(i)));

				if((int)average_uniqueHits_exp.get(i).doubleValue()==0 && (int)average_uniqueHits_ctl.get(i).doubleValue()==0)
				{
					fisher_insertion_pvalue.add(1.0);
					fisher_insertion_pvalue2.add(1.0);
				}
				else
				{
					FishersExactTest fisher = new FishersExactTest(table,H1.NOT_EQUAL);
					fisher_insertion_pvalue.add(fisher.getApproxSP());
					fisher_insertion_pvalue2.add(fisher.getApproxSP());
				}
			}
		}
		//correct for multiple testing (BH)
		Collections.sort(fisher_insertion_pvalue);

		double counter1=1.0;
		for(int i = 0;i<fisher_insertion_pvalue.size();i++)
		{
			//adj_fisher_insertion_pvalue.add(Math.min(1.0,(((double)fisher_insertion_pvalue.size()/counter1)*(double)fisher_insertion_pvalue.get(i))));
			adj_fisher_insertion_pvalue.add(Math.min(1.0,(((double)(fisher_insertion_pvalue.size()-test_ins_count)/counter1)*(double)fisher_insertion_pvalue.get(i))));
			if(counter1<(fisher_insertion_pvalue.size()-test_ins_count))			
				counter1++;
		}
	returnValue.add(fisher_insertion_pvalue);//0
	returnValue.add(fisher_insertion_pvalue2);//1
	returnValue.add(adj_fisher_insertion_pvalue);//2
	return returnValue;
}
////////////////////////////////////METHODS END/////////////////////////////////////////////////////////
}
