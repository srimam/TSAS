package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import java.math.*;
import jsc.distributions.*;
/*
Code to analyze TnSeq data using only one sample (or averaging over replicates of one condition)

usage: java TnSeq_analysis_oneSample
*/

public class TnSeq_analysis_oneSample
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


//Reading parameters file
		FileInputStream f3 = new FileInputStream(new File("Parameters.txt"));
		BufferedReader b3 = new BufferedReader(new InputStreamReader(f3));


//Processing Parameters.txt file
		String lines="",genome_fasta="",gff="",min_hits="",clippings="",capping="",weights="",format="",result_format="";
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
				else if(result[0].trim().equalsIgnoreCase("Treatment")) 
				{ 
					String[] result1;
					if(result.length>1) result1 = p7.split(result[1].trim());
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
						System.out.println("Result format : Short");
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
			System.out.println("\n\n"+gff+ " (GFF file) not found... Exiting run!!");
			return;
		}

		//Mapped file format
		if(!(format.equalsIgnoreCase("Bowtie") || format.equalsIgnoreCase("SOAP") || format.equalsIgnoreCase("Eland")))
		{
			System.out.println("Unrecognized mapped file format ("+format+")... Exiting run!!");
			return;
		}

		//Aligned reads Treatment
		if(mapped_reads_reps_exp.size()==0)
		{
			System.out.println("Please provide a valid mapped reads file (treatment)... Exiting run!!");
			return;
		}
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			newFile = new File(mapped_reads_reps_exp.get(i));
			if(!newFile.exists())
			{
				System.out.println("\n\n"+mapped_reads_reps_exp.get(i)+ " (mapped reads files) not found... Exiting run!!");
				return;
			}
		}
		


//Parsing fasta file for genome sequence data...

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
		ArrayList<Integer> gene_length = new ArrayList<Integer>();
		lines = "";
		int tracker=0;
		while ((lines = b1.readLine())!=null)
		{
			String[] result = p1.split(lines);
			if(result.length>2 && result[2].equals("gene"))
			{
				gene_length.add(Integer.parseInt(result[4])-Integer.parseInt(result[3])+1);
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
		
		
//Parsing and processing alignment file...

		System.out.println("Processing aligned reads...\n");
		ArrayList<ArrayList> unique_locations=new ArrayList<ArrayList>();
		ArrayList<ArrayList> hits_per_unique_location=new ArrayList<ArrayList>();
		int[] total_reads=new int[mapped_reads_reps_exp.size()];
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			FileInputStream f4 = new FileInputStream(new File(mapped_reads_reps_exp.get(i)));
			BufferedReader b4 = new BufferedReader(new InputStreamReader(f4));
			ArrayList<Object> treatment_results =  TnSeq_analysis_oneSample.process_mapped_reads(b4,replicon_names,replicon_seqs,mapped_reads_reps_exp.get(i),min_hits,format);
			if(treatment_results.size()==0) return;
			else
			{
				unique_locations.add((ArrayList<String>)treatment_results.get(0));
				hits_per_unique_location.add((ArrayList<Integer>)treatment_results.get(1));
				total_reads[i] = (Integer)treatment_results.get(3);
			}
			f4.close();b4.close();
		}


//Identifying hits per gene...

		System.out.println("Analyzing hits per gene...\n");
		ArrayList<ArrayList> TotalReads_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> uniqueHits_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> norm_uniqueHits_per_gene=new ArrayList<ArrayList>();
		ArrayList<ArrayList> GeneSize=new ArrayList<ArrayList>();
		int[] max =new int[mapped_reads_reps_exp.size()];
		int[] total_uniqueHits_genes=new int[mapped_reads_reps_exp.size()];
		int[] total_reads_genes=new int[mapped_reads_reps_exp.size()];
		int[] total_uniqueHits_genome=new int[mapped_reads_reps_exp.size()];
        double[] total_norm_uniqueHits_per_gene=new double[mapped_reads_reps_exp.size()];
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			ArrayList<Object> treatment_results2 = TnSeq_analysis_oneSample.gene_level_analysis(min_hits,clippings,mapped_reads_reps_exp.get(i),end,start,unique_locations.get(i),hits_per_unique_location.get(i),replicon,geneID,annotation);

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
				total_uniqueHits_genome[i] = (Integer)treatment_results2.get(8);
			}
		}
		
//Averaging across replicates (if any)...

		double ave_unique_locations=0.0;
		for(int i=0;i<mapped_reads_reps_exp.size();i++)
		{
			ave_unique_locations+=total_uniqueHits_genome[i];
		}
		ave_unique_locations=ave_unique_locations/(double)mapped_reads_reps_exp.size();

		ArrayList<Double> ave_uniqueHits_per_gene=new ArrayList<Double>();
		ArrayList<Double> ave_TotalReads_per_gene=new ArrayList<Double>();
		for(int j = 0;j<geneID.size();j++)
		{
			double temp=0.0,temp2=0.0;
			for(int i=0;i<mapped_reads_reps_exp.size();i++)
			{
				temp+=(Double)uniqueHits_per_gene.get(i).get(j);
				temp2+=(Double)TotalReads_per_gene.get(i).get(j);
			}
			ave_uniqueHits_per_gene.add(temp/(double)mapped_reads_reps_exp.size());
			ave_TotalReads_per_gene.add(temp2/(double)mapped_reads_reps_exp.size());
		}

//Calculating statistical significance using binomial distribution and adjusting pvalues using BH approach...

		ArrayList<Double> pvalues = new ArrayList<Double>();//X<=x
		ArrayList<Double> pvalues2 = new ArrayList<Double>();//X<=x
		ArrayList<Double> adj_pvalues = new ArrayList<Double>();//X<=x
		ArrayList<Double> bonferroni_adj_pvalues = new ArrayList<Double>();//X<=x

		ArrayList<Double> pvalues3 = new ArrayList<Double>();//X>x
		ArrayList<Double> pvalues4 = new ArrayList<Double>();//X>x
		ArrayList<Double> adj_pvalues2 = new ArrayList<Double>();//X>x
		ArrayList<Double> bonferroni_adj_pvalues2 = new ArrayList<Double>();//X>x
		for(int i = 0;i<geneID.size();i++)
		{
			double temp_hold=TnSeq_analysis_oneSample.binomial_essentiality((Integer)GeneSize.get(0).get(i),(ave_unique_locations/(double)genomeSize),(int)ave_uniqueHits_per_gene.get(i).doubleValue());
			pvalues.add(Math.min(1.0,temp_hold));
			pvalues2.add(Math.min(1.0,temp_hold));

			pvalues3.add(1.0-Math.min(1.0,temp_hold));
			pvalues4.add(1.0-Math.min(1.0,temp_hold));
		}
		Collections.sort(pvalues);
		Collections.sort(pvalues3);
		double counter1=1;
		for(int i = 0;i<pvalues.size();i++)
		{
			adj_pvalues.add(Math.min(1.0,(((double)pvalues.size()/counter1)*(double)pvalues.get(i))));
			bonferroni_adj_pvalues.add(Math.min(1.0,(((double)pvalues.size())*(double)pvalues.get(i))));

			adj_pvalues2.add(Math.min(1.0,(((double)pvalues3.size()/counter1)*(double)pvalues3.get(i))));
			bonferroni_adj_pvalues2.add(Math.min(1.0,(((double)pvalues3.size())*(double)pvalues3.get(i))));
			counter1++;
		}

//Printing out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Essential_genes.txt"));
		int counter=0;
		for(int i = 0;i<geneID.size();i++)
		{
			if(result_format.equalsIgnoreCase("Long"))
			{
				counter++;
				if(counter==1)
					buf.write("Gene ID\tAnnotation\tGene length (bp)\tNo. of Unique hits\tNormalize Unique hits (hits/bp)\tTotal number of reads\tPvalue (Essential)\tAdj. Pvalue (Essential)\tFWER (Essential)\tPvalue (Improved fitness)\tAdj. Pvalue (Improved fitness)\tFWER (Improved fitness)\n");
				buf.write(geneID.get(i)+"\t"+annotation.get(i)+"\t"+ gene_length.get(i)+"\t"+ave_uniqueHits_per_gene.get(i).doubleValue()+"\t"+(ave_uniqueHits_per_gene.get(i)/(Integer)GeneSize.get(0).get(i)) +"\t"+ave_TotalReads_per_gene.get(i)+"\t"+pvalues2.get(i)+"\t"+ adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\t"+ bonferroni_adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\t"+pvalues4.get(i)+"\t"+ adj_pvalues2.get(pvalues3.indexOf(pvalues4.get(i)))+"\t"+ bonferroni_adj_pvalues2.get(pvalues3.indexOf(pvalues4.get(i)))+"\n");
			}
			else
			{
				counter++;
				if(counter==1)
					buf.write("Gene ID\tAnnotation\tGene length (bp)\tNo. of Unique hits\tNormalize Unique hits (hits/bp)\tTotal number of reads\tAdj. Pvalue (Essential)\tFWER (Essential)\tAdj. Pvalue (Improved fitness)\tFWER (Improved fitness)\n");
				buf.write(geneID.get(i)+"\t"+annotation.get(i)+"\t"+ gene_length.get(i)+"\t"+ave_uniqueHits_per_gene.get(i).doubleValue()+"\t"+(ave_uniqueHits_per_gene.get(i)/(Integer)GeneSize.get(0).get(i)) +"\t"+ave_TotalReads_per_gene.get(i)+"\t"+ adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\t"+ bonferroni_adj_pvalues.get(pvalues.indexOf(pvalues2.get(i)))+"\t"+ adj_pvalues2.get(pvalues3.indexOf(pvalues4.get(i)))+"\t"+ bonferroni_adj_pvalues2.get(pvalues3.indexOf(pvalues4.get(i)))+"\n");
			}
		}
		buf.flush();buf.close();

		System.out.println("Done...");
		s = Calendar.getInstance();
		System.out.println(s.getTime());
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}

/////////End of analysis....


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
			if(replicon.get(i).equals(result[0])&& Integer.parseInt(result[1])>=(start.get(i)+clipping) && Integer.parseInt(result[1])<=(end.get(i)-clipping) && hits_per_unique_location.get(j)>=uniqueLocation_hits_threshold)//chech if unique hit within gene i
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
	returnValue.add(total_uniqueHits_genome);//8
	return returnValue;
}

//method to calculate lower tail significance ie P[X<=x] using binomial distribution
private static double binomial_essentiality(int length_of_gene, double probability_of_hit_in_genome,int number_of_unique_hits_per_gene)
{
	int n = length_of_gene;//length of gene
	double p = probability_of_hit_in_genome;//probablity of hit... genome-wide number of insertions per bp
	int k = number_of_unique_hits_per_gene;//number of unique hits in gene...

	double pbinom= 0.0;
	double nfact= 1.0;
	double kfact= 1.0;
	if(k <=230)// this algorithm appears to break down after 230 insertions
	{
		for(int j=k;j>=0;j--)
		{
			for(int i=n;i>(n-j);i--)
			{
				nfact+=Math.log10(i);
			}
			for(int i=j;i>1;i--)
			{
				kfact+=Math.log10(i);
			}
			if(Math.pow(p,j)==0.0||Math.pow(1-p,n-j)==0.0){}
			else
				pbinom+=Math.pow(10,(nfact-kfact)+Math.log10(Math.pow(p,j))+Math.log10(Math.pow(1-p,n-j)));
			nfact= 1.0;
			kfact= 1.0;
		}
	}
	else//this function has lower precision but works ok at greater than 230 insertions
	{
		Binomial binom = new Binomial(length_of_gene,probability_of_hit_in_genome);
		pbinom = binom.cdf(number_of_unique_hits_per_gene);
	}
	return pbinom;
}
}
