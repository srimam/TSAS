import bin.*;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to analyze TnSeq data. This class calls either one or two sample analysis depending on the input in the Parameters.txt file

usage: java TSAS
*/

public class TSAS
{
	public static void main (String[] args)
	{
		try
		{
			BufferedReader b = new BufferedReader(new FileReader("Parameters.txt"));

			String lines="",analysis="";
			Pattern p6 = Pattern.compile("=");
			
			while ((lines = b.readLine())!=null)
			{
				if(lines.length()>0 && lines.charAt(0)!='#')
				{
					String[] result = p6.split(lines);
					if(result[0].trim().equalsIgnoreCase("Analysis_type")) 
					{ 
						if(result.length>1) analysis = result[1].trim();
						else analysis = "";
						if(analysis.equals("2")) 
						{
							System.out.println("Analysis type : Two-sample analysis\n\n");
							TnSeq_analysis_twoSample.execute();
						}
						else
						{
							System.out.println("Analysis type : One-sample analysis\n\n");
							TnSeq_analysis_oneSample.execute();
						}
					}
				}
			}
			b.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
