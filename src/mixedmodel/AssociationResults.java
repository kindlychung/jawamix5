package mixedmodel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;


public class AssociationResults {
	
	public static String[] supported_types={"emmax","local"};
	
	public String type;
	public int sample_size;
	public int num_of_var=0;
	public int[] chr;
	public int[] location;
	public int[] indexes_in_genotype;
	public double[] pvalue;
	public double[] AdjustedR2;
	public int[] MAF_count;
	
	/*
	 * Constructor: the input pvalue file should be a .csv file and it is fine if separated by "," or ", ";
	 * the first line is skipped, and the second line is the header.
	 */
	public AssociationResults(String pvaule_file, double pvalue_threshold){
		try{
			String sep=",";
			BufferedReader br = new BufferedReader(new FileReader(pvaule_file));
			String line=br.readLine();
			if(line.startsWith("#SampleSize="))
			this.sample_size=Integer.parseInt(line.split("; ")[0].split("=")[1]);
			line=br.readLine();
			if(line.equals("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count")){ // emmax file
				 this.type="emmax";
			}else if(line.equals("chr, position, p-value, variance_local, variance_global, h1_heritabilities")){ //local-kinship file
				this.type="local";
				sep=", ";
			}else{
				System.out.println("Header line is not correct!");
				return;
			}
			line=br.readLine();
			while(line!=null){
				if(Double.parseDouble(line.split(sep)[2])<pvalue_threshold)this.num_of_var++;
				line=br.readLine();
			}
			this.chr=new int[this.num_of_var];
			this.location=new int[this.num_of_var];
			this.pvalue=new double[this.num_of_var];
			this.AdjustedR2=new double[this.num_of_var];
			this.MAF_count=new int[this.num_of_var];
			br = new BufferedReader(new FileReader(pvaule_file));
			line=br.readLine();line=br.readLine();line=br.readLine();
			int p_index=0;
			while(line!=null){
				String temp[]=line.split(sep);
				double pvalue=Double.parseDouble(temp[2]);
				if(pvalue<pvalue_threshold){
					this.chr[p_index]=Integer.parseInt(temp[0]);
					this.location[p_index]=Integer.parseInt(temp[1]);
					this.pvalue[p_index]=pvalue;
					this.AdjustedR2[p_index]=Double.parseDouble(temp[3]);
					if(this.type.equals("emmax"))this.MAF_count[p_index]=Integer.parseInt(temp[6]);
					p_index++;
				}
				line=br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public AssociationResults(String pvaule_file, double pvalue_threshold, double maf_threshold){
		try{
			String sep=",";
			BufferedReader br = new BufferedReader(new FileReader(pvaule_file));
			String line=br.readLine();
			if(line.startsWith("#SampleSize="))
			this.sample_size=Integer.parseInt(line.split("; ")[0].split("=")[1]);
			double maf_count_threshold=this.sample_size*2*maf_threshold;
			line=br.readLine();
			if(line.equals("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count")|| // emmax file
					line.equals("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count,Region_id,syn_loc")|| //Rare + Syn
					line.equals("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count,Region_id")){ // Rare file
				 this.type="emmax";
			}else if(line.equals("chr, position, p-value, variance_local, variance_global, h1_heritabilities")){ //local-kinship file
				this.type="local";
				sep=", ";
			}else{
				System.out.println("Header line is not correct!\n"+line);
				return;
			}
			line=br.readLine();
			while(line!=null){
				String temp[]=line.split(sep);
				double pvalue=Double.parseDouble(temp[2]);
				double maf_count= Double.parseDouble(temp[6]);
				if(pvalue<pvalue_threshold && maf_count>=maf_count_threshold){
					this.num_of_var++;
				}
				line=br.readLine();
			}
			this.chr=new int[this.num_of_var];
			this.location=new int[this.num_of_var];
			this.pvalue=new double[this.num_of_var];
			this.AdjustedR2=new double[this.num_of_var];
			this.MAF_count=new int[this.num_of_var];
			br = new BufferedReader(new FileReader(pvaule_file));
			line=br.readLine();line=br.readLine();line=br.readLine();
			int p_index=0;
			while(line!=null){
				String temp[]=line.split(sep);
				double pvalue=Double.parseDouble(temp[2]);
				if(pvalue<pvalue_threshold && Double.parseDouble(temp[6])>=maf_count_threshold){
					this.chr[p_index]=Integer.parseInt(temp[0]);
					this.location[p_index]=Integer.parseInt(temp[1]);
					this.pvalue[p_index]=pvalue;
					this.AdjustedR2[p_index]=Double.parseDouble(temp[3]);
					if(this.type.equals("emmax"))this.MAF_count[p_index]=Integer.parseInt(temp[6]);
					p_index++;
				}
				line=br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void generating_indexes(VariantsDouble genotype){
		this.indexes_in_genotype=new int[this.num_of_var];
		int[] num_var_per_char=new int[genotype.num_chrs];
		for(int var=0;var<this.num_of_var;var++){
			num_var_per_char[this.chr[var]-1]++;
		}
		for(int var=0;var<this.num_of_var;var++){
			int chr=this.chr[var]-1;
			int index=Arrays.binarySearch(genotype.locations[chr], this.location[var]);
			if(index<0){
				System.out.println("Error: Variant location not in the list. Failed to set up indexes!");
				return;
			}
			this.indexes_in_genotype[var]=index;
		}
		
	}
}
