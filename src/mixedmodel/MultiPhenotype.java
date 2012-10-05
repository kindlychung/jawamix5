package mixedmodel;

import java.io.BufferedReader;
import java.io.FileReader;

public class MultiPhenotype {
	
	public int num_of_pheno;
	public int sample_size_all_pheno;
//	String[] pheno_names;
	public Phenotype[] phenotypes;
	public String[] sample_ids_all_pheno;
	
	public MultiPhenotype(String phenotype_file){
		try{
			BufferedReader br = new BufferedReader(new FileReader(phenotype_file));
			String line=br.readLine(); // the header.
			String[] temp=line.split("\t");
			this.num_of_pheno=temp.length-1;
			String[] pheno_names=new String[this.num_of_pheno];
			this.phenotypes=new Phenotype[this.num_of_pheno];
			for(int i=0;i<this.num_of_pheno;i++){
				pheno_names[i]=temp[i+1];
			}		
			line=br.readLine();
			while(line!=null){
				this.sample_size_all_pheno++;
				line=br.readLine();
			}
			this.sample_ids_all_pheno=new String[this.sample_size_all_pheno];
			double[][] full_phenotype=new double[num_of_pheno][this.sample_size_all_pheno];
			br = new BufferedReader(new FileReader(phenotype_file));
			line=br.readLine();
			line=br.readLine();
			int sample_index=0;
			while(line!=null){
				temp=line.split("\t");
				this.sample_ids_all_pheno[sample_index]=temp[0];					
				for(int i=0;i<this.num_of_pheno;i++){
					if(temp[i+1].equals("NA")||temp[i+1].equals("NaN")||temp[i+1].equals("*")||temp[i+1].equals("N")){
						full_phenotype[i][sample_index]=Double.NaN;
					}else{
						full_phenotype[i][sample_index]=Double.parseDouble(temp[i+1]);
					}
				}
				sample_index++;
				line=br.readLine();
			}
			for(int k=0;k<this.num_of_pheno;k++){
				this.phenotypes[k]=new Phenotype(pheno_names[k], sample_ids_all_pheno, full_phenotype[k]);
			}
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	
}
