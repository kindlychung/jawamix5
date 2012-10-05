package mixedmodel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import myMathLib.StatFuncs;

import org.apache.commons.math3.stat.StatUtils;


public class Phenotype {
	public String phe_id;
	public String[] sample_ids;
	public double[] values;
	double[] new_Y_4emmax;
	String transform_approach="NoTransform";
	double normality_pvaule=-1;
	
	public Phenotype(String infile){
		ArrayList<String> ids=new ArrayList<String>();
		ArrayList<Double> values=new ArrayList<Double>();
		try{
			BufferedReader br =new BufferedReader(new FileReader(infile));
			String line=br.readLine();//header
			this.phe_id= line.split("\t")[1];
			line=br.readLine();
			while(line!=null){
				String temp[] = line.split("\t");
				ids.add(temp[0]);
				if(temp[1].equals("NA")||temp[1].equals("NaN")||temp[1].equals("N")||temp[1].equals("*")){
					values.add(Double.NaN);
				}else{
					values.add(Double.parseDouble(temp[1]));
				}				
				line=br.readLine();
			}
			this.sample_ids=new String[ids.size()];
			this.values=new double[ids.size()];
			for(int k=0;k<ids.size();k++){
				this.sample_ids[k]=ids.get(k);
				this.values[k]=values.get(k);
			}
		}catch(Exception e){e.printStackTrace();}
		remove_NaN();
	}
	
	public Phenotype(String phe_id, String[] sample_ids, double[] values){
		this.phe_id=phe_id;
		this.sample_ids=sample_ids;
		this.values=values;
		remove_NaN();
	}
	
	public void transform(){
		this.normality_pvaule=StatFuncs.normality_test_ks(this.values);
		if(normality_pvaule>0.5){
			return;
		}
		System.out.println("Phenotype may not follow normal distribution. K-S pvalue="+normality_pvaule);
		String[] potential_transforms={"\"log\"","\"exp\"","\"sqrt\"","\"sqr\""};
		double[] pvalues=new double[potential_transforms.length];
		double[][] new_values=new double[potential_transforms.length][this.values.length];
		double min=StatUtils.min(values);
		if(min<=0){
			for(int k=0;k<this.values.length;k++){
				new_values[0][k]=Math.log(values[k]-min+1);
				new_values[1][k]=Math.exp(values[k]);
				new_values[2][k]=Math.sqrt(values[k]-min);
				new_values[3][k]=(values[k]-min)*(values[k]-min);
			}
		}else{
			for(int k=0;k<this.values.length;k++){
				new_values[0][k]=Math.log(values[k]);
				new_values[1][k]=Math.exp(values[k]);
				new_values[2][k]=Math.sqrt(values[k]);
				new_values[3][k]=(values[k])*(values[k]);
			}
		}		
		for(int i=0;i<potential_transforms.length;i++){
			pvalues[i]=StatFuncs.normality_test_ks(new_values[i]);
			System.out.println("Trying "+potential_transforms[i]+" transformation. P-value="+pvalues[i]+".");
		}
		int best_trans=-1;
		double the_pvalue=normality_pvaule;
		for(int i=0;i<potential_transforms.length;i++){
			if(pvalues[i]>the_pvalue){
				best_trans=i;
				the_pvalue=pvalues[i];
			}
		}
		if(best_trans!=-1){
			System.out.println("Used "+potential_transforms[best_trans]+" transformation. New K-S P-value="+pvalues[best_trans]);
			this.transform_approach=potential_transforms[best_trans];
			this.values=new_values[best_trans];
			normality_pvaule=pvalues[best_trans];
		}
	}
	
	public void substract_pedigree_mean_for_NAM(){
		HashMap<String, Integer> counts=new HashMap<String, Integer>();
		HashMap<String, Double> sum=new HashMap<String, Double>();
		for(int i=0;i<this.sample_ids.length;i++){
			String pedigree_id=this.sample_ids[i].split("_")[0];
			myFileFunctions.FileFunc.add2_counts_hashmap(counts, pedigree_id, 1);
			myFileFunctions.FileFunc.add2_sum_hashmap(sum, pedigree_id, this.values[i]);
		}for(int i=0;i<this.sample_ids.length;i++){
			String pedigree_id=this.sample_ids[i].split("_")[0];
			this.values[i]=this.values[i]-sum.get(pedigree_id)/counts.get(pedigree_id);
		}
		System.out.println("Finished substracting pedigree mean for NAM: "+this.phe_id);
	}
	
	void remove_NaN(){
		int nan=0;
		for(int i=0;i<sample_ids.length;i++){
			if(Double.isNaN(this.values[i])){
				nan++;
			}
		}
		String[] new_ids=new String[sample_ids.length-nan];
		double[] new_values=new double[sample_ids.length-nan];
		int index=0;
		for(int i=0;i<sample_ids.length;i++){
			if(!Double.isNaN(this.values[i])){
				new_ids[index]=sample_ids[i];
				new_values[index]=values[i];
				index++;
			}
		}this.sample_ids=new_ids;
		this.values=new_values;
	}
	
	/*
	 * 	generate new phenotype array by applying decomposition array, 
	 * only samples indexed in sampel_pheno_index[] applied
	 */
	public void generate_new_Y_by_multiplying(double[][] decomposed_array){	
		this.new_Y_4emmax=new double[sample_ids.length];
		if(decomposed_array.length!=values.length){
			System.out.println("Decomposition array is inconsistent to phenotype array!");
			return;
		}
		for(int i=0;i<values.length;i++){
			for(int j=0;j<values.length;j++){
				this.new_Y_4emmax[i]=this.new_Y_4emmax[i]+decomposed_array[i][j]*values[j];
			}
		}
	}
}
