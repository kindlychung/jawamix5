package mixedmodel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;

import myPlotLab.MyHeatMap;
import myPlotLab.MyHistogram;
import myPlotLab.MyManhattanPlot;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

public class EMMAX {
	
	public static void emmax_analysis(String genotype_hdf5_file, String phe_file, String kinship_file, String output_folder, int round,
			double p_value_after_correction, int min_sample_size, double maf_plot_threshold){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
				EMMA emma=new EMMA(phenotypeS.phenotypes[phe_index], genotype, genotype.read_in_kinship(kinship_file));				
				if(emma.sample_size>=min_sample_size){
					System.out.println("Running Phenotype:"+ phe_index+":"+emma.phenotype.phe_id);
					String outputfile =output_folder+"/EMMAX."+phe_index+"_"+emma.phenotype.phe_id+".top";
					long startTime = System.currentTimeMillis();
					run_emmax_or_stepwise_for_one_pheno_block(outputfile, p_value_after_correction,emma, maf_plot_threshold, genotype_hdf5_file, round);
					System.out.println("Time Consumed: "+ (System.currentTimeMillis() - startTime)/1000+" seconds.\n");	
				}else{
					System.out.println("Not running due to too small sample size:" +emma.sample_size);
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}	
	
	public static void emmax_analysis(String genotype_hdf5_file, String phe_file, String kinship_file, String output_folder, int round,
			double p_value_after_correction, int min_sample_size, int phe_index, double maf_plot_threshold){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			EMMA emma=new EMMA(phenotypeS.phenotypes[phe_index], genotype, genotype.read_in_kinship(kinship_file));				
			if(emma.sample_size>=min_sample_size){
				System.out.println("Running Phenotype:"+ phe_index+":"+emma.phenotype.phe_id);
				String outputfile =output_folder+"/EMMAX."+phe_index+"_"+emma.phenotype.phe_id+".top";
				long startTime = System.currentTimeMillis();
				run_emmax_or_stepwise_for_one_pheno_block(outputfile, p_value_after_correction,emma, maf_plot_threshold, genotype_hdf5_file, round);
				System.out.println("Time Consumed: "+ (System.currentTimeMillis() - startTime)/1000+" seconds.\n");	
			}else{
				System.out.println("Not running due to too small sample size:" +emma.sample_size);
			}					
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void run_emmax_or_stepwise_for_one_pheno_block(String single_marker_outfile, double pvalue_threshold, EMMA emma, 
			double maf_plot_threshold, String genotype_hdf5_file, int round){
		// run REML
		int[] indexes={}, chrs={};
		emma.REMLE(emma.phenotype.values, emma.cbind(chrs, indexes), null);	
		System.out.println("Sample size: "+emma.sample_size);
		System.out.println("remle_REML="+emma.remle_REML);
		System.out.println("remle_delta="+emma.remle_delta);
		System.out.println("remle_ve="+emma.remle_ve);
		System.out.println("remle_vg="+emma.remle_vg);
		System.out.println("remle_pseudo-heritability="+emma.heritability);		
		double[][] decomposed_array = emma.reml_decompositioned_array();				
		emma.phenotype.generate_new_Y_by_multiplying(decomposed_array); 
		double[] intsept=new double[emma.sample_size];
		for(int i=0;i<emma.sample_size;i++){
			for(int j=0;j<emma.sample_size;j++){
				intsept[i]+=decomposed_array[i][j];
			}
		}
		// p-value threshold
		double total_num_sites=0;
		for(int chr=0;chr<emma.genotype.num_chrs;chr++)total_num_sites+=emma.genotype.num_sites[chr];
		double corrected_threshold=pvalue_threshold/total_num_sites;			
		try{			
			BufferedWriter bw = new BufferedWriter(new FileWriter(single_marker_outfile));
			bw.write("#SampleSize="+emma.sample_size+ "; REML="+emma.remle_REML+"; pseudo-heritability="+emma.heritability
					+"; Transform="+emma.phenotype.transform_approach+"(P="+emma.phenotype.normality_pvaule+")"+
					"; Phenotype="+emma.phenotype.phe_id+"; GenotypeFile="+genotype_hdf5_file+"\n");
			bw.write("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count\n");
			System.out.println("start EMMAX");
			
			for(int chr=0;chr<emma.genotype.num_chrs;chr++){
				int snp=0;
				for (HDF5MDDataBlock<MDDoubleArray> block : emma.genotype.position_fast_blocks[chr]){
					double[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){
						double[] X_ori=new double[emma.sample_size];
						for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
							X_ori[sample_index] = data4thisblock[var_index][emma.indexes_with_phenotype_in_genotype[sample_index]];
						}
						boolean the_same=true;
						for(int sample_index=1; sample_index<emma.sample_size; sample_index++){
							if(Double.compare(X_ori[sample_index], X_ori[0])!=0)the_same=false;
						}				
						if(the_same)continue;
						double[][] Xs_after=new double[emma.sample_size][2];						
						for(int i=0;i<emma.sample_size;i++){
							Xs_after[i][0]=intsept[i];
							for(int j=0;j<emma.sample_size;j++){
								Xs_after[i][1]+=(decomposed_array[i][j]*X_ori[j]);
							}
						}
						double[][] result=EMMAX.reg2results_emmax(emma.phenotype.new_Y_4emmax, Xs_after);
						if(result!=null){
							result[3]=new double[1];
							result[3][0]=EMMAX.mafc(X_ori);
							String out_res=EMMAX.results2string(result, corrected_threshold);
							if(out_res!=null){
								bw.write((chr+1)+","+emma.genotype.locations[chr][snp+var_index]+","+out_res+"\n");
							}
						}						
					}snp+=data4thisblock.length;
				}
			}			
			bw.close();
			make_plot_one_phen(single_marker_outfile,maf_plot_threshold);
			if(round>=2){				 
				// set the initial variant
				AssociationResults candidates=new AssociationResults(single_marker_outfile, corrected_threshold, maf_plot_threshold);
				candidates.generating_indexes(emma.genotype);
				int[] best_markers_indexes=new int[round];
				int[] best_markers_chrs=new int[round];
				int[] best_markers_locs=new int[round];		
//				OLSMultipleLinearRegression best_reg=null;
				double[] var_explained=new double[round];
				int best_P_index=0;
				// find the best P-value
				double best_var=0;
				for(int var_index=0;var_index<candidates.num_of_var;var_index++){
					if(Double.compare(best_var, candidates.AdjustedR2[var_index])<0){
						best_var=candidates.AdjustedR2[var_index];
						best_P_index=var_index;
					}
				}
				HashSet<Integer> best_indexes_set =new HashSet<Integer>();
				best_markers_indexes[0]=best_P_index;
				var_explained[0]=candidates.AdjustedR2[best_P_index];
				best_indexes_set.add(best_P_index);
				best_markers_chrs[0]=candidates.chr[best_P_index]-1;
				best_markers_locs[0]=candidates.location[best_P_index];				
				for(int round_index=2; round_index<=round;round_index++){
					// prepare the array including existing best markers
					double[][] Xs_ori=new double[emma.sample_size][round_index+1];
					for(int k=1;k<round_index;k++){
						int chr=best_markers_chrs[k-1], location=best_markers_locs[k-1];
						double[] the_full_list=emma.genotype.load_one_variant_by_location(chr, location);
						for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
							Xs_ori[sample_index][k] = the_full_list[emma.indexes_with_phenotype_in_genotype[sample_index]];
						}
					}
					double[][] Xs_after=new double[emma.sample_size][round_index+1];						
					for(int i=0;i<emma.sample_size;i++){
						Xs_after[i][0]=intsept[i];
						for(int var=1;var<round_index;var++){
							for(int j=0;j<emma.sample_size;j++){
								Xs_after[i][var]+=(decomposed_array[i][j]*Xs_ori[j][var]);
							}
						}						
					}
					// add the new marker and do regression.
//					best_P=1;
					best_P_index=0;
					double best_vc=0;						
					System.out.println("Start Step-wise. Round="+round_index);
					for(int var_index=0;var_index<candidates.num_of_var;var_index++){
						if(best_indexes_set.contains(var_index))continue;
						int chr=candidates.chr[var_index]-1;	
						double[] the_full_list=emma.genotype.load_one_variant_by_index(chr, candidates.indexes_in_genotype[var_index]);
						for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
							Xs_ori[sample_index][round_index] = the_full_list[emma.indexes_with_phenotype_in_genotype[sample_index]];
						}
						for(int i=0;i<emma.sample_size;i++){							
							for(int j=0;j<emma.sample_size;j++){
								Xs_after[i][round_index]+=(decomposed_array[i][j]*Xs_ori[j][round_index]);
							}													
						}
						double[][] result=EMMAX.reg2results_emmax(emma.phenotype.new_Y_4emmax, Xs_after);						
//						double[][] result=EMMAX.reg2results_lm(reg1, emma.sample_size);
						if(result!=null && Double.compare(best_vc, result[2][0])<0){
							best_P_index=var_index;
							best_vc=result[2][0];
						}
					}							
					// fill the structure that records existing selections
					var_explained[round_index-1]=best_vc;
					best_indexes_set.add(best_P_index);
					best_markers_chrs[round_index-1]=candidates.chr[best_P_index]-1;
					best_markers_locs[round_index-1]=candidates.location[best_P_index];
				}
				// output the results 
				System.out.println("Start writing stepwise results.");
				BufferedWriter bw_stw=new BufferedWriter(new FileWriter(single_marker_outfile+".stepwise.csv"));
				bw_stw.write("#Stepwise regression using variants reported in"+single_marker_outfile+"\n");
				bw_stw.write("#Round,chr,location,AdjustedR2\n");
				for(int k=0;k<round;k++){
					bw_stw.write((k+1)+","+(best_markers_chrs[k]+1)+","+best_markers_locs[k]+","+var_explained[k]+"\n");
				}
//				bw_stw.write("\nRegressionParameters:\n");
//				double[] coef_last=best_reg.estimateRegressionParameters();
//				for(int i=0;i<round;i++)bw_stw.write(coef_last[i]+",");
//				bw_stw.write(coef_last[round]+"\n");
				bw_stw.close();
			
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * This function makes two plots immediately after calling 
	 * public static void run_emmax_for_one_pheno(String output_pvalue_file, double pvalue_threshold, EMMA emma, double maf_plot_threshold)
	 *	(1) Manhattan plot for p-value.
	 *	(2) Manhattan plot for R2.
	 *	  
	 */
	public static void make_plot_one_phen(String emmax_results_file, double maf_plot_threshold){	
		double max_logp=30;
		String phe_name=emmax_results_file.split("/")[emmax_results_file.split("/").length-1];
		try{
			if((new File(emmax_results_file)).exists()){
				AssociationResults mlr=new AssociationResults(emmax_results_file, 2, maf_plot_threshold);
				double[] breaks=LocalKinshipAnalyzer.generate_breaks(LocalKinshipAnalyzer.lengths_of_chrs);
				double[] locations=new double[mlr.location.length];
				double[] pvalues=new double[mlr.location.length];
				double[] var=new double[mlr.location.length];
				for(int i=0;i<locations.length;i++){
					locations[i]=mlr.location[i]+breaks[mlr.chr[i]-1];		
					if(Double.compare(mlr.pvalue[i],0)>0){
						pvalues[i]=-Math.log10(mlr.pvalue[i]);
					}else if(Double.compare(mlr.pvalue[i],0)==0){
						pvalues[i]=max_logp;
					}else{System.out.println("P-value<0: Something wrong...");}
					var[i]=mlr.AdjustedR2[i]*100.0;
				}	
				MyManhattanPlot plot=new MyManhattanPlot(phe_name, "Locations", "-log10(p-value)", mlr.chr,
						locations, pvalues, 2000, 400, emmax_results_file+".pvalue.png");
				plot=new MyManhattanPlot(phe_name, "Locations", "Variance Explained (%)", mlr.chr,
						locations, var, 2000, 400, emmax_results_file+".VarR2.png");
			}else{
				System.out.println("NOFILE: "+emmax_results_file);
			}
			
			
		}catch(Exception e){e.printStackTrace();}
	}
		
	public static double[][] reg2results_lm(OLSMultipleLinearRegression reg1, int sample_size){
		double[][] results=new double[4][];	
		try{
			results[0]=reg1.estimateRegressionParameters();
			results[1]=myMathLib.StatFuncs.multi_reg_pvalues(results[0], reg1.estimateRegressionParametersStandardErrors(), sample_size);
			results[2]=new double[2];
			results[2][0]=reg1.calculateAdjustedRSquared();
			results[2][1]=reg1.estimateRegressionStandardError();
		}catch(Exception e){
			//e.printStackTrace();
			return null;
		} 				
		return results;
	}
	
	public static double[][] reg2results_emmax(double[] Y, double[][] Xs){
		int sample_size=Y.length;
		int var_num=Xs[0].length;
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		reg1.setNoIntercept(true);
		reg1.newSampleData(Y, Xs);		
		double[][] results=new double[4][];
		try{
		results[0]=reg1.estimateRegressionParameters();
		results[1]=myMathLib.StatFuncs.multi_reg_pvalues(results[0], reg1.estimateRegressionParametersStandardErrors(), sample_size);
		results[2]=new double[2];		
		results[2][1]=reg1.estimateRegressionStandardError();
		
		OLSMultipleLinearRegression reg2=new OLSMultipleLinearRegression();	
		reg2.setNoIntercept(true);
		double[] new_Y=new double[sample_size];
		double[][] new_X=new double[sample_size][var_num-1];
		for(int i=0;i<sample_size;i++){
			new_Y[i]=Y[i]-results[0][0]*Xs[i][0];
			for(int k=0;k<var_num-1;k++){
				new_X[i][k]=Xs[i][k+1];
			}
		}
		reg2.newSampleData(new_Y, new_X);
		results[2][0]=reg2.calculateAdjustedRSquared();
		}catch(Exception e){
			//e.printStackTrace();
			return null;
		} 
		return results;
	}
	
	/*
	 * load a few particular (may not continues) variants, specified by indeexes.
	 */
	public static double[][] load_data_according_phenotype(int[] chr_indexes, int[] snp_indexes, EMMA emma){
		double[][] Xs = new double[snp_indexes.length][emma.sample_size];	
		for(int snp=0; snp<snp_indexes.length; snp++){
			boolean the_same=true;
			double[] the_data=emma.genotype.load_one_variant_by_index(chr_indexes[snp], snp_indexes[snp]);
			for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
				Xs[snp][sample_index] = the_data[emma.indexes_with_phenotype_in_genotype[sample_index]];
				if(Double.compare(Xs[snp][sample_index], Xs[snp][0])!=0)the_same=false;
			}				
			if(the_same)
				return null;
		}return Xs;	
	}
	
	/*
	 * load a region of continues variants, specified by location.
	 */
	public static double[][] load_data_according_phenotype(int chr, int start_location, int end_location, EMMA emma){
		double[][] ori_data=emma.genotype.load_variants_in_region(chr, start_location, end_location);
		double[][] Xs = new double[ori_data.length][emma.sample_size];	
		boolean[] the_sames=new boolean[ori_data.length];
		int not_same=ori_data.length;
		for(int snp=0; snp<ori_data.length; snp++){
			boolean the_same=true;
			double[] the_data=ori_data[snp];
			for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
				Xs[snp][sample_index] = the_data[emma.indexes_with_phenotype_in_genotype[sample_index]];
				if(Double.compare(Xs[snp][sample_index], Xs[snp][0])!=0)the_same=false;
			}				
			if(the_same){
				the_sames[snp]=true;
				not_same--;
			}
		}
		if(not_same==0)return null;
		double[][] informative_Xs=new double[not_same][];
		int index=0;
		for(int snp=0; snp<ori_data.length; snp++){
			if(!the_sames[snp]){
				informative_Xs[index]=Xs[snp];
				index++;
			}
		}		
		return informative_Xs;	
	}
	
	/*
	 * return maf count
	 */
	public static int mafc(double[] X){
		HashMap<Integer, Integer> counts=new HashMap<Integer, Integer>();
		for(int k=0;k<X.length;k++){ 
			int the_allele=(int)(X[k]+0.5);
			if(counts.containsKey(the_allele)){
				int new_count=counts.get(the_allele)+1;
				counts.put(the_allele, new_count);
			}else{counts.put(the_allele, 1);}
		}
		int[] counts_array=new int[counts.size()];
		int[] allele_array=new int[counts.size()];
		int index=0;
		for(int allele:counts.keySet()){
			allele_array[index]=allele;
			counts_array[index++]=counts.get(allele);
		}if(counts_array.length==2){
			return ((counts_array[0]>counts_array[1])?counts_array[1]:counts_array[0])*2;
		}else if(counts_array.length==3){
			int zeor=0, other=0;
//			System.out.println("THREE");
			for(int n=0;n<3;n++){
				if(allele_array[n]==0){zeor+=2*counts_array[n];}
				else if(allele_array[n]==1){zeor+=counts_array[n]; other+=counts_array[n];}
				else if(allele_array[n]==2){other+=2*counts_array[n];}
				else{return -1;}
			}
			return (zeor<other)?zeor:other;
		}else{
			return -1;
		}
	}
	
	/*
	 * return maf count and its allele
	 * results[0]= allele, [1]= count
	 */
	public static int[] mafc_with_allele(double[] X){
		int[] results=new int[2];
		HashMap<Integer, Integer> counts=new HashMap<Integer, Integer>();
		for(int k=0;k<X.length;k++){ 
			int the_allele=(int)(X[k]+0.5);
			if(counts.containsKey(the_allele)){
				int new_count=counts.get(the_allele)+1;
				counts.put(the_allele, new_count);
			}else{counts.put(the_allele, 1);}
		}
		int[] counts_array=new int[counts.size()];
		int[] allele_array=new int[counts.size()];
		int index=0;
		for(int allele:counts.keySet()){
			allele_array[index]=allele;
			counts_array[index++]=counts.get(allele);
		}if(counts_array.length==2){
			if(counts_array[0]>counts_array[1]){
				results[0]=allele_array[1];
				results[1]=counts_array[1]*2;
			}else{
				results[0]=allele_array[0];
				results[1]=counts_array[0]*2;
			}return results;
		}else if(counts_array.length==3){
			int zero=0, other=0;
//			System.out.println("THREE");
			for(int n=0;n<3;n++){
				if(allele_array[n]==0){zero+=2*counts_array[n];}
				else if(allele_array[n]==1){zero+=counts_array[n]; other+=counts_array[n];}
				else if(allele_array[n]==2){other+=2*counts_array[n];}
				else{return null;}
			}
			if (zero<other){
				results[0]=0;
				results[1]=zero;
			}else{
				results[0]=1;
				results[1]=other;
			}return results;
		}else{
			return null;
		}
	}
	
	/*
	 * Apply decomposition array on one SNP and add the transformed intsept 
	 */
	public static double[][] apply_array(double[] intsept, double[] X_ori, double[][] decomposed_array){
		if(X_ori.length!=decomposed_array.length)return null;
		int sample_size=X_ori.length;
		boolean the_same=true;
		for(int sample_index=1; sample_index<sample_size; sample_index++){
			if(Double.compare(X_ori[sample_index], X_ori[0])!=0)the_same=false;
		}				
		if(the_same)return null;
		double[][] Xs_after=new double[sample_size][2];						
		for(int i=0;i<sample_size;i++){
			Xs_after[i][0]=intsept[i];
			for(int j=0;j<sample_size;j++){
				Xs_after[i][1]+=(decomposed_array[i][j]*X_ori[j]);
			}
		}
		return Xs_after;
	}
	
	/*
	 * extract the needed info to print out
	 */
	
	public static String results2string(double[][] array_from_regression, double corrected_threshold){
		if(array_from_regression==null)
			return null;
		//"pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count\n"
		if(array_from_regression[1][1]<=corrected_threshold){
			String result=array_from_regression[1][1]+","+array_from_regression[2][0]+","+array_from_regression[0][1]+","+
						array_from_regression[2][1]+","+(int)array_from_regression[3][0];			
			return result;
		}else{
			return null;
		}
	}
	/**
	 * @param args
	 */
//	public static void main(String[] args) {
//		long start=System.currentTimeMillis();
//		String genotype_csv_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv";
//		String phe_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/phen_raw_092910.tsv";
//		String kinship_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv.kinship";
//		String output_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/check_programs/the_new/";
//		String hdf5_file=genotype_csv_file+".hdf5";
//		double p_value_before_correction=1000;
//		int min_sample_size=40;
//	//	VariantsDouble mydata=VariantsDouble.importCSV(genotype_csv_file, hdf5_file);
//		
//		double maf_plot_threshold=0.05;
//		
//		emmax_analysis(hdf5_file, phe_file, kinship_file, output_folder,p_value_before_correction,min_sample_size, maf_plot_threshold);
//		
//		
//		
////		start=System.currentTimeMillis();
////		VariantsDouble my_data_again=new VariantsDouble(hdf5_file);
////		
////		System.out.println((System.currentTimeMillis()-start)/1000);	
////		for(int i=0;i<5;i++){
////			for(int k=0;k<my_data_again.num_sites[i];k++){
////				double[] x=my_data_again.load_one_variant_by_index(i, k);
////			}
////			System.out.println("finished chr"+(i+1));
////		}
////		System.out.println((System.currentTimeMillis()-start)/1000);
////		my_data_again.output2csv(genotype_csv_file+".compare.csv");
//
//	}

}
