package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;


import myMathLib.Test;
import myPlotLab.MyHeatMap;
import myPlotLab.MyHistogram;
import myPlotLab.MyManhattanPlot;


public class LocalKinshipAnalyzer {
	VariantsDouble variants;
	int tiling_win_size;
	int[][][] gene_coordinates; // chr x num_of_genes_per_chr x 2[start, end]
	String[][] gene_ids; 	// chr x num_of_genes_per_chr
	double[][] global_kinship;
	public static int[] lengths_of_chrs={30427671,19698289,23459830,18585056,26975502};
	
	/*
	 * Constructor.
	 * the "gene_file" actually can be any file specifying interested regions in the format:
	 * id\chr\tstart\itend\n
	 * 
	 * the chr has to be numbered as "1","2,"... instead of strings like "Chr1" etc.
	 */
	public LocalKinshipAnalyzer(String hdf5_file, int win_size, String gene_file){
		try{
			this.variants=new VariantsDouble(hdf5_file);
			this.tiling_win_size=win_size;
			if(gene_file!=null){
				BufferedReader br=new BufferedReader(new FileReader(gene_file));
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * This is just a template of doing anything that need to go through the tiling windows.
	 */
	public void go_through_tiling_windows(){
		try{
//			BufferedWriter bw=new BufferedWriter(new FileWriter(outfile));
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=tiling_win_size;
				while(start<last_position){
					double[][] data=variants.load_variants_in_region(chr, start, end);					
					for(int k=0;k<data.length;k++){
						//do something
					}
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}
			}//bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * calculate local-kinship 
	 * scale = the biggest difference between genotypes.
	 */
	public void calculating_kinship_tiling_windows(String output_folder, double scale){		
		System.out.println("Local_kinship for win-size="+tiling_win_size+".");
		try{
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){
					double[][] data=variants.load_variants_in_region(chr, start, end);
					String the_kin_file=output_folder+"kinship."+variants.sample_size+".chr"
							+(chr+1)+"."+start+"."+end;
					BufferedWriter bw=new BufferedWriter(new FileWriter(the_kin_file+".raw.ibs"));
					double[][] kinship=new double[variants.sample_size][variants.sample_size];
					for(int i=0;i<variants.sample_size;i++){
						kinship[i][i]=1;
						for(int j=i+1;j<variants.sample_size;j++){
							for(int k=0;k<data.length;k++){
								kinship[i][j]+=(scale-Math.abs(data[k][i]-data[k][j]));
							}
							kinship[i][j]=kinship[i][j]/scale/data.length;
							kinship[j][i]=kinship[i][j];
						}
					}
					for(int i=0;i<variants.sample_size;i++){
						for(int j=0;j<variants.sample_size-1;j++){
							bw.write(kinship[i][j]+",");
						}bw.write(kinship[i][variants.sample_size-1]+"\n");
					}
					bw.close();
					VariantsDouble.re_scale_kinship_matrix(the_kin_file+".raw.ibs", the_kin_file+".rescaled.ibs");
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void local_VO(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
			String mlr_output_file, double step){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			System.out.println(phenotype.phe_id);
			
			EMMA emma_global=new EMMA(phenotype, genotype, genotype.read_in_kinship(global_kinship_file));
//			emma_global.phenotype.transform();
			emma_global.REMLE_null();
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			BufferedWriter bw_dist=new BufferedWriter(new FileWriter(mlr_output_file+".dist.csv"));
			bw.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			bw_dist.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			bw.write("chr, position, p-value, variance_local, variance_global, h1_heritabilities\n");
			bw_dist.write("chr, position");
			for(double alpha=0;alpha<1;alpha+=step){bw_dist.write(", "+alpha);}
			bw_dist.write("\n");
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){
					bw_dist.write((chr+1)+", "+start);
	//				double[][] data=variants.load_variants_in_region(chr, start, end);
					String the_kin_file=local_kinship_folder+"kinship."+variants.sample_size+".chr"
							+(chr+1)+"."+start+"."+end+".rescaled.ibs";
					EMMA emma_local=new EMMA(emma_global.phenotype, genotype, genotype.read_in_kinship(the_kin_file));	
					EMMA emma_new=new EMMA(emma_global.phenotype, genotype,  genotype.read_in_kinship(the_kin_file));		
					if(emma_global.sample_size!=emma_local.sample_size){
						System.out.println("emma_global.sample_size!=emma_local.sample_size");
					}
					double[] best_h={1,0,emma_global.heritability}; // p-value, h_local, h_global
					for(double alpha=0;alpha<1;alpha+=step){
						double[][] new_kinship=new double[emma_global.sample_size][emma_global.sample_size];
						for(int i=0;i<emma_global.sample_size;i++){
							for(int j=0;j<emma_global.sample_size;j++){
								new_kinship[i][j]=emma_local.kinship_matrix[i][j]*alpha+emma_global.kinship_matrix[i][j]*(1-alpha);
							}
						}
						emma_new.clear();
						emma_new.kinship_matrix=VariantsDouble.re_scale_kinship_matrix(new_kinship);
//						check_kinship(emma_new.kinship_matrix);
//						System.out.println("\nwin"+start+":alpha;"+alpha);
						emma_new.REMLE_null();
						bw_dist.write(", "+emma_new.remle_REML+"/"+emma_new.heritability);
						if(!Double.isNaN(emma_new.remle_REML)){
							double p_value=Test.chi2pr(-2.0*(emma_global.remle_REML-emma_new.remle_REML), 1);
							if(p_value<best_h[0]){
								best_h[0]=p_value;
								best_h[1]=emma_new.heritability*alpha;
								best_h[2]=emma_new.heritability*(1-alpha);
							}
						}
					}					
					bw.write((chr+1)+", "+start+", "+best_h[0]+", "+best_h[1]+", "+best_h[2]+", "+(best_h[1]+best_h[2])+"\n");
					bw_dist.write("\n");
					System.out.println("finsihed win"+start);
//					bw.flush();
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}
			}bw.close();
			bw_dist.close();
			make_plot_one_phen(mlr_output_file, mlr_output_file+".dist.csv");
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double[] generate_breaks(int[] lengths_of_chrs){
		double[] breaks=new double[lengths_of_chrs.length];
		breaks[0]=0;
		for(int chr=1;chr<lengths_of_chrs.length;chr++)breaks[chr]=breaks[chr-1]+lengths_of_chrs[chr-1];
		return breaks;
	}
	/*
	 * This function makes five plots immediately after calling 
	 * public void local_VO(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
	 *		String output_file, double step)
	 *	(1) Manhattan plot for p-value of likelihood ratio results.
	 *	(2) Manhattan plot for R2 of likelihood ratio results.
	 *	(3) Heatmap for Bayesian result
	 *	(4) Histogram for regional VC distribution
	 *	(5) Histogram for total variance explained by chromosome, and chr-length as the control		  
	 */
	public static void make_plot_one_phen(String mlr_file, String bayesian_file){
		
		String phe_name=mlr_file.split("/")[mlr_file.split("/").length-1];
		try{
			if((new File(mlr_file)).exists()){
				AssociationResults mlr=new AssociationResults(mlr_file, 2);
				double[] breaks=generate_breaks(lengths_of_chrs);
				double[] locations=new double[mlr.location.length];
				double[] pvalues=new double[mlr.location.length];
				double[] var=new double[mlr.location.length];
				for(int i=0;i<locations.length;i++){
					locations[i]=mlr.location[i]+breaks[mlr.chr[i]-1];		
					pvalues[i]=-Math.log10(mlr.pvalue[i]);
					var[i]=mlr.AdjustedR2[i]*100.0;
				}	
				MyManhattanPlot plot=new MyManhattanPlot(phe_name, "Locations", "-log10(p-value)", mlr.chr,
						locations, pvalues, 2000, 400, mlr_file+".pvalue.png");
				plot=new MyManhattanPlot(phe_name, "Locations", "Variance Explained (%)", mlr.chr,
						locations, var, 2000, 400, mlr_file+".VarR2.png");
			}else{
				System.out.println("NOFILE: "+mlr_file);
			}
			if((new File(bayesian_file)).exists()){
				RegionalDistribution dist=new RegionalDistribution(bayesian_file);
				MyHeatMap heatmap=new MyHeatMap(dist.distribution, phe_name, "Locations (Mb)", 
						"Variance (%)", 120, 100, bayesian_file+".overall.png");
				MyHistogram hist=new MyHistogram(dist.regional_contribution, dist.regional_category,phe_name, "Variance Contributed (%)", "Proportion of Regions (%)",
						700, 400, bayesian_file+".regional.png");
				hist=new MyHistogram(dist.chr_contribution, dist.chr_names, (new String[] {"Variance Component Explained", "Chromosome Length"}), phe_name, null, "Proportion", 
						700, 400, bayesian_file+".chr.png");
			}else{
				System.out.println("NOFILE: "+bayesian_file);
			}
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	
	public static void check_kinship(double[][] kinship_matrix){
		for(int i=0;i<kinship_matrix.length;i++){
			for(int j=i+1;j<kinship_matrix.length;j++){
				if(Double.compare(kinship_matrix[i][j],kinship_matrix[j][i])!=0){
					System.out.println("kinship wrong, i,j");
				}
			}
		}
	}
	
	/*
	 * OLD program trying whether the output are the same to input
	 */
//	public void go_through_windows_no_tiling(String outfile){
//		try{
//			BufferedWriter bw=new BufferedWriter(new FileWriter(outfile));
//			for(int chr=0;chr<variants.num_chrs;chr++){
//				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
//				int start=1, end=tiling_win_size;
//				while(start<last_position){
//					double[][] data=variants.load_variants_in_region(chr, start, end);
//					start=end+1;end=start+tiling_win_size-1;
//					for(int k=0;k<data.length;k++){
//					for(int i=0;i<data[k].length;i++){
//						bw.write(data[k][i]+",");
//						}bw.write("\n");
//					}
//				}
//			}bw.close();
//		}catch(Exception e){e.printStackTrace();}		
//	}
	
	
	
	public static void main(String[] args) {
		long start=System.currentTimeMillis();
		String genotype_csv_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv";
		String phe_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/phen_raw_092910.tsv";
//		String global_kinship_file_ori="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv.kinship";
		String global_kinship_file_rescaled="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv.kinship.rescaled";
		String hdf5_file=genotype_csv_file+".mafc.hdf5";

		String local_kinship_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/check_programs/local_kinship_files/";
		String local_kinship_VO_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/check_programs/local_kinship_VO/";
		
		String to_plot_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/local_results/";
		String sample_file1=to_plot_folder+"DNA_local/200k/Local_VO.29.2052_germ_284_t0.w200000.csv";
		String sample_file2=sample_file1+".dist.csv";
//		make_plot_one_phen(sample_file1, sample_file2);
		
//		int win_size=50000;
//		LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(hdf5_file, win_size, null);
//		double scale=2.0;
//		local_k.calculating_kinship_tiling_windows(local_kinship_folder, scale);
//		
//		MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
//		for(int phe_index=4;phe_index<phenotypeS.num_of_pheno;phe_index++){
//			String out_phe_file=local_kinship_VO_folder+"Local_VO."+phe_index+"_"+phenotypeS.phenotypes[phe_index].phe_id+".csv";
//			local_k.local_VO(phenotypeS.phenotypes[phe_index], hdf5_file, global_kinship_file_rescaled, local_kinship_folder,
//					out_phe_file, 0.01);
//		}
//		System.out.println((System.currentTimeMillis()-start)/1000);


	}
}
