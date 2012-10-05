package mixedmodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

import nam.Cross_info;
import nam.Founder;
import nam.Imputation;
import nam.RILs;


public class JAWAMix5 {

	static String mit="Copyright (C) 2012, QUAN LONG and QINGRUN ZHANG" +
			"\n\nPermission is hereby granted, free of charge, to any person obtaining a copy " +
			"\nof this software and associated documentation files (the \"Software\"), " +
			"\nto deal in the Software without restriction, including without limitation the " +
			"\nrights to use, copy, modify, merge, publish, distribute, sublicense, " +
			"\nand/or sell copies of the Software, and to permit persons to whom the " +
			"\nSoftware is furnished to do so, subject to the following conditions:" +
			"\n\nThe above copyright notice and this permission notice shall be included " +
			"\nin all copies or substantial portions of the Software." +
			"\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS " +
			"\nOR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, " +
			"\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE " +
			"\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER " +
			"\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING " +
			"\nFROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS " +
			"\nIN THE SOFTWARE.";
	public static void main(String[] args) {
		if(args.length==0){
			System.out.println("==========================================\n"+mit+"\n==========================================\n");
			System.out.println("JAWAMix5: JAva implementation of Whole-genome Association studies " +
					"using Mixed model based on HDF5 (r1.0.0).\n" +
					"Developer: Quan LONG & Qingrun ZHANG\n" +
					"Usage: java -Xmx2g -jar jawamix5.jar [function]");
			System.out.println("Supported functions:" +
					"\n\temmax" +
					"\n\tlm" +
					"\n\temmax_stepwise" +
					"\n\tlm_stepwise" +
					"\n\tlocal" +
					"\n\trare" +
					"\n\tnam"+
					"\n\timport" +
					"\n\tchar2num" +
					"\n\tkinship");
			System.exit(0);
		}
		String function=args[0];
		if(function.equals("kinship")){
			if(args.length==1){
				System.out.println("Compute IBS kinship metrix for the other analysis");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"[-w\ttiling_window_size(bp) (df=WG)]\n\t" +
						"[-scale\tmax_genotype_coding (df=2)]\n\t");
				System.exit(0);
			}else{
				String input=null, output_folder=null;
				int tiling_window_size=-1; 
				double scale=2.0;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-w"))tiling_window_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-scale"))scale=Double.parseDouble(args[k+1]);
					}
				}if(input==null||output_folder==null){
					System.out.println("Input or output-folder can't be null!");
				}else{
					if(tiling_window_size!=-1){
						LocalKinshipAnalyzer analyzer= new LocalKinshipAnalyzer(input, tiling_window_size, null);
						analyzer.calculating_kinship_tiling_windows(output_folder, scale);
					}else{ //	tiling_window_size==-1, i.e., whole-genome				
						VariantsDouble calculator=new VariantsDouble(input);
						System.out.println("Calculating global IBS kinship for "+input);
						calculator.calculate_raw_kinship(output_folder+"kinship.raw.ibs", scale);
						VariantsDouble.re_scale_kinship_matrix(output_folder+"kinship.raw.ibs", output_folder+"kinship.rescaled.ibs");
					}
				}
			}
		}else if(function.equals("char2num")){
			if(args.length==1){
				System.out.println("Convert char-coded genotype CSV file to number-coded genotype CSV file ready for \"import\"");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\toutput_file>\n\t");
				System.exit(0);
			}else{
				String input=null, output=null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
			}
				}if(input==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					VariantsDouble.char2num(input, output);
				}
			}
		}else if(function.equals("import")){
			if(args.length==1){
				System.out.println("Import genotype from .csv file to .hdf5 file");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\tout_put_hdf5_file>\n\t" +
						"[-b\tblock_size (df=5000)]\n\t");
				System.exit(0);
			}else{
				String input=null, output=null;
				int block_size=5000;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
						else if(args[k].equals("-b"))block_size=Integer.parseInt(args[k+1]);
					}
				}if(input==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					VariantsDouble.importCSV(input, output, block_size);
				}
			}
		}else if(function.equals("emmax")){
			if(args.length==1){
				System.out.println("Run EMMAX for phenotype(s)");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");
				System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null, kinship=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100, round=1;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik"))kinship=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| kinship==null||output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1)
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					else{
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}
		}else if(function.equals("emmax_stepwise")){
			if(args.length==1){
				System.out.println("Run EMMAX-based stepwise regression for phenotype(s)");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"[-r\t<round> (df=2)]\n\t"+
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");
				System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null, kinship=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100, round=2;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik"))kinship=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| kinship==null||output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1)
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					else{
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}
		}else if(function.equals("lm_stepwise")){
			if(args.length==1){
				System.out.println("Run Stepwise regression for phenotype(s) without population structure accounted.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"[-r\t<round> (df=2)]\n\t"+
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");
				System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100, round=2;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1){
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					}else{
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}
		}else if(function.equals("lm")){
			if(args.length==1){
				System.out.println("Run linear regression for phenotype(s) without population structure accounted.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100;
				final int round=1;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1){
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					}else{
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}
		}else if(function.equals("local")){
			if(args.length==1){
				System.out.println("Run local variance component analysis");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-w\ttiling_window_size(bp)>\n\t" +
						"<-ik_l\tlocal_kinship_files_folder>\n\t" +
						"<-ik_g\tglobal_kinship_file>\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-step\t<grid_size> (df=0.01)]");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_folder=null, 
						local_kinship_files_folder=null, global_kinship_file=null;
				int win_size=-1;
				int the_phe_index=-1;
				double step=0.01;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik_l"))local_kinship_files_folder=args[k+1];
						else if(args[k].equals("-ik_g"))global_kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-w"))win_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-step"))step=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| local_kinship_files_folder==null||output_folder==null 
						||global_kinship_file==null || win_size==-1){
					System.out.println("Window size or input or output can't be null!");
				}else{
					LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(input_geno, win_size, null);
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1){
						System.out.println("Running all phenotypes? It is suggested to specify a phenotype index. Otherwise it may be slow." +
								"\nLet us have a try!");
						for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
							String out_phe_file=output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".w"+win_size+".csv";
							local_k.local_VO(phenotypeS.phenotypes[phe_index], input_geno, global_kinship_file, local_kinship_files_folder,
									out_phe_file, step);
						}
					}else{
						String out_phe_file=output_folder+"Local_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".w"+win_size+".csv";
						local_k.local_VO(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, local_kinship_files_folder,
								out_phe_file, step);
					}
				}
			}
		}else if(function.equals("rare")){
			if(args.length==1){
				System.out.println("Run rare variants analysis leveraging (or not) potential sythetic associations, " +
						"with population structure acounted (or not)");
				System.out.println("Usage:\n\t" +
						"<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"<-is\tinput_single_marker_results (preferably by emmax)>\n\t" +
						"<-ir\tregions_grouping_file>\n\t" +
						"<-o\toutput_prefix>\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-ld\t<LD_threshold(D-prime)> (df=0.8)]\n\t" +
						"[-rare\t<rare_threshold> (df=0.01)]\n\t" +
						"[-syn_pvalue\t<potential syntehtic association pvalue threshold> (df=0.00001)]\n\t" +
						"[-syn_maf\t<potential syntehtic association MAF threshold> (df=0.00001)]\n\t" +
						"[-dist\t<distance2region>(df=100000)]");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_prefix=null, kinship_file=null, input_single_marker_results_file=null, 
						input_region_file=null;
				double ld_threshold=0.8, rare_threshold=0.01, synthetic_pvalue_threshold=0.00001, synthetic_maf_threhold=0.1;
				int the_phe_index=-1, distance2region=100000;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-is"))input_single_marker_results_file=args[k+1];
						else if(args[k].equals("-ir"))input_region_file=args[k+1];
						else if(args[k].equals("-ik"))kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_prefix=args[k+1];
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-dist"))distance2region=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-ld"))ld_threshold=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-rare"))rare_threshold=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-syn_pvalue"))synthetic_pvalue_threshold=Double.parseDouble(args[k+1]);						
						else if(args[k].equals("-syn_maf"))synthetic_maf_threhold=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| kinship_file==null||output_prefix==null || input_single_marker_results_file==null
						|| input_region_file==null){
					System.out.println("Input or output can't be null!");
				}else{
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1)
						for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
							RareAnalyzer rare_analyzer=new RareAnalyzer(input_single_marker_results_file, 
									synthetic_pvalue_threshold, synthetic_maf_threhold, 
									input_geno, input_region_file, distance2region, phenotypeS.phenotypes[phe_index], kinship_file,
									ld_threshold, rare_threshold);
							rare_analyzer.rare_association_4in1(rare_threshold, output_prefix);
						}
					else{
						RareAnalyzer rare_analyzer=new RareAnalyzer(input_single_marker_results_file, 
								synthetic_pvalue_threshold, synthetic_maf_threhold, 
								input_geno, input_region_file, distance2region, phenotypeS.phenotypes[the_phe_index], kinship_file,
								ld_threshold, rare_threshold);
						rare_analyzer.rare_association_4in1(rare_threshold, output_prefix);
					}
				}
			}
		}else if(function.equals("nam")){
			System.out.println("Nested Association Mapping"+"\n");	
			if(args.length<10){
				System.out.println("Usage: \n\t"
						+"<-p\tlist of RIL_pedigreee files folder>\n\t"
						+ "<-c\tCross_ID>\n\t" 
						+ "<-fg\tfounder whole genome genotype file>\n\t"
						+"<-ril\tRIL phenotype>\n\t"
						+ "<-o\tOutput-prefix>\n\t"
						+"[-b\tblock size in hdf5 (df=5000)]\n\t"
						+"[-index\tphenotype_index (df=ALL, start from zero)]\n\t"
						+"[-r\tround (df=2)]\n\t"
						+"[-p\tpvalue_after_multi.correct.> (df=1000)]\n\t"
						+"[-maf\tmaf_threshold_plot (df=0.05)]\n\t");			
				System.exit(0);
			}						
			String pedigreefiles_folder =null;	
			String crossid =null;
			String RILpheno =null;
			String geno250k =null;	
			String output_prefix=null;
			int block_size=5000, the_phe_index=-1;
			
			double p_after_corr=1000;
			double maf_threshold_plot=0.05; 
			int round=2;
			for(int k=0; k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-p")) pedigreefiles_folder=args[k+1];
					else if(args[k].equals("-c"))  crossid=args[k+1];
					else if(args[k].equals("-ril")) RILpheno=args[k+1];
					else if(args[k].equals("-fg")) geno250k=args[k+1];
					else if(args[k].equals("-o"))  output_prefix=args[k+1];	
					else if(args[k].equals("-b"))  block_size=Integer.parseInt(args[k+1]);	
					else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
					else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
					else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
					else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
				}			
			}if(pedigreefiles_folder==null||crossid==null||RILpheno==null||geno250k==null||output_prefix==null){
				System.out.println("Input or output paths can't be empty!");				
				System.exit(0);
			}
			VariantsDouble.run_nam_imputation(pedigreefiles_folder, crossid, RILpheno, geno250k, output_prefix, block_size);
			LM.stepwise_reg_nam(output_prefix+".geno.num.csv.hdf5", output_prefix+".pheno.tsv", the_phe_index,
					output_prefix, p_after_corr, maf_threshold_plot, round);
		}else{
			System.out.println("Typo? \""+function+"\" is not a supported function");
			System.exit(0);
		}

	}

}
