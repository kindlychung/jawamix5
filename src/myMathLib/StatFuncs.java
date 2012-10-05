package myMathLib;

import java.util.Arrays;

import org.apache.commons.math3.distribution.KolmogorovSmirnovDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

public class StatFuncs {
	
	/*
	 *  Kolmogorov–Smirnov test
	 */
	public static void ks_test(double[] x, double[] y){
		double[] xx=x.clone();
		double[] yy=y.clone();
		Arrays.sort(xx);
		Arrays.sort(yy);
	}
	
	/*
	 * normality_test using K-S. 
	 */
	public static double normality_test_ks(double[] x){
		double[] data=x.clone(); // don't modify the original data x!
		Arrays.sort(data);
		double mean=StatUtils.mean(data);
		double var=StatUtils.variance(data,mean);
		NormalDistribution normal=new NormalDistribution(mean, Math.sqrt(var));
		double n=data.length;
		double Dn=0;
		for(int i=0;i<n;i++){
			double Fx=normal.cumulativeProbability(data[i]);
			double Fnx=(i+1)/n;			
			if(Math.abs(Fnx-Fx)>Dn)Dn=Math.abs(Fnx-Fx);
		}
		KolmogorovSmirnovDistribution ks_dist=new KolmogorovSmirnovDistribution((int)n);
		double pvalue=1-ks_dist.cdf(Dn);	
		return pvalue;
	}
	
	public static boolean variable_or_not(double[] data){
		for(int k=1;k<data.length;k++){
			if(Double.compare(data[0],data[k])!=0)return true;
		}return false;
	}
	
	public static double[] multi_reg_pvalues(double[] b, double[] sde, int n){
		TDistribution t_dis=new TDistribution(n-b.length);
		double[] t=new double[b.length];
		for(int k=0;k<b.length;k++)t[k]=Math.abs(b[k]/sde[k]);
		double[] pvalues=new double[b.length];
		for(int k=0;k<b.length;k++)pvalues[k]=2*(1-t_dis.cumulativeProbability(t[k]));
		return pvalues;
	}
	
	public static double correlationPearsons(double[] v1, double[] v2){
		PearsonsCorrelation calculator=new PearsonsCorrelation();
		return calculator.correlation(v1, v2);
	}
	
	public static double correlationSpearmans(double[] v1, double[] v2){
		SpearmansCorrelation calculator=new SpearmansCorrelation();
		return calculator.correlation(v1, v2);
	}
}
