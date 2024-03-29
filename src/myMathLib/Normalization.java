package myMathLib;

import java.util.ArrayList;
import java.util.Collections;

public class Normalization {
	
	public static void normalize_colomns(double[][] data){
		double sums[]=new double[data[0].length];
		for(int k=0;k<data.length;k++){
			for(int i=0;i<data[0].length;i++){
				sums[i]+=data[k][i];
			}
		}
		for(int k=0;k<data.length;k++){
			for(int i=0;i<data[0].length;i++){
				data[k][i]/=sums[i];
			}
		}
	}
	
	public static double median(ArrayList<Double> data){
		Collections.sort(data);
		int index=data.size()/2;
		if(data.size()%2==1){
			return data.get(index);
		}else{
			return (data.get(index)+data.get(index-1))/2.0;
		}
	}
	
	public static double[][] byte2double(byte[][] input){
		if(input==null)return null;
		double[][] out=new double[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=new double[input[i].length];
			for(int j=0;j<input[i].length;j++){
				out[i][j]=input[i][j];
			}
		}
		return out;
	}
	
	public static double[] byte2double(byte[] input){
		if(input==null)return null;
		double[] out=new double[input.length];
		for(int i=0;i<input.length;i++){			
			out[i]=input[i];			
		}
		return out;
	}
	
}
