package myMathLib;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;

import org.apache.commons.math3.stat.StatUtils;





public class TryPrograms {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double[] data={0.03077972, -0.66081015,  0.21684005,
				-1.41882095,  0.13929080, -0.79641312,  1.17466404,  0.76569981, -0.56422328, -0.22769212};
//		System.out.println(StatUtils.mean(data));
//		System.out.println(StatUtils.variance(data));
		
		System.out.println(StatFuncs.normality_test_ks(data));
	}

}
