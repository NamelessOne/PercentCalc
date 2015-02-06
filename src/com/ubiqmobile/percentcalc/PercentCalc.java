package com.ubiqmobile.percentcalc;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import umontreal.iro.lecuyer.functionfit.SmoothingCubicSpline;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Created by NamelessOne
 * 03.02.2015.
 */
public class PercentCalc {
    private static final double THRESHOLD_STRETCH = 0.8;
    private static final double THRESHOLD_TIME = 0.01;

    private static final double AMPLITUDE_MIN = 1.1;
    private static final double AMPLITUDE_MAX = 1.6;

    private static final double DURATION_MIN = 0.03;
    private static final double DURATION_MAX = 0.11;

    private static final int L = 6;

    private static final double DT = 0.001;

    //public static final String FILENAME = "AccTest_Dec.csv";    //есть
    //public static final String FILENAME = "AccTestbefore_Dec.csv";  //есть
    //public static final String FILENAME = "AccTest_Nov.csv";    //есть
    public static final String FILENAME = "AccTestbefore_Nov.csv";  //есть

    private List<Double> t = new ArrayList<Double>();
    private List<Double> acc_tt = new ArrayList<Double>();
    private List<Integer> ind_max = new ArrayList<Integer>();
    private List<Integer> ind_min = new ArrayList<Integer>();
    private List<Double> stretch = new ArrayList<Double>();
    //private List<Double> mean_stretch = new ArrayList<Double>();
    private List<Double> var_stretch = new ArrayList<Double>();
    private double percent;

    public double percent_calc() throws IOException {
        List<Double> acc_norm = read_acc();
        /*if(FILENAME.equals("AccTest_Nov.csv"))
        {
            t.set(77, (t.get(76) + t.get(78)) / 2);
        }*/
        //Один конкретный баг в файле данных. В джаве не нужно.
        List<Double> tt = new ArrayList<Double>();
        double tValue = 0.0;
        while (tValue < t.get(t.size() - 1)) {
            tt.add(tValue);
            tValue += DT;
        }
        double[] a = new double[t.size()];
        double[] b = new double[acc_norm.size()];
        int i = 0;
        for (Double value : t) {
            a[i] = value;
            i++;
        }
        i = 0;
        for (Double value : acc_norm) {
            b[i] = value;
            i++;
        }
        PolynomialSplineFunction f = new SplineInterpolator().interpolate(a, b);
        for (double value : tt) {
            acc_tt.add(f.value(value));
            //System.out.println(f.value(value));
        }

        for (int j = 1; j < acc_tt.size() - 1; j++) {
            if (acc_tt.get(j) >= acc_tt.get(j + 1) && acc_tt.get(j) > acc_tt.get(j - 1))
                ind_max.add(j);
            if (acc_tt.get(j) <= acc_tt.get(j + 1) && acc_tt.get(j) < acc_tt.get(j - 1))
                ind_min.add(j);
        }

        if (ind_min.size() > ind_max.size()) {
            ind_min.remove(0);
        }
        if (ind_max.size() > ind_min.size()) {
            ind_max.remove(0);
        }

        for (int j = 0; j < ind_max.size(); j++) {
            stretch.add(acc_tt.get(ind_max.get(j)) - acc_tt.get(ind_min.get(j)));
        }
        double[] mean_stretch = new double[ind_max.size() - 2 * L];
        for (int i_max = L; i_max < ind_max.size() - L; i_max++) {
            DescriptiveStatistics ds = new DescriptiveStatistics();
            for (double value : stretch.subList(i_max - L, i_max + L + 1))
                ds.addValue(value);
            mean_stretch[i_max - L] = ds.getMean();
            var_stretch.add(ds.getStandardDeviation());
        }
        double[] t_stretch = new double[ind_max.size() - 2 * L];
        for (int j = L; j < ind_max.size() - L; j++) {
            t_stretch[j - L] = tt.get(ind_max.get(j));
        }
        SmoothingCubicSpline fit = new SmoothingCubicSpline(t_stretch, mean_stretch, 0.5);
        double[] mean_stretch_smooth = new double[t_stretch.length];
        for (int j = 0; j < t_stretch.length; j++) {
            mean_stretch_smooth[j] = fit.evaluate(t_stretch[j]);
        }
        boolean[] chi_ampl_stretch = new boolean[mean_stretch_smooth.length + 2 * L];
        for (int j = 0; j < chi_ampl_stretch.length; j++) {
            if (j < L || j > chi_ampl_stretch.length - L - 1)
                chi_ampl_stretch[j] = false;
            else if (mean_stretch_smooth[j - L] < 3)
                chi_ampl_stretch[j] = false;
            else
                chi_ampl_stretch[j] = true;
        }
        double[] der_mean_stretch = new double[mean_stretch_smooth.length - 1];
        for (int j = 0; j < der_mean_stretch.length; j++) {
            der_mean_stretch[j] = mean_stretch_smooth[j + 1] - mean_stretch_smooth[j];
        }
        boolean[] chi_der_stretch = new boolean[der_mean_stretch.length + 1 + 2 * L];
        for (int j = 0; j < chi_der_stretch.length; j++) {
            if (j <= L || j >= chi_der_stretch.length - L - 1)
                chi_der_stretch[j] = false;
            else {
                chi_der_stretch[j] = (Math.abs(der_mean_stretch[j - L]) < THRESHOLD_STRETCH / L) & (Math.abs(der_mean_stretch[j - L - 1]) < THRESHOLD_STRETCH / L);
            }
        }
        boolean[] chi_der_stretch_temp = chi_der_stretch.clone();
        for (int j = 2; j < chi_der_stretch.length - 2; j++) {
            chi_der_stretch_temp[j] = chi_der_stretch[j - 2] & chi_der_stretch[j - 2] & chi_der_stretch[j] & chi_der_stretch[j + 1] & chi_der_stretch[j + 2];
        }
        chi_der_stretch = chi_der_stretch_temp;
        double[] len_max = new double[ind_max.size() - 1];
        for (int j = 0; j < ind_max.size() - 1; j++) {
            len_max[j] = tt.get(ind_max.get(j + 1)) - tt.get(ind_max.get(j));
        }
        DelOutliersStruct dos_struct = del_outliers(len_max, 3);
        int[] ind_cluster = dos_struct.ind_out;
        double sigma = dos_struct.sigma;
        double[] len_max_sort = dos_struct.y;
        double[] ind_max_sort = new double[ind_cluster.length];
        for (int j = 0; j < ind_max_sort.length; j++) {
            ind_max_sort[j] = ind_max.get(ind_cluster[j]);
        }
        double[] mean_time_sort = new double[len_max_sort.length - 2 * L];
        double[] var_time_sort = new double[len_max_sort.length - 2 * L];
        for (int i_max = L; i_max < len_max_sort.length - L; i_max++) {
            DescriptiveStatistics ds = new DescriptiveStatistics();
            for (int j = i_max - L; j < i_max + L + 1; j++)
                ds.addValue(len_max_sort[j]);
            mean_time_sort[i_max - L] = ds.getMean();
            //mean_time_sort[i_max - L] = mean(len_max_sort, i_max - L, i_max + L);
            // var_time_sort[i_max - L] = getStdDev(mean_time_sort[i_max - L], len_max_sort, i_max - L, i_max + L);
            var_time_sort[i_max - L] = ds.getStandardDeviation();
        }
        double[] var_time_max = new double[var_time_sort.length + 2 * L];
        for (int j = 0; j < var_time_max.length; j++) {
            if (j < L || j > var_time_max.length - L - 1) {
                var_time_max[j] = 0;
            } else {
                var_time_max[j] = var_time_sort[j - L];
            }
        }
        double[] t_time = new double[ind_cluster.length - 2 * L];
        for (int j = 0; j < t_time.length; j++) {
            t_time[j] = tt.get(ind_max.get(j + L));
        }
        SmoothingCubicSpline scs = new SmoothingCubicSpline(t_time, mean_time_sort, 0.5);
        double[] mean_time_smooth = new double[t_time.length];
        for (int j = 0; j < t_time.length; j++) {
            mean_time_smooth[j] = scs.evaluate(t_time[j]);
        }
        double[] der_mean_time = new double[mean_time_smooth.length - 1];
        for (int j = 0; j < der_mean_time.length; j++) {
            der_mean_time[j] = mean_time_smooth[j + 1] - mean_time_smooth[j];
        }
        boolean[] chi_der_time = new boolean[der_mean_time.length + 1 + 2 * L];
        for (int j = 0; j < chi_der_time.length; j++) {
            if (j <= L || j >= chi_der_time.length - L - 1)
                chi_der_time[j] = false;
            else {
                chi_der_time[j] = (Math.abs(der_mean_time[j - L]) < THRESHOLD_TIME / L) & (Math.abs(der_mean_time[j - L - 1]) < THRESHOLD_TIME / L);
            }
        }
        boolean[] chi_time = new boolean[ind_max.size() - 1];
        int ind = 0;
        for (int value : ind_cluster) {
            chi_time[value] = chi_der_time[ind];
            ind++;
        }
        double[] var_time = new double[ind_max.size() - 1];
        ind = 0;
        for (int value : ind_cluster) {
            var_time[value] = var_time_max[ind];
            ind++;
        }
        if (sum(ind_max) > sum(ind_min)) {
            boolean[] temp_chi_time = new boolean[chi_time.length + 1];
            temp_chi_time[0] = chi_time[0];
            for (int j = 1; j < chi_time.length; j++) {
                temp_chi_time[j] = chi_time[j - 1];
            }
            chi_time = temp_chi_time;
            double[] temp_var_time = new double[var_time.length + 1];
            temp_var_time[0] = var_time[0];
            for (int j = 1; j < var_time.length; j++) {
                temp_var_time[j] = var_time[j - 1];
            }
            var_time = temp_var_time;
        } else {
            boolean[] temp_chi_time = new boolean[chi_time.length + 1];
            temp_chi_time[chi_time.length] = chi_time[chi_time.length - 1];
            for (int j = 1; j < chi_time.length; j++) {
                temp_chi_time[j] = chi_time[j];
            }
            chi_time = temp_chi_time;
            double[] temp_var_time = new double[var_time.length + 1];
            temp_var_time[var_time.length] = var_time[var_time.length - 1];
            for (int j = 1; j < var_time.length; j++) {
                temp_var_time[j] = var_time[j];
            }
            var_time = temp_var_time;
        }
        List<Integer> ind_full = new ArrayList<Integer>();
        for(int j = 0; j < chi_time.length; j++)
        {
            if(chi_time[j]&chi_ampl_stretch[j]&chi_der_stretch[j])
            {
                ind_full.add(j);
            }
        }
        if(ind_full.isEmpty())
            throw new IOException("Not enough information");
        else
        {
            DescriptiveStatistics ds = new DescriptiveStatistics();
            for(int value: ind_full)
            {
                ds.addValue(var_stretch.get(value-L));
            }
            double amplitude = ds.getMean();
            ds.clear();
            for(int value: ind_full)
            {
                ds.addValue(var_time[value]);
            }
            double duration = ds.getMean();
            if((amplitude>AMPLITUDE_MAX)||(duration>DURATION_MAX))
            {
                percent = 100;
            }
            else
            {
                if((amplitude < AMPLITUDE_MIN)||(duration<DURATION_MIN))
                {
                    percent = 0;
                }
                else 
                {
                    percent = Math.max((amplitude-AMPLITUDE_MIN)/(AMPLITUDE_MAX-AMPLITUDE_MIN), (duration-DURATION_MIN)/(DURATION_MAX-DURATION_MIN));
                }
            }
        }
        return percent;
    }

    private List<Double> read_acc() throws IOException {
        BufferedReader fid = new BufferedReader(new FileReader(FILENAME));
        fid.readLine(); //skip first line
        int i_line = 0;
        String t_line;
        int ind_bracket;
        List<Double[]> acc = new ArrayList<Double[]>();

        while ((t_line = fid.readLine()) != null) {
            ind_bracket = t_line.indexOf('(');
            List<Integer> ind_comma = new ArrayList<Integer>();
            int i = 0;
            while (t_line.indexOf(',', i) != -1) {
                ind_comma.add(t_line.indexOf(',', i));
                i = t_line.indexOf(',', i) + 1;
            }
            if (ind_comma.size() == 3) {
                t.add(i_line, Double.valueOf(t_line.substring(5, ind_bracket)) / 1000);
                Double[] d = new Double[3];
                d[0] = Double.valueOf(t_line.substring(ind_comma.get(0) + 2, ind_comma.get(1)));
                d[1] = Double.valueOf(t_line.substring(ind_comma.get(1) + 2, ind_comma.get(2)));
                d[2] = Double.valueOf(t_line.substring(ind_comma.get(2) + 2));
                acc.add(i_line, d);
                i_line = i_line + 1;
            }
        }
        double t_start = t.get(0);
        for (int i = 0; i < t.size(); i++) {
            t.set(i, t.get(i) - t_start);
        }
        List<Double> res = new ArrayList<Double>();
        for (Double[] d : acc) {
            //System.out.println(Math.pow(d[0] + d[1] + d[2], 2.2));
            res.add(Math.sqrt(Math.pow(d[0], 2) + Math.pow(d[1], 2) + Math.pow(d[2], 2)));
        }
        return res;
    }

    private DelOutliersStruct del_outliers(double[] x, int k) {
        List<Double> y = new ArrayList<Double>();
        for (double d : x) {
            y.add(d);
        }
        List<Double> y_out = new ArrayList<Double>();
        for (int i = 0; i < y.size() + 1; i++) {
            y_out.add(0.0);
        }
        //double[] y = x.clone();
        //double[] y_out = new double[y.length+1];
        List<Integer> ind_out = new ArrayList<Integer>();
        for (int i = 0; i < x.length; i++) {
            ind_out.add(i);
        }
        double sigma = 0;
        while (y_out.size() > y.size()) {
            y_out = new ArrayList<Double>(y);
            DescriptiveStatistics ds = new DescriptiveStatistics();
            for (double value : y)
                ds.addValue(value);
            double m = ds.getMean();
            sigma = ds.getStandardDeviation();
            List<Integer> ind = new ArrayList<Integer>();
            for (int i = 0; i < y.size(); i++) {
                if ((y.get(i) - m) * (y.get(i) - m) < k * k * sigma * sigma)
                    ind.add(i);
            }
            for (int i = 0; i < y.size(); i++) {
                if (!ind.contains(i))
                    y.remove(i);
            }
            for (int i = 0; i < ind_out.size(); i++) {
                if (!ind.contains(i))
                    ind_out.remove(i);
            }
        }
        double[] res_y = new double[y.size()];
        for (int i = 0; i < y.size(); i++) {
            res_y[i] = y.get(i);
        }
        int[] res_ind_out = new int[ind_out.size()];
        for (int i = 0; i < ind_out.size(); i++) {
            res_ind_out[i] = ind_out.get(i);
        }
        return new DelOutliersStruct(res_y, res_ind_out, sigma);
    }

    private class DelOutliersStruct {
        public double sigma;
        public int[] ind_out;
        public double[] y;

        public DelOutliersStruct(double[] y, int[] ind_out, double sigma) {
            this.sigma = sigma;
            this.ind_out = ind_out;
            this.y = y;
        }
    }

    private double sum(List<Integer> array) {
        double sum = 0;
        for (double value : array) {
            sum += value;
        }
        return sum;
    }
}
