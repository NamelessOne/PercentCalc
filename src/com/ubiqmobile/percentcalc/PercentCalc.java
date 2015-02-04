package com.ubiqmobile.percentcalc;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Created by NamelessOne
 * 03.02.2015.
 */
public class PercentCalc {
    private static final double THRESHOLD_STRETCH = 0.8;
    private static final double THRESHOLD_TIME = 0.01;

    private static final double AMPLITUDE_MAN = 1.1;
    private static final double AMPLITUDE_MAX = 1.6;

    private static final double DURATION_MIN = 0.03;
    private static final double DURATION_MAX = 0.11;

    private static final int L = 6;

    private static final double DT = 0.001;

    //public static final String FILENAME = "AccTest_Dec.csv";
    //public static final String FILENAME = "AccTestbefore_Dec.csv";
    //public static final String FILENAME = "AccTest_Nov.csv";
    public static final String FILENAME = "AccTestbefore_Nov.csv";

    private List<Double> t = new ArrayList<Double>();
    private List<Double> acc_norm;
    private List<Double> acc_tt = new ArrayList<Double>();
    private List<Integer> ind_max = new ArrayList<Integer>();
    private List<Integer> ind_min = new ArrayList<Integer>();
    private List<Double> stretch = new ArrayList<Double>();


    public void percentCalc() throws IOException
    {
        acc_norm = readAcc();
        List<Double> tt = new ArrayList<Double>();
        double tValue = 0.0;
        while(tValue<t.get(t.size()-1))
        {
            tt.add(tValue);
            tValue += DT;
        }
        double[] a = new double[t.size()];
        double[] b = new double[acc_norm.size()];
        int i = 0;
        for(Double value: t)
        {
            a[i] = value;
            i++;
        }
        i = 0;
        for(Double value: acc_norm)
        {
             b[i] = value;
             i++;
        }
        PolynomialSplineFunction f = new SplineInterpolator().interpolate(a, b);
        for(double value: tt)
        {
            acc_tt.add(f.value(value));
            //System.out.println(f.value(value));
        }

        for (int j=1; j<acc_tt.size()-1; j++)
        {
            if(acc_tt.get(j)>=acc_tt.get(j+1)&&acc_tt.get(j)>acc_tt.get(j-1))
                ind_max.add(j);
            if(acc_tt.get(j)<=acc_tt.get(j+1)&&acc_tt.get(j)<acc_tt.get(j-1))
                ind_min.add(j);
        }

        if(ind_min.size()> ind_max.size())
        {
            ind_min.remove(0);
        }
        if(ind_max.size()> ind_min.size())
        {
            ind_max.remove(0);
        }

        for (int j = 0; j<ind_max.size(); j++)
        {
            System.out.print("min = " + acc_tt.get(ind_min.get(j)) + "_" + ind_min.get(j));
            System.out.println(" max = " + acc_tt.get(ind_max.get(j)) + "_" + ind_max.get(j));
            stretch.add(acc_tt.get(ind_max.get(j)) - acc_tt.get(ind_min.get(j)));
            System.out.println("stretch = " + String.valueOf(acc_tt.get(ind_max.get(j)) - acc_tt.get(ind_min.get(j))));
        }


        //CubicSpline cs = new CubicSpline(a, b);
        //Один конкретный баг в файле данных. В джаве не нужно.
        /*if(FILENAME.equals("AccTest_Nov.csv"))
        {
            //t(77) = (t(76) + t(78))/2;
        }*/
    }

    private List<Double> readAcc() throws IOException
    {
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
                t.add(i_line, Double.valueOf(t_line.substring(6, ind_bracket)) / 1000);
                Double[] d = new Double[3];
                d[0] = Double.valueOf(t_line.substring(ind_comma.get(0) + 2, ind_comma.get(1) - 1));
                d[1] = Double.valueOf(t_line.substring(ind_comma.get(1) + 2, ind_comma.get(2) - 1));
                d[2] = Double.valueOf(t_line.substring(ind_comma.get(2) + 2));
                acc.add(i_line, d);
                i_line = i_line + 1;
            }
        }
        double t_start = t.get(0);
        for(int i = 0; i < t.size(); i++)
        {
            t.set(i,  t.get(i) - t_start);
        }
        List<Double> res = new ArrayList<Double>();
        for (Double[] d: acc)
        {
            //System.out.println(Math.pow(d[0] + d[1] + d[2], 2.2));
            res.add(Math.sqrt(Math.pow(d[0], 2) + Math.pow(d[1], 2) + Math.pow(d[2], 2)));
        }
        return res;
    }

    private void delOutliers()
    {

    }
}
