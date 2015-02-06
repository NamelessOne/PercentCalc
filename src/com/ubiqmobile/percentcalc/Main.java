package com.ubiqmobile.percentcalc;

import java.io.IOException;

public class Main {
    public static void main(String[] args) {
	    try
        {
            System.out.println(new PercentCalc().percent_calc());
        }catch (IOException e){
            e.printStackTrace();
        }
    }
}
