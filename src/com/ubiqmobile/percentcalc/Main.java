package com.ubiqmobile.percentcalc;

import java.io.IOException;

public class Main {
    public static void main(String[] args) {
	    try
        {
            new PercentCalc().percentCalc();
        }catch (IOException e){
            e.printStackTrace();
        }
    }
}
