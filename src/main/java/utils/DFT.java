package utils;

import compute.Complex;

import java.util.ArrayList;
import java.util.List;

public class DFT {
    public static Complex[] dft (Complex[] values){

        if (values.length == 1){
            return values;
        }

        Complex[] dftTransform = new Complex[values.length];

        for (int i = 0; i < dftTransform.length; i++){
            dftTransform[i] = new Complex();

            for (int j = 0; j < values.length; j++){
                Complex valueOfExpotential = new Complex(0, -2 * Math.PI * j * i / values.length);
                Complex w = Complex.expotential(valueOfExpotential);
                dftTransform[i] = dftTransform[i].add(values[j].multiply(w));
            }
        }

        return dftTransform;
    }
}
