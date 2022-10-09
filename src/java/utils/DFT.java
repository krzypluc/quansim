package utils;


import org.apache.commons.math3.complex.Complex;

public class DFT {
    public static Complex[] forwardDFT (Complex[] values){

        if (values.length == 1){
            return values;
        }

        Complex[] dftTransform = new Complex[values.length];

        for (int i = 0; i < dftTransform.length; i++){
            dftTransform[i] = Complex.valueOf(0, 0);
            
            for (int j = 0; j < values.length; j++){
                Complex valueOfExpotential = Complex.valueOf(0, -2 * Math.PI * j * i / values.length);
                Complex w = valueOfExpotential.exp();
                dftTransform[i] = dftTransform[i].add(values[j].multiply(w));
            }
        }

        return dftTransform;
    }

    public static Complex[] inverseDFT(Complex[] values){
        Complex[] dftTransform = new Complex[values.length];
        return dftTransform;
    }
}
