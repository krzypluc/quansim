package utils;

import org.apache.commons.numbers.complex.Complex;

public class DFT {
    public static Complex[] dft (Complex[] values){

        if (values.length == 1){
            return values;
        }

        Complex[] dftTransform = new Complex[values.length];

        for (int i = 0; i < dftTransform.length; i++){
            dftTransform[i] = Complex.ofCartesian(0, 0);

            for (int j = 0; j < values.length; j++){
                Complex valueOfExpotential = Complex.ofCartesian(0, -2 * Math.PI * j * i / values.length);
                Complex w = valueOfExpotential.exp();
                dftTransform[i] = dftTransform[i].add(values[j].multiply(w));
            }
        }

        return dftTransform;
    }
}
