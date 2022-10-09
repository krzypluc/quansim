package utils;

import org.apache.commons.numbers.complex.Complex;
import org.apache.commons.numbers.complex.Complex;

public class functions {
    public static Complex waveFunction(Complex x) {
        // valueOfExpotential = -1 * x^2
        Complex valueOfExpontetial = x.multiply(x).multiply(-1);

        // returns e^((-1) * x^2)
        return valueOfExpontetial.exp();
    }
}
