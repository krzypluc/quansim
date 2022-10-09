package utils;

import compute.Complex;

public class functions {
    public static Complex waveFunction(Complex x){
        Complex arg = x.multiply(x);
        arg = arg.multiply(-1);

        double epowa = Math.exp(arg.getReal());
        double cosb = Math.cos(arg.getImag());
        double sinb = Math.sin(arg.getImag());

        return new Complex(epowa * cosb, epowa * sinb);
    }
}
