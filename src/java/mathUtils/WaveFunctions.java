package mathUtils;

import org.apache.commons.math3.complex.Complex;
import org.yaml.snakeyaml.Yaml;

import java.io.*;
import java.util.*;


public class WaveFunctions {

    public enum WaveFunction {
        GAUSSIAN,
        DIRACSDELTA
    }

    public static Complex waveFunction(double x, WaveFunction func) {
        Complex value;

        switch (func) {
            case GAUSSIAN:
                // valueOfExpotential = -1 * x^2
                value = Complex.valueOf(Math.pow(x, 2) * (-1));

                // returns e^((-1) * x^2)
                value = value.exp();
                return value;

            case DIRACSDELTA:
                if (-2 > x){
                    return Complex.ZERO;
                }
                else{
                    value = Complex.I.multiply(x).exp();
                    return value;
                }


            default:
                throw new RuntimeException("There is no such implemented WaveFunction.");
        }
    }
}

