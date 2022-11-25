package mathUtils;

import org.apache.commons.math3.analysis.function.HarmonicOscillator;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.yaml.snakeyaml.Yaml;

import java.io.*;
import java.util.*;


public class Functions {
    public enum Potential {
        HARMONIC,
        NONE
    }

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

    public static double potential(double x){
        double value;
        // harmonic
        // value = Math.pow(x, 2) / 2;

        // no potential
        value = 0;

        return value;
    }

    public static HashMap<String, Object> initYandX(Complex[] y, double[] x, double period, double dx) {
        HashMap<String, Object> ret = new HashMap<String, Object>();

        Complex sumOfValues = Complex.ZERO;
        Complex firstValue = Complex.ZERO;
        Complex lastValue = Complex.ZERO;

        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;
        double[] potential = new double[x.length];

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            x[i] = (-period / 2) + dx * i;
            y[i] = waveFunction(x[i], WaveFunction.GAUSSIAN);
            sumOfValues = sumOfValues.add(y[i].multiply(y[i]));
            potential[i] = potential(x[i]);

            // Saving first value - use in integral
            if (i == procID * lengthOfpiece) {
                firstValue = y[i];
            }

            // Saving last value - use in integral
            if (i == (procID + 1) * lengthOfpiece - 1) {
                lastValue = y[i];
            }
        }

        Complex integral = sumOfValues
                .multiply(2.0)
                .subtract(lastValue)
                .subtract(firstValue)
                .multiply(dx / 2);

        ret.put("x", x);
        ret.put("y", y);
        ret.put("integral", integral);
        ret.put("potential", potential);
        return ret;
    }
    
    public static Map<String, Object> loadConfigFromYaml(String path) throws IOException {
        InputStream inputStream = new FileInputStream(new File(path));

        Yaml yaml = new Yaml();
        return yaml.load(inputStream);
    }
}

