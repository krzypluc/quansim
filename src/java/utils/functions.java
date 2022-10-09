package utils;

import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import scala.collection.Map;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;


public class functions {
    public static Complex waveFunction(Complex x) {
        // valueOfExpotential = -1 * x^2
        Complex valueOfExpontetial = x.multiply(x).multiply(-1);

        // returns e^((-1) * x^2)
        return valueOfExpontetial.exp();
    }

    public static Object[] initYandX(Complex[] y, double[] x, double period, double dx) {
        Complex sumOfValues = Complex.ZERO;
        Complex firstValue = Complex.ZERO;
        Complex lastValue = Complex.ZERO;

        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            x[i] = (-period / 2) + dx * i;
            y[i] = waveFunction(Complex.valueOf(x[i], 0.0));
            sumOfValues = sumOfValues.add(y[i]);

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

        return new Object[]{y, x, integral};
    }
}
