package utils;


import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;

public class Integrator {
    public static Complex TrapezoidComplex1D(Complex[] y, double[] x, double dx) {
        Complex sumOfValues = Complex.ZERO;
        Complex firstValue = Complex.ZERO;
        Complex lastValue = Complex.ZERO;

        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
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

        return integral;
    }

    public static double TrapezoidDouble1D(double[] y, double[] x) {
        double sum = 0;
        double dx = x[1] - x[0];
        double val;

        for (int i = 0; i < (y.length - 1); i++) {
            val = (y[i] + y[i + 1]) * dx / 2;
            sum += val;
        }

        return sum;
    }
}
