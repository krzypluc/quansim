package utils;


import org.apache.commons.math3.complex.Complex;

public class integrator {
    public static Complex TrapezoidComplex1D(Complex[] y, double[] x) {
        double dx = x[1] - x[0];

        double sumReal = 0, sumImag = 0;
        double valReal, valImag;

        for (int i = 0; i < (y.length - 1); i++) {
            valReal = (y[i].getReal() + y[i + 1].getReal()) * dx / 2;
            valImag = (y[i].getImaginary() + y[i + 1].getImaginary()) * dx / 2;;
            sumReal += valReal;
            sumImag += valImag;
        }

        return Complex.valueOf(sumReal, sumImag);
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
