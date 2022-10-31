package mathUtils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;

public class BesselFunctions {
    private static final int maxEval = (int) Math.pow(2, 6);

    public static double Besselnx(int n, double x) {
        double startOfIntegral = 0.0;
        double endOfIntegral = Math.PI;
        double multiplier = 1 / Math.PI;

        IntegralBesselFunction bessFunction = new IntegralBesselFunction(n, x);

        SimpsonIntegrator integrator = new SimpsonIntegrator();
        double integral = integrator.integrate(maxEval, bessFunction, startOfIntegral, endOfIntegral);

        integral = integral * multiplier;

        return integral;
    }
}

class IntegralBesselFunction implements UnivariateFunction {
    private int n;
    private double x;

    public IntegralBesselFunction(int n, double x) {
        this.n = n;
        this.x = x;
    }

    public double value(double tau){
        return Math.cos(n * tau - x * Math.sin(tau));
    }
}
