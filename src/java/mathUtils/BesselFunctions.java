package mathUtils;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;

public class BesselFunctions {
    private static final int maxEval = (int) Math.pow(2, 36);

    public static double Besselnx(int n, double x) {
        double startOfIntegral = 0.0;
        double endOfIntegral = Math.PI;
        double multiplier = 1 / Math.PI;

        IntegralBesselFunction bessFunction = new IntegralBesselFunction(n, x);

        SimpsonIntegrator integrator = new SimpsonIntegrator(
                0.01,
                0.0001,
                SimpsonIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
        double integral = integrator.integrate(maxEval, bessFunction, startOfIntegral, endOfIntegral);
        integral = integral * multiplier;

        return integral;
    }
}

class IntegralBesselFunction implements UnivariateFunction {
    private final int n;
    private final double x;

    public IntegralBesselFunction(int n, double x) {
        this.n = n;
        this.x = x;
    }

    public double value(double tau) {
        return Math.cos(n * tau - x * Math.sin(tau));
    }
}
