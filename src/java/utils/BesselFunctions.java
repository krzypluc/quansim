package utils;

import org.apache.commons.math3.complex.Complex;
import scala.Int;

import java.io.IOException;
import java.util.Map;

import static utils.Functions.loadConfigFromYaml;

public class BesselFunctions {
    private Map<String, Object> config;
    private Map<String, Integer> constants;

    public BesselFunctions() throws IOException {
        config = loadConfigFromYaml("config/config.yml");
        constants = (Map<String, Integer>) config.get("constants");
    }

    public double Besselnx(int n, double x) {
        double startOfIntegral = 0.0;
        double endOfIntegral = Math.PI;
        double multiplier = 1 / Math.PI;

        double firstValue = 0;
        double lastValue = 0;

        int resExp = (int) constants.get("besselResolutionExpotential");
        int resolution = (int) Math.pow(2, resExp);
        double[] y = new double[resolution + 1];

        double sum = 0;
        double dtau = (endOfIntegral - startOfIntegral) / resolution;
        double tau = startOfIntegral;

        for (int i = 0; i <= resolution; i++) {
            y[i] = Math.cos(n * tau - x * Math.sin(tau));
            sum += y[i];
            tau += dtau;

            // Saving first value - use in integral
            if (i == 0) {
                firstValue = y[i];
            }

            // Saving last value - use in integral
            if (i == (resolution - 1)) {
                lastValue = y[i];
            }
        }
        double integral = ((sum * 2) - lastValue - firstValue) * (dtau / 2);
        integral = integral * multiplier;

        return integral;
    }
}
