package mathUtils.chebyshev;

import mathUtils.BesselFunctions;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.complex.Complex;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Misc {
    public static double[] getChebyshevConstants(double dt, double[] potential, double mass, double dx) {
        Double[] pot = ArrayUtils.toObject(potential);
        List<Double> list = Arrays.asList(pot);
        double Vmin = Collections.min(list);
        double Vmax = Collections.max(list);

        double R = dt * ((Math.pow(Math.PI, 2) / (2 * mass * Math.pow(dx, 2))) + Vmax - Vmin) / 2;
        double G = Vmin * dt;
        
        double Pmax = Math.PI / dx;
        double Emin = Vmin;
        double Emax = Vmax + Math.pow(Pmax, 2) / (2 * mass);
        double deltaE = Emax - Emin;
        
        return new double[]{R, G, Vmin, deltaE};
    }

    public static Complex getAk(int k, double R, double G) {
        // exp(i * (R + G))
        Complex exp = Complex.I.multiply(R + G).exp();

        // Ck = 1 if k == 1, else: Ck = 2
        int Ck = (k == 1 ? 1 : 2);

        // e^(i * (R + G)) * Ck * Jk(R)
        exp = exp.multiply(Ck).multiply(BesselFunctions.Besselnx(k, R));

        return exp;
    }
}
