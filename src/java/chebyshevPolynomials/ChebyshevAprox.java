package chebyshevPolynomials;

import mathUtils.DerivativeFFT;
import mathUtils.SavitzkyGolay;
import org.apache.commons.math3.complex.Complex;

public class ChebyshevAprox {
    public static Complex[][][] aproximate(double[] x, Complex[] y, double[] potential, int timesteps, double dt, double mass, double planckConstant, double alfa, double dx){

        double[] RandG = Misc.getChebyshevConstants(dt, potential, mass, dx);
        double R = RandG[0];
        double G = RandG[1];
        double Vmin = RandG[2];
        double deltaE = RandG[3];

        int N = (int) (((deltaE * dt) / 2) * alfa);
        //int N = 183;

        // Hamiltonian
        double momentumConst = (-1) * Math.pow(planckConstant, 2) / (2 * mass);

        // Polynomial matrix
        Complex[][] chebyshevPolynomials = new Complex[N][x.length];

        // Matrix containing history of timesteps
        Complex[][] yHistory = new Complex[timesteps][y.length];
        Complex[][] yDerHistory = new Complex[timesteps][y.length];
        yHistory[0] = y;

        // Get Bessel values
        Complex a0 = Misc.getAk(0, R, G);
        Complex a1 = Misc.getAk(1, R, G);

        // Matrix containing previous Chebyshev polynomials

        Complex[] ySecondDerivative = new Complex[y.length];
        Complex momentum;
        Complex potentialChebPart;
        double normalizationPart;
        Complex momentumNormalizationPart;

        Complex[] yDer;
        Complex[] chebPrevPrev;
        Complex[] chebPrev;
        Complex ak;
        Complex[] sumOfChebPolynomials = new Complex[y.length];
        Complex chebIntegral = Complex.ZERO;

        Complex value1 = Complex.ZERO;
        Complex value2 = Complex.ZERO;
        Complex sumator = Complex.ZERO;

        Complex energyIntegral = Complex.ZERO;
        Complex[] energyCalc = new Complex[y.length];
        // Interating over timesteps
        // timesteps - 1: because in every cicle we are calculating wave function for the next step.

        for (int h = 0; h < (timesteps - 1); h++) {

            y = yHistory[h];
            ySecondDerivative = DerivativeFFT.derivativeComplex(y, x);

            for (int i = 0; i < sumOfChebPolynomials.length; i++) {
                sumOfChebPolynomials[i] = Complex.ZERO;
            }

            for (int i = 0; i < chebyshevPolynomials[0].length; i++) {
                // 1 * y * a0
                chebyshevPolynomials[0][i] = y[i].multiply(a0);

                // Add to the sum
                sumOfChebPolynomials[i] = sumOfChebPolynomials[i].add(chebyshevPolynomials[0][i]);

                // Momentum - (-h^2/2m) * y'')
                momentumNormalizationPart = Complex.I.multiply(-1);
                momentum = ySecondDerivative[i].multiply(momentumConst).multiply(momentumNormalizationPart);

                // potential * phi
                potentialChebPart = y[i].multiply(potential[i]);

                // Norm part - (deltaE * y[i]) / 2 + Vmin * y[i]
                // normalizationPart = y[i].multiply(( - deltaE / 2) + Vmin);
                normalizationPart = dt / R;

                // Numerator -
                chebyshevPolynomials[1][i] = momentum.add(potentialChebPart).multiply(normalizationPart);
                energyCalc[i] = momentum.add(potentialChebPart);

                // Denominator
                //chebyshevPolynomials[1][i] = chebyshevPolynomials[1][i].divide(y[i].multiply(deltaE));

                // Multiply by 2
                //chebyshevPolynomials[1][i] = chebyshevPolynomials[1][i].multiply(2);

                // Multiply by a1
                //chebyshevPolynomials[1][i] = chebyshevPolynomials[1][i].multiply(a1);

                // Add to the sum
                sumOfChebPolynomials[i] = sumOfChebPolynomials[i].add(chebyshevPolynomials[1][i].multiply(a1));
            }

            // Calulcate energy
//            energyIntegral = Complex.ZERO;
//            for (int i = 0; i < energyCalc.length - 1; i++) {
//                energyIntegral = energyIntegral.add(
//                        energyCalc[i].multiply(y[i].conjugate())
//                                .add(
//                                        energyCalc[i + 1].multiply(y[i + 1].conjugate())
//                                )
//                                .multiply(dx)
//                                .divide(2)
//                );
//            }
//            System.out.println(energyIntegral);

            for (int i = 2; i < N; i++) {
                yDer = DerivativeFFT.derivativeComplex(chebyshevPolynomials[i - 1], x);

                // Savitzky-Golay Filter
                for (int j = 0; j < yDer.length; j++) {
                    yDer[j] = SavitzkyGolay.smoothWindow7(yDer, j);
                }

                chebPrevPrev = chebyshevPolynomials[i - 2];
                chebPrev = chebyshevPolynomials[i - 1];
                ak = Misc.getAk(i, R, G);

                for (int j = 0; j < chebyshevPolynomials[i].length; j++) {
                    // --- Calculate (-2i) * Hnorm * cheb[i-1]
                    // Momentum normalization part is -i
                    // Momentum - (-h^2/2m) * y'')

                    momentumNormalizationPart = Complex.I.multiply(-1);
                    momentum = yDer[j].multiply(momentumConst).multiply(momentumNormalizationPart);

                    // potential * phi
                    potentialChebPart = chebPrev[j].multiply(potential[j]);

                    // Norm part - (deltaE * y[i]) / 2 + Vmin * y[i]
                    // normalizationPart = chebPrev[j].multiply((deltaE / 2) + Vmin);

                    // Normalisation part
                    normalizationPart = dt / R;

                    // Numerator
                    chebyshevPolynomials[i][j] = momentum.add(potentialChebPart).multiply(normalizationPart);

                    // Denominator
                    // chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].divide(chebPrev[j].multiply(deltaE));

                    // Multiply by 2
                    chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].multiply(2);

                    // Multiply by -2i
                    // chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].multiply(Complex.I).multiply(-2);

                    // --- Add cheb[i - 2]
                    chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].add(chebPrevPrev[j]);

                    // Multiply by ak
                    //chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j];

                    // Add polynomial
                    sumOfChebPolynomials[j] = sumOfChebPolynomials[j].add(chebyshevPolynomials[i][j].multiply(ak));
                }
            }

            // Smooth with Savitzky-Golay filter
            for (int i = 0; i < sumOfChebPolynomials.length; i++) {
                sumOfChebPolynomials[i] = SavitzkyGolay.smoothWindow7(sumOfChebPolynomials, i);
            }

            chebIntegral = Complex.ZERO;
            // Integrate
            for (int i = 0; i < sumOfChebPolynomials.length - 1; i++) {
                // y1 ** 2
                value1 = sumOfChebPolynomials[i].multiply(sumOfChebPolynomials[i]);

                // y2 ** 2
                value2 = sumOfChebPolynomials[i + 1].multiply(sumOfChebPolynomials[i + 1]);

                // Calculate area
                sumator = value1.add(value2).multiply(dx).divide(2);

                // Add sumator
                chebIntegral = chebIntegral.add(sumator);
            }

            // Normalize
            for (int i = 0; i < sumOfChebPolynomials.length; i++) {
                sumOfChebPolynomials[i] = sumOfChebPolynomials[i].divide(chebIntegral);
            }

            // Saving history
            for (int i = 0; i < y.length; i++) {
                yHistory[h + 1][i] = sumOfChebPolynomials[i];
                yDerHistory[h][i] = ySecondDerivative[i];
            }

            // for the last timestep
            if (h == (timesteps - 2)){
                yDerHistory[h + 1] = ySecondDerivative;
            }

        }

        Complex[][][] ret = new Complex[2][][];
        ret[0] = yHistory;
        ret[1] = yDerHistory;

        return ret;
    }
}
