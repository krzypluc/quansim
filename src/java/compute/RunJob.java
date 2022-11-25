package compute;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import mathUtils.SavitzkyGolay;
import mathUtils.chebyshev.Misc;
import miscellaneous.GroupName;
import miscellaneous.hdf5Handler;
import org.apache.commons.math3.complex.Complex;
import org.pcj.*;
import mathUtils.FFTDerivative;

import static java.lang.Math.PI;
import static mathUtils.Functions.*;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    public enum SharedRunJob {
        x, y, integral, transformedFromTheNextProcess, transformed, potential
    }

    private double[] x;
    private Complex[] y;
    private double[] potential;
    private Complex integral = Complex.ZERO;
    private Complex[] transformedFromTheNextProcess;
    private Complex[] transformed;

    private static Map<String, Object> config;
    private static Map<String, String> filePaths;
    private static Map<String, Integer> domain;
    private static Map<String, Double> constants;
    private static Map<String, Double> time;

    @Override
    public void main() {
        // Seting values from config
        // by default config is stored in config folder
        int RESOLUTION_EXPOTENTIAL = (int) domain.get("resolutionExpotential");
        int PERIOD_NUMBER = (int) domain.get("period");
        int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

        // Time
        double dt = (double) time.get("dt");
        double startTime = (double) time.get("startTime");
        double endTime = (double) time.get("endTime");

        // Constants
        double mass = (double) constants.get("mass");
        double alfa = (double) constants.get("alfa");
        double planckConstant = (double) constants.get("planckConstant");

        double period = PERIOD_NUMBER * PI;

        // One dimensional wave function: x - arguments, y - values
        x = new double[RESOLUTION];
        y = new Complex[RESOLUTION];
        integral = Complex.ZERO;
        transformedFromTheNextProcess = new Complex[0];

        // Use constants as process id and number of processes
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;
        double dx = period / (RESOLUTION - 1);

        // Calculate values of initial wave function.
        // Getting and casting all values from returned array.
        // Calulcate integral in all processes.
        HashMap<String, Object> initVal = initYandX(y, x, period, dx);

        // Casting, because initVal returns objects, not certain types.
        y = (Complex[]) initVal.get("y");
        x = (double[]) initVal.get("x");
        integral = (Complex) initVal.get("integral");
        potential = (double[]) initVal.get("potential");

        // ---------------- PARALLEL INTEGRAL
        // Collect integral and share to processes final value.
        if (procID == 0) {
            Complex coll = PCJ.reduce(
                    // Lambda operator / reductor
                    (subtotal, element) -> subtotal.add(element),
                    SharedRunJob.integral);
            // Broadcast to all processes
            PCJ.asyncBroadcast(coll, SharedRunJob.integral);
        }

        PCJ.waitFor(SharedRunJob.integral, 1);

        // Normalize wave function
        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {

            y[i] = y[i].divide(integral);
        }

        // Collect all calculations
        if (procID == 0) {
            int procIncr = 1;
            for (int i = lengthOfpiece; i < y.length; i++) {
                x[i] = PCJ.get(procIncr, SharedRunJob.x, i);
                y[i] = PCJ.get(procIncr, SharedRunJob.y, i);
                potential[i] = PCJ.get(procIncr, SharedRunJob.potential, i);

                if ((i + 1) % lengthOfpiece == 0) {
                    procIncr++;
                }
            }
            PCJ.asyncBroadcast(y, SharedRunJob.y);
            PCJ.asyncBroadcast(x, SharedRunJob.x);
            PCJ.asyncBroadcast(potential, SharedRunJob.potential);
        }

        PCJ.waitFor(SharedRunJob.y);

        double[] RandG = Misc.getChebyshevConstants(dt, potential, mass, dx);
        double R = RandG[0];
        double G = RandG[1];
        double Vmin = RandG[2];
        double deltaE = RandG[3];

        int N = (int) (((deltaE * dt) / 2) * alfa);
        //int N = 183;
        int timesteps = (int) ((endTime - startTime) / dt) + 1;

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
            ySecondDerivative = FFTDerivative.derivativeComplex(y, x);

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
                yDer = FFTDerivative.derivativeComplex(chebyshevPolynomials[i - 1], x);

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

        // Saving to HDF5
        if (procID == 0) {
            hdf5Handler.saveYandDerivative(x, yHistory, yDerHistory, timesteps, filePaths);
        }
    }

    public static void main(String[] args) throws IOException {
        String configPath;
        try {
            configPath = args[0];
        } catch (ArrayIndexOutOfBoundsException e) {
            System.out.println("You need to provide path to config.");
            throw e;
        }

        config = loadConfigFromYaml(configPath);
        filePaths = (Map<String, String>) config.get("filePaths");
        domain = (Map<String, Integer>) config.get("domain");
        constants = (Map<String, Double>) config.get("constants");
        time = (Map<String, Double>) config.get("time");

        String nodesFile = filePaths.get("nodesFile");


        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}
