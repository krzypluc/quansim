package compute;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import mathUtils.chebyshev.Misc;
import miscellaneous.GroupName;
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
        int RESOLUTION_EXPOTENTIAL = (int) domain.get("resolutionExpotential");
        int PERIOD_NUMBER = (int) domain.get("period");
        int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

        // Time
        double dt = (double) time.get("dt");

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

        int N = (int) (R * alfa);

        // Hamiltonian
        double momentumConst = (-1) * Math.pow(planckConstant, 2) / (2 * mass);

        // Polynomial 0, and 1 - or n-2 and n-1
        Complex[][] chebyshevPolynomials = new Complex[N][x.length];

        // k = n + 1
        Complex a0 = Misc.getAk(0, deltaE, Vmin, dt);
        Complex a1 = Misc.getAk(1, deltaE, Vmin, dt);

        Complex[] ySecondDerivative = FFTDerivative.derivativeComplex(y, x);
        Complex momentum;
        Complex potentialChebPart;
        Complex normalizationPart;
        for (int i = 0; i < chebyshevPolynomials[0].length; i++) {
            // 1 * y * a0
            chebyshevPolynomials[0][i] = y[i].multiply(a0);

            // Momentum - (-h^2/2m) * y'')
            momentum = ySecondDerivative[i].multiply(momentumConst);

            // potential * phi
            potentialChebPart = y[i].multiply(potential[i]);

            // Norm part - (deltaE * y[i]) / 2 + Vmin * y[i]
            normalizationPart = y[i].multiply((deltaE / 2) + Vmin);

            // Numerator -
            chebyshevPolynomials[1][i] = momentum.add(potentialChebPart).add(normalizationPart);

            // Denominator
            chebyshevPolynomials[1][i] = chebyshevPolynomials[1][i].divide(y[i].multiply(deltaE));

            // Multiply by 2 * ak
            chebyshevPolynomials[1][i] = chebyshevPolynomials[1][i].multiply(2).multiply(a1);
        }

        Complex[] yDer;
        Complex[] cheb2;
        Complex[] cheb1;
        for (int i = 2; i < N; i++) {
            yDer = FFTDerivative.derivativeComplex(chebyshevPolynomials[i - 1], x);
            cheb2 = chebyshevPolynomials[i - 2];
            cheb1 = chebyshevPolynomials[i - 1];

            for (int j = 0; j < chebyshevPolynomials[i].length; j++) {
                // --- Calculate (-2i) * Hnorm * cheb[i-1]
                // Momentum - (-h^2/2m) * y'')
                momentum = yDer[i].multiply(momentumConst);

                // potential * phi
                potentialChebPart = cheb1[j].multiply(potential[i]);

                // Norm part - (deltaE * y[i]) / 2 + Vmin * y[i]
                normalizationPart = cheb1[j].multiply((deltaE / 2) + Vmin);

                // Numerator -
                chebyshevPolynomials[i][j] = momentum.add(potentialChebPart).add(normalizationPart);

                // Denominator
                chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].divide(cheb1[j].multiply(deltaE));

                // Multiply by -2i
                chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j].multiply(Complex.I).multiply(-2);

                // --- Add cheb[i - 2]
                chebyshevPolynomials[i][j] = chebyshevPolynomials[i][j] + cheb2[j];
            }
        }

        // Saving to HDF5
        if (procID == 0) {
            String groupName = GroupName.getGroupName();
            String hdf5FileName = (String) filePaths.get("hdf5JavaFile");

            IHDF5Writer writer = HDF5Factory.open(hdf5FileName);
            writer.writeDoubleArray(groupName + "/x", x);

            double[][] yTransformedDouble = new double[y.length][2];
            for (int i = 0; i < y.length; i++) {
                yTransformedDouble[i][0] = ySecondDerivative[i].getReal();
                yTransformedDouble[i][1] = ySecondDerivative[i].getImaginary();
            }

            double[][] yDouble = new double[y.length][2];
            for (int i = 0; i < y.length; i++) {
                yDouble[i][0] = y[i].getReal();
                yDouble[i][1] = y[i].getImaginary();
            }

            writer.writeDoubleMatrix(groupName + "/y_double_derr", yTransformedDouble);
            writer.writeDoubleMatrix(groupName + "/y", yDouble);

            writer.close();
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