package compute;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import mathUtils.chebyshev.ChebyshevAprox;
import mathUtils.chebyshev.Misc;
import miscellaneous.hdf5Handler;
import org.apache.commons.math3.complex.Complex;
import org.pcj.*;

import static java.lang.Math.PI;
import static mathUtils.Functions.*;


@SuppressWarnings("unused")
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
        int timesteps = (int) ((endTime - startTime) / dt) + 1;

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

        // Chebyshev aproxymation
        Complex[][][] chebOutput = ChebyshevAprox.aproximate(x, y, timesteps, dt, potential, mass, planckConstant, alfa, dx);
        Complex[][] yHistory = chebOutput[0];
        Complex[][] yDerHistory = chebOutput[1];

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
