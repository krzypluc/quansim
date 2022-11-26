package compute;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import mathUtils.InitializeWaveFunction;
import mathUtils.ParallelIntegrator;
import mathUtils.chebyshev.ChebyshevAprox;
import miscellaneous.HDF5Handler;
import org.apache.commons.math3.complex.Complex;
import org.pcj.*;

import static java.lang.Math.PI;
import static mathUtils.Functions.*;


@SuppressWarnings("unused")
@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    public enum SharedRunJob {
        x, y, integral, potential
    }

    private double[] x;
    private Complex[] y;
    private double[] potential;
    private Complex integral = Complex.ZERO;

    private static Map<String, Object> config;
    private static Map<String, String> filePaths;
    private static Map<String, Integer> domain;
    private static Map<String, Double> constants;
    private static Map<String, Double> time;

    @Override
    public void main() throws IOException {
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
        potential = new double[RESOLUTION];
        integral = Complex.ZERO;

        // Use constants as process id and number of processes
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;
        double dx = period / (RESOLUTION - 1);

        // Calculate values of initial wave function.
        // Getting and casting all values from returned array.
        // Calulcate integral in all processes.
        InitializeWaveFunction waveFunctionInit = new InitializeWaveFunction(y, x, period, dx);
        waveFunctionInit.main();

        // Integrate parallely - putting value on "itegral"
        ParallelIntegrator parallelIntegrator = new ParallelIntegrator();
        parallelIntegrator.main();

        PCJ.waitFor(SharedRunJob.y);

        // Normalize wave function
        for (int i = 0; i < y.length; i++) {
            y[i] = y[i].divide(integral);
        }

        // Chebyshev aproxymation
        Complex[][][] chebOutput = ChebyshevAprox.aproximate(
                x,
                y,
                potential,
                timesteps,
                dt,
                mass,
                planckConstant,
                alfa,
                dx
        );

        Complex[][] yHistory = chebOutput[0];
        Complex[][] yDerHistory = chebOutput[1];

        // Saving to HDF5
        if (procID == 0) {
            HDF5Handler.saveYandDerivative(x, yHistory, yDerHistory, timesteps, filePaths);
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
