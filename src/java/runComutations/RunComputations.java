package runComutations;

import java.io.*;
import java.util.Map;

import initializeWaveFunction.InitializeWaveFunction;
import initializeWaveFunction.InitializeWaveFunctionStorage;
import mathUtils.Potential;
import mathUtils.WaveFunctions;
import miscellaneous.YamlLoader;
import parallelFFT.ParallelFFT;
import parallelIntegrator.ParalellComplexIntegratorStorage;
import parallelIntegrator.ParallelComplexIntegrator;
import chebyshevPolynomials.ChebyshevAprox;
import miscellaneous.HDF5Handler;
import org.apache.commons.math3.complex.Complex;
import org.pcj.*;

import static java.lang.Math.PI;

@RegisterStorage(RunComputations.SharedRunJob.class)
public class RunComputations implements StartPoint {

    @Storage(RunComputations.class)
    public enum SharedRunJob {
        // RunComputations storage
        x,
        y,
        integral,
        potential,

        // ParallelComplexIntegrator storage
        parallelIntegralStorage,

        // InitializeWaveFunction storage
        initWFStorage;
    }

    // SharedRunJob storage init
    double[] x;
    Complex[] y;
    double[] potential;
    Complex integral;

    // ParallelComplexIntegrator storage
    ParalellComplexIntegratorStorage parallelIntegralStorage;

    // InitializeWaveFunction storage
    InitializeWaveFunctionStorage initWFStorage;

    static Map<String, Object> config;
    static Map<String, String> filePaths;
    static Map<String, Integer> domain;
    static Map<String, Double> constants;
    static Map<String, Double> time;
    static Map<String, String> functions;

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

        // Wave function and potential load
        WaveFunctions.WaveFunction waveFunctionType = WaveFunctions.WaveFunction.valueOf(functions.get("waveFunction"));
        Potential.PotentialType potentialType = Potential.PotentialType.valueOf(functions.get("potential"));

        // Calculating time array
        double[] time = new double[timesteps];
        time[0] = startTime;
        for (int i = 1; i < timesteps; i++){
            time[i] = time[i - 1] + dt;
        }

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

        //ParallelComplexIntegrator storage init

        // Use constants as process id and number of processes
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;
        double dx = period / (RESOLUTION - 1);

        // Calculate values of initial wave function.
        // Getting and casting all values from returned array.
        // Calulcate integral in all processes.
        InitializeWaveFunction waveFunctionInit = new InitializeWaveFunction(y, x, potential, period, dx, waveFunctionType, potentialType);
        waveFunctionInit.main();

        // Wait for calculations
        PCJ.waitFor(SharedRunJob.initWFStorage);

        // Update values from initialization of wave function
        y = initWFStorage.getY();
        x = initWFStorage.getX();
        potential = initWFStorage.getPotential();

        // Integrate parallely - putting value on "itegral"
        ParallelComplexIntegrator parallelIntegrator = new ParallelComplexIntegrator(y, x);
        parallelIntegrator.main();

        PCJ.waitFor(SharedRunJob.parallelIntegralStorage, 2);

        // Update values from integral
        integral = parallelIntegralStorage.getIntegralValue();

        // Normalize wave function
        for (int i = 0; i < y.length; i++) {
            y[i] = y[i].divide(integral);
        }

        // Not working yet - testing purpose
        ParallelFFT fft = new ParallelFFT(y, x);
        fft.main();

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
            HDF5Handler.saveYandDerivative(x, yHistory, yDerHistory, potential, time, timesteps, filePaths);
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

        config = YamlLoader.loadConfigFromYaml(configPath);
        filePaths = (Map<String, String>) config.get("filePaths");
        domain = (Map<String, Integer>) config.get("domain");
        constants = (Map<String, Double>) config.get("constants");
        time = (Map<String, Double>) config.get("time");
        functions = (Map<String, String>) config.get("functions");

        String nodesFile = filePaths.get("nodesFile");
        
        PCJ.executionBuilder(RunComputations.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}
