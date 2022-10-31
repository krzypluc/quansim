package compute;

import java.io.*;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.pcj.*;
import utils.BesselFunctions;
import utils.FFTDerivative;
import utils.Integrator;

import static java.lang.Math.PI;
import static utils.Functions.initYandX;
import static utils.Functions.loadConfigFromYaml;

// import org.apache.commons.numbers.complex.Complex;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    public enum SharedRunJob {
        x, y, integral, transformedFromTheNextProcess, transformed
    }

    private double[] x;
    private Complex[] y;
    private Complex integral = Complex.ZERO;
    private Complex[] transformedFromTheNextProcess;
    private Complex[] transformed;

    private static Map<String, Object> config;

    @Override
    public void main() throws IOException {
        // Set values form config
        Map<String, Double> time = (Map<String, Double>) config.get("time");
        Map<String, Integer> domain = (Map<String, Integer>) config.get("domain");

        // Seting values from config
        int RESOLUTION_EXPOTENTIAL = (int) domain.get("resolutionExpotential");
        int PERIOD_NUMBER = (int) domain.get("period");
        int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

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
        Object[] initVal = initYandX(y, x, period, dx);

        // Casting, because initVal returns objects, not certain types.
        y = (Complex[]) initVal[0];
        x = (double[]) initVal[1];
        integral = (Complex) initVal[2];
        
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

        PCJ.waitFor(SharedRunJob.integral);
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

                if ((i + 1) % lengthOfpiece == 0) {
                    procIncr++;
                }
            }
            PCJ.asyncBroadcast(y, SharedRunJob.y);
            PCJ.asyncBroadcast(x, SharedRunJob.x);
        }
        PCJ.waitFor(SharedRunJob.y);

        Complex[] y_derivative = FFTDerivative.derivativeComplex(y, x);

        if (procID == 0) {

            double jk = BesselFunctions.Besselnx(1, 2);

            System.out.println(jk);
        }

    }

    public static void main(String[] args) throws IOException {
        config = loadConfigFromYaml("config/config.yml");
        String nodesFile = (String) config.get("nodesFileName");
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}