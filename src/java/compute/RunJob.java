package compute;

import java.io.*;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.*;
import org.yaml.snakeyaml.Yaml;
import utils.FFTDerivative;
import utils.FTUtils;

import static java.lang.Math.PI;
import static utils.functions.initYandX;
import static utils.functions.loadConfigFromYaml;

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

    @Override
    public void main() throws IOException {
        // Load config
        Map<String, Object> config = loadConfigFromYaml("config/config.yml");

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

        // Collect integral and share to processes final value.
        if (procID == 0) {
            Complex coll = PCJ.reduce(
                    // Lambda operator / reductor
                    (subtotal, element) -> subtotal.add(element),
                    SharedRunJob.integral);
            // Broadcast to all processes
            PCJ.broadcast(coll, SharedRunJob.integral);
        }

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
            PCJ.asyncBroadcast(x, SharedRunJob.x);
            PCJ.asyncBroadcast(y, SharedRunJob.y);
        }
        PCJ.waitFor(SharedRunJob.y);
        

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

        Complex[] y_derivative = FFTDerivative.derivativeComplex(y, x);

        if (procID == 0){
            System.out.println(y.length);
        }


    }

    public static void main(String[] args) throws IOException {
        String nodesFile = (String) loadConfigFromYaml("config/config.yml").get("nodesFileName");
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}