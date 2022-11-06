package compute;

import java.io.*;
import java.time.LocalDateTime;
import java.util.Map;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import mathUtils.chebyshev.Misc;
import miscellaneous.GroupName;
import org.apache.commons.math3.complex.Complex;
import org.pcj.*;
import mathUtils.FFTDerivative;

import static java.lang.Math.PI;
import static mathUtils.Functions.initYandX;
import static mathUtils.Functions.loadConfigFromYaml;


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
    private static Map<String, String> filePaths;
    private static Map<String, Integer> domain;
    private static Map<String, Double> constants;
    private static Map<String, Double> time;

    @Override
    public void main() throws IOException {
        // Seting values from config
        int RESOLUTION_EXPOTENTIAL = (int) domain.get("resolutionExpotential");
        int PERIOD_NUMBER = (int) domain.get("period");
        int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

        double dt = (double) time.get("dt");

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

                if ((i + 1) % lengthOfpiece == 0) {
                    procIncr++;
                }
            }
            PCJ.asyncBroadcast(y, SharedRunJob.y);
            PCJ.asyncBroadcast(x, SharedRunJob.x);
        }

        PCJ.waitFor(SharedRunJob.y);

        Complex[] y_derivative = FFTDerivative.derivativeComplex(y, x);

        // Saving tp HDF5
        if (procID == 0) {

            System.out.println(Misc.getR(dt, potential, mass));

            String groupName = GroupName.getGroupName();
            String hdf5FileName = (String) filePaths.get("hdf5JavaFile");

            IHDF5Writer writer = HDF5Factory.open(hdf5FileName);
            writer.writeDoubleArray(groupName + "/x", x);

            double[][] yTransformedDouble = new double[y.length][2];
            for (int i = 0; i < y.length; i++) {
                yTransformedDouble[i][0] = y_derivative[i].getReal();
                yTransformedDouble[i][1] = y_derivative[i].getImaginary();
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
        config = loadConfigFromYaml("config/config.yml");
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