package compute;

import java.io.*;

import org.pcj.*;

import static java.lang.Math.PI;
import static utils.functions.waveFunction;
import static utils.integrator.TrapezoidComplex1D;

import org.apache.commons.numbers.complex.Complex;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 6;
    static final int PERIOD_EXPOTENTIAL = 6;
    static final int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

    static final double period = PERIOD_EXPOTENTIAL * PI;

    @Storage(RunJob.class)
    protected enum SharedRunJob {
        x, y, integral
    }

    // One dimensional wave function: x - arguments, y - values
    private double[] x = new double[RESOLUTION];
    private Complex[] y = new Complex[RESOLUTION];
    private Complex integral = Complex.ZERO;

    @Override
    public void main() throws Throwable {
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;
        double dx = period / (RESOLUTION - 1);

        // Calculate values of initial wave function
        Complex sumOfValues = Complex.ZERO;
        Complex firstValue = Complex.ZERO;
        Complex lastValue = Complex.ZERO;

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            x[i] = (-period / 2) + dx * i;
            y[i] = waveFunction(Complex.ofCartesian(x[i], 0));
            sumOfValues = sumOfValues.add(y[i]);

            // Saving first value - use in integral
            if (i == procID * lengthOfpiece) {
                firstValue = y[i];
            }

            // Saving last value - use in integral
            if (i == (procID + 1) * lengthOfpiece - 1) {
                lastValue = y[i];
            }
        }

        // Caculate integral
        integral = sumOfValues
                .multiply(2.0)
                .subtract(lastValue)
                .subtract(firstValue)
                .multiply(dx / 2);

        // Calculate integral and share to processes
        if (procID == 0) {
            Complex coll = PCJ.reduce(
                    // Lambda operator / reductor
                    (subtotal, element) -> subtotal.add(element),
                    SharedRunJob.integral);

            PCJ.broadcast(coll, SharedRunJob.integral);
        }

        // Normalize wave function
        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            y[i] = y[i].divide(integral);
        }

        // Collect calculations
        if (procID == 0) {

            int procIncr = 1;
            for (int i = lengthOfpiece; i < y.length; i++) {
                x[i] = PCJ.get(procIncr, SharedRunJob.x, i);
                y[i] = PCJ.get(procIncr, SharedRunJob.y, i);

                if ((i + 1) % lengthOfpiece == 0) {
                    procIncr++;
                }
            }
        }

    }

    public static void main(String[] args) throws IOException {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}