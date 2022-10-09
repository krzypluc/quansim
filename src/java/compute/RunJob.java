package compute;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.*;

import static java.lang.Math.PI;
import static utils.functions.initYandX;
import static utils.functions.waveFunction;
import static utils.integrator.TrapezoidComplex1D;
import static utils.parallelFFT.FFTComplex1D;

// import org.apache.commons.numbers.complex.Complex;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 6;
    static final int PERIOD_EXPOTENTIAL = 6;
    static final int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

    static final double period = PERIOD_EXPOTENTIAL * PI;

    @Storage(RunJob.class)
    public enum SharedRunJob {
        x, y, integral, nextYFragment
    }

    // One dimensional wave function: x - arguments, y - values
    private double[] x = new double[RESOLUTION];
    private Complex[] y = new Complex[RESOLUTION];
    private Complex integral = Complex.ZERO;
    private ArrayList<Complex> nextYFragment = new ArrayList<Complex>();

    @Override
    public void main() throws Throwable {
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;
        double dx = period / (RESOLUTION - 1);

        // Calculate values of initial wave function.
        // Getting and casting all values from returned array.
        Object[] initVal = initYandX(y, x, period, dx);
        y = (Complex[]) initVal[0];
        x = (double[]) initVal[1];
        integral = (Complex) initVal[2];

        // Calculate integral and share to processes
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
            PCJ.broadcast(x, SharedRunJob.x);
            PCJ.broadcast(y, SharedRunJob.y);
        }

        // FFT
        PCJ.waitFor(SharedRunJob.y);
        Complex[] transformed = FFTComplex1D(y);

    }

    public static void main(String[] args) throws IOException {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}