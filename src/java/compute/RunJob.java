package compute;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.*;
import utils.DFT;

import static compute.RunJob.SharedRunJob.nextYFragment;
import static java.lang.Math.PI;
import static utils.functions.initYandX;
import static utils.functions.waveFunction;
import static utils.integrator.TrapezoidComplex1D;

// import org.apache.commons.numbers.complex.Complex;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 5;
    static final int PERIOD_EXPOTENTIAL = 6;
    public static final int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

    static final double period = PERIOD_EXPOTENTIAL * PI;

    @Storage(RunJob.class)
    public enum SharedRunJob {
        x, y, integral, nextYFragment
    }

    // One dimensional wave function: x - arguments, y - values
    private double[] x = new double[RESOLUTION];
    private Complex[] y = new Complex[RESOLUTION];
    private Complex integral = Complex.ZERO;
    public Complex[] nextYFragment;

    @Override
    public void main() throws Throwable {
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


        /*
        ------------- FFT LOOP
         */

        // Fragment of wave function which is handled by a process.
        Complex[] yFragment = new Complex[lengthOfpiece];

        // Distribute fragments of wave funtion to processes.
        int j = 0;
        for (int i = 0; i < y.length; i++) {
            // If process is in first half
            if (procID < procCount / 2){
                if ((i % procCount) == procID * 2) {
                    yFragment[j] = y[i];
                    j++;
                }
            }

            // If process is in the second half
            else {
                if ((i % procCount) == (procID - procCount / 2) * 2 + 1) {
                    yFragment[j] = y[i];
                    j++;
                }
            }

        }
        // Hash set with numbers of processes. Hash set is used for speed.
        HashSet<Integer> processesNumbers = new HashSet<Integer>();

        // Add a number of processes to hash set.
        for (int i = 0; i < procCount; i++) {
            processesNumbers.add(i);
        }

        // Initiate FFT in each process.
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

        // Calculate FFT of a fragment.
        Complex[] transformed = fft.transform(yFragment, TransformType.FORWARD);
        // Complex[] transformed = DFT.forwardDFT(yFragment);

        // Fragment buffer is used for buffing a fragment of wave function previously used by a process.
        Complex[] yFragmentBuffer;

        // "Worked previously" value is used to determine, if the certain process was working in the previous loop.
        // Due to the fact, that FFT is not stationary, only in the first loop, all of the processes are working.
        // The next loop is always working with a 2 times less processes due to the reducing.
        boolean workedPreviously = true;

        // Exp loop number is a numbers of 2^n, where n is the number of the current loop.
        // Is is used to determine to which process should current process transfer the data.
        int expLoopNumber = 1;

        // Buffer used in the
        HashSet<Integer> processesNumbersBuffer = new HashSet<Integer>();
        Complex p;
        Complex q;
        int N;
        int k;

        // While is checking if the number of processes currently working is more than 1.
        while (processesNumbers.size() > 1) {
            // First step is to remove the processes that won't be working in the next loop.
            // In this case it's always a not even process.
            // So in the firs loop processes with the ID 1, 3, 5 ... will be turned off.
            processesNumbersBuffer = new HashSet<Integer>(processesNumbers);
            j = 0;
            for (int proc : processesNumbersBuffer) {
                if ((j % 2) == 1) {
                    processesNumbers.remove(proc);
                }
                j++;
            }

            // If process is working in this loop, it is waiting for the next process, to transfer the data.
            if (processesNumbers.contains(procID)) {
                nextYFragment = transformed;
                // Process is waiting for the data.
                // PCJ.waitFor(SharedRunJob.nextYFragment);
                nextYFragment = (Complex[]) PCJ.get(procID + expLoopNumber, SharedRunJob.nextYFragment);

                // Creating a new array with twice of the previous length.
                // Process will save its own data and concatenate a data from the next working process.
                yFragmentBuffer = new Complex[transformed.length * 2];

                // Elements from the even process are added as a even element in the new array.
                for (int i = 0; i < transformed.length; i++) {
                    p = transformed[i];
                    q = Complex.I.multiply(-2 * PI * i / (transformed.length * 2)).multiply(nextYFragment[i]);

                    yFragmentBuffer[i] = p.add(q);
                    yFragmentBuffer[i + transformed.length] = p.subtract(q);
                }

                // The new transformed function array
                transformed = yFragmentBuffer;
            }

            else if (!processesNumbers.contains(procID) && workedPreviously) {
                workedPreviously = false;
                nextYFragment = transformed;
                // PCJ.put(nextYFragment, procID - expLoopNumber, SharedRunJob.nextYFragment);
                System.out.println(procID + " ----> " + (procID - expLoopNumber));
            }

            else {
                break;
            }

            expLoopNumber = expLoopNumber * 2;
        }

        PCJ.barrier();

        FastFourierTransformer fftCheck = new FastFourierTransformer(DftNormalization.STANDARD);
        // Calculate FFT of a fragment.
        Complex[] transformedCheck = fftCheck.transform(y, TransformType.FORWARD);

        if(procID == 0){
            for(int i = 0; i < transformedCheck.length; i++){
                System.out.println("Real: " + transformedCheck[i]);
                System.out.println("Parallel: " + transformed[i]);
            }
        }

        PCJ.barrier();
    }

    public static void main(String[] args) throws IOException {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}