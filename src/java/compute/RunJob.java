package compute;

import java.io.*;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.*;
import org.yaml.snakeyaml.Yaml;
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

        /*
        ------------- FFT LOOP


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

        // Initiate FFT in each process.
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

        // Calculate FFT of a fragment.
        transformed = fft.transform(yFragment, TransformType.FORWARD);
        // Complex[] transformed = DFT.forwardDFT(yFragment);

        // Fragment buffer is used for buffing a fragment of wave function previously used by a process.
        Complex[] yFragmentBuffer;

        Complex p;
        Complex q;
        int N = y.length;
        int k;

        int numberOfLoops = (int) (log(procCount) / log(2));

        int sizeOfGroup = 1;
        int expLoopCounter = 1;
        double groupIdBuffer;
        int groupId;
        int communicateWith = -1;

        for (int i = 0; i < numberOfLoops; i++){
            groupIdBuffer = (double) (procID / sizeOfGroup);
            groupId = (int) groupIdBuffer;

            // Even communicate with next processes
            if (groupId % 2 == 0) {
                communicateWith = procID + expLoopCounter;
            }

            // Odd processes communicate with previous processes.
            else {
                communicateWith = procID - expLoopCounter;
            }

            transformedFromTheNextProcess = transformed;

            PCJ.asyncPut(transformed, communicateWith, SharedRunJob.transformedFromTheNextProcess);
            PCJ.waitFor(SharedRunJob.transformedFromTheNextProcess, i);
            transformedFromTheNextProcess = (Complex[]) PCJ.get(communicateWith, SharedRunJob.transformedFromTheNextProcess);

            if (groupId % 2 == 0) {
                for (int l = 0; l < transformed.length; l++){
                    p = transformed[l];
                    q = Complex.I.multiply(-2 * PI * l / N).exp().multiply(transformedFromTheNextProcess[l]);
                    transformed[l] = p.add(q);
                }
            } else {
                for (int l = 0; l < transformed.length; l++) {
                    p = transformedFromTheNextProcess[l];
                    q = Complex.I.multiply(-2 * PI * l / N).exp().multiply(transformed[l]);
                    transformed[l] = p.subtract(q);
                }
            }

            sizeOfGroup = sizeOfGroup * 2;
            expLoopCounter = expLoopCounter * 2;
        }

        // Old fft fragment
        /*


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

        if (procID == 0){
            Complex[] buffer;
            Complex[] outcome = new Complex[y.length];
            for (int l = 1; l < procCount; l++){
                buffer = (Complex[]) PCJ.get(l, SharedRunJob.transformed);

                j = 0;
                for (int m = (l * (y.length / procCount)); m < ((l + 1) * (y.length / procCount)); m++){
                    outcome[m] = buffer[j];
                    j++;
                }
            }

            Complex[] transformedCheck = fft.transform(y, TransformType.FORWARD);

            for (Complex nbr : transformedCheck){
                System.out.println(nbr);
            }
            System.out.println("-----");
            for (Complex nbr : outcome){
                System.out.println(nbr);
            }

            */

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

        Complex[] y_transformed = fft.transform(y, TransformType.FORWARD);
        double[] freq = FTUtils.freq(x);

        if (procID == 0){
            System.out.println(time.get("endTime"));
        }

//        if (procID == 0) {
//            String groupName = LocalDateTime.now().toString();
//            IHDF5Writer writer = HDF5Factory.open("hdf5/myfile.h5");
//            writer.writeDoubleArray(groupName + "/x", x);
//
//            double[][] yTransformedDouble = new double[y.length][2];
//            for (int i = 0; i < y.length; i++){
//                yTransformedDouble[i][0] = y_transformed[i].getReal();
//                yTransformedDouble[i][1] = y_transformed[i].getImaginary();
//            }
//
//            writer.writeDoubleMatrix(groupName + "/y", yTransformedDouble);
//            writer.close();
//        }
    }

    public static void main(String[] args) throws IOException {
        String nodesFile = (String) loadConfigFromYaml("config/config.yml").get("nodesFileName");
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}