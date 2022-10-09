package compute;

import java.io.*;

import org.pcj.*;

import static utils.functions.waveFunction;


@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 6;
    static final int PERIOD_EXPOTENTIAL = 6;
    static final int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

    static final double PI = Math.PI;
    static final double period = PERIOD_EXPOTENTIAL * Math.PI;

    @Storage(RunJob.class)
    private enum SharedRunJob {
        x, y
    }

    private double[] x = new double[RESOLUTION];
    private Complex[] y = new Complex[RESOLUTION];

    @Override
    public void main() throws Throwable {
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        // Length of the array per process
        int lengthOfpiece = (int) y.length / procCount;

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++){
            x[i] = (-period / 2) + period / (RESOLUTION - 1) * i;
            y[i] = waveFunction(new Complex(x[i], 0));
        }

        if (procID == 0){
            int procIncr = 1;
            for (int i = lengthOfpiece; i < y.length; i++){
                y[i] = PCJ.get(procIncr, SharedRunJob.y, i);

                if ((i + 1) % lengthOfpiece == 0){
                    procIncr++;
                }
            }

            for (Complex cmp : y){
                System.out.println(cmp);
            }
        }

    }

        public static void main (String[]args) throws IOException {
            String nodesFile = "nodes.txt";
            PCJ.executionBuilder(RunJob.class)
                    .addNodes(new File(nodesFile))
                    .start();
        }
    }