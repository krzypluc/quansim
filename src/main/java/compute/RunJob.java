package compute;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;
import utils.DFT;
import utils.FFT;
import utils.Utilities;

@RegisterStorage(RunJob.SharedRunJob.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 6;
    static final int PERIOD_EXPOTENTIAL = 3;
    static final int RESOLUTION = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);

    static double PI = Math.PI;
    static double period = PERIOD_EXPOTENTIAL * Math.PI;

    @Storage(RunJob.class)
    protected
    enum SharedRunJob {
        x
    }

    private int[] x;

    @Override
    public void main() throws Throwable {
        System.out.println("Hello there!");
    }

    public static void main(String[] args) throws IOException {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}