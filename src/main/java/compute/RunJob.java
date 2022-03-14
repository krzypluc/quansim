package compute;

import java.io.File;
import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;

/*
author: Krzysztof Płuciennik
*/

@RegisterStorage(RunJob.Shared.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    enum Shared { nAll, A, B, C, startTime}
    int nAll = 1048576;
    double[] A = new double[nAll];
    double[] B = new double[nAll];
    double[] C = new double[nAll];
    long startTime;

    @Override
    public void main() {
        PCJ.barrier();



    }

    public static void main(String[] args) throws Throwable {
        String nodesFile  = "nodes.txt";
        PCJ.executionBuilder (RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}