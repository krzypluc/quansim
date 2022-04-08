package compute;

import java.io.File;
import java.util.Random;

import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;
import utils.Complex;

@RegisterStorage(RunJob.Shared.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    enum Shared { nAll, A}
    int nAll = 10480000;
    int[] A = new int[nAll];

    @Override
    public void main() {

        if (PCJ.myId() == 0){
            Random rand = new Random();

            for (int i = 0; i < nAll; i++){
                A[i] = rand.nextInt();
            }
            PCJ.broadcast(A, Shared.A);
        }

    }

    public static void main(String[] args) throws Throwable {
        String nodesFile  = "nodes.txt";
        PCJ.executionBuilder (RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}