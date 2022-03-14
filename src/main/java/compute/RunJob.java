package compute;

import java.io.File;
import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;
import utils.Complex;

@RegisterStorage(RunJob.Shared.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    enum Shared { nAll, A}
    long nAll = 1048576L;
    Complex[] A = new Complex[Math.toIntExact(nAll)];

    @Override
    public void main() {
        PCJ.barrier();

        for (int i =0; i < nAll; i++){
            A[i] = new Complex(1, 2);
        }
    }

    public static void main(String[] args) throws Throwable {
        String nodesFile  = "nodes.txt";
        PCJ.executionBuilder (RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}