package compute;

import java.io.File;

import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;

@RegisterStorage(RunJob.Shared.class)
public class RunJob implements StartPoint {

    @Storage(RunJob.class)
    enum Shared {nAll, threadCount, array}

    int nAll = 10480000;
    int threadCount =  PCJ.threadCount();
    Complex array;

    @Override
    public void main() {
        double PI = Math.PI;

        if (PCJ.myId() == 0) {
            array = new Complex(1, 1);
//            for (int i=0; i< nAll; i++){
//                array[i] = new Complex(i, i);
//            }

            PCJ.broadcast(array, Shared.array);
        }

        PCJ.barrier();

        if (PCJ.myId() == 3) System.out.println(array.add(array));

    }

    public static void main(String[] args) throws Throwable {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}