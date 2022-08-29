package compute;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;

@RegisterStorage(RunJob.Shared.class)
public class RunJob implements StartPoint {

    static final int RESOLUTION_EXPOTENTIAL = 6;
    static final int PERIOD_EXPOTENTIAL = 3;

    static double PI = Math.PI;
    static double period = PERIOD_EXPOTENTIAL * Math.PI;

    @Storage(RunJob.class)
    enum Shared {resolution, threadCount, waveFunction, x}

    int resolution = (int) Math.pow(2, RESOLUTION_EXPOTENTIAL);
    int threadCount =  PCJ.threadCount();
    ArrayList<Complex> waveFunction = new ArrayList<>();
    ArrayList<Complex> x = new ArrayList<>();

    @Override
    public void main() {

        if (PCJ.myId() == 0) {

            double step = period / resolution;
            for (double i = -(period / 2); i <= (period / 2); i += step){
                Complex elem = new Complex(i, 0);
                x.add(new Complex(i, 0));
                waveFunction.add(utils.functionsMisc.waveFunction(elem));
            }

            PCJ.broadcast(x, Shared.x);
            PCJ.broadcast(waveFunction, Shared.waveFunction);
        }

        PCJ.barrier();

        if (PCJ.myId() == 3) System.out.println(waveFunction.get(31));

    }

    public static void main(String[] args) throws Throwable {
        String nodesFile = "nodes.txt";
        PCJ.executionBuilder(RunJob.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}