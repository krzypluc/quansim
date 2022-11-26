package mathUtils;

import compute.RunJob;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.pcj.RegisterStorage;
import org.pcj.StartPoint;

import java.util.HashMap;

@RegisterStorage(RunJob.SharedRunJob.class)
public class InitializeWaveFunction implements StartPoint {

    Complex[] y;
    double[] x;
    double period;
    double dx;

    public InitializeWaveFunction(Complex[] y, double[] x, double period, double dx) {
        this.period = period;
        this.x = x;
        this.y = y;
        this.dx = dx;
    }

    @Override
    public void main() {
        double[] potential = new double[x.length];

        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;

        // All processes calulculates their part of x, y and potential
        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            x[i] = (-period / 2) + dx * i;
            y[i] = Functions.waveFunction(x[i], Functions.WaveFunction.GAUSSIAN);
            potential[i] = Functions.potential(x[i]);

            // Put all variables on processor nr 0
            PCJ.asyncPut(y[i], 0, RunJob.SharedRunJob.y, i);
            PCJ.asyncPut(x[i], 0, RunJob.SharedRunJob.x, i);
            PCJ.asyncPut(potential[i], 0, RunJob.SharedRunJob.potential, i);
        }


        if (procID == 0) {
            y = PCJ.get(0, RunJob.SharedRunJob.y);
            x = PCJ.get(0, RunJob.SharedRunJob.x);
            potential = PCJ.get(0, RunJob.SharedRunJob.potential);

            PCJ.asyncBroadcast(y, RunJob.SharedRunJob.y);
            PCJ.asyncBroadcast(x, RunJob.SharedRunJob.x);
            PCJ.asyncBroadcast(potential, RunJob.SharedRunJob.potential);
        }
    }
}
