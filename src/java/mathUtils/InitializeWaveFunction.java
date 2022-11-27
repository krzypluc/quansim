package mathUtils;

import compute.RunComputations;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.pcj.RegisterStorage;
import org.pcj.StartPoint;

@RegisterStorage(RunComputations.SharedRunJob.class)
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
            PCJ.asyncPut(y[i], 0, RunComputations.SharedRunJob.y, i);
            PCJ.asyncPut(x[i], 0, RunComputations.SharedRunJob.x, i);
            PCJ.asyncPut(potential[i], 0, RunComputations.SharedRunJob.potential, i);
        }


        if (procID == 0) {
            y = PCJ.get(0, RunComputations.SharedRunJob.y);
            x = PCJ.get(0, RunComputations.SharedRunJob.x);
            potential = PCJ.get(0, RunComputations.SharedRunJob.potential);

            PCJ.asyncBroadcast(y, RunComputations.SharedRunJob.y);
            PCJ.asyncBroadcast(x, RunComputations.SharedRunJob.x);
            PCJ.asyncBroadcast(potential, RunComputations.SharedRunJob.potential);
        }
    }
}
