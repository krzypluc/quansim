package initializeWaveFunction;

import runComutations.RunComputations;
import mathUtils.Potential;
import mathUtils.WaveFunctions;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.pcj.RegisterStorage;
import org.pcj.StartPoint;

@RegisterStorage(RunComputations.SharedRunJob.class)
public class InitializeWaveFunction implements StartPoint {

    Complex[] y;
    double[] x;
    double[] potential;
    double period;
    double dx;
    WaveFunctions.WaveFunction waveFunctionType;
    Potential.PotentialType potentialType;

    public InitializeWaveFunction(Complex[] y, double[] x, double[] potential, double period, double dx, WaveFunctions.WaveFunction waveFunctionType, Potential.PotentialType potentialType) {
        this.y = y;
        this.x = x;
        this.potential = potential;
        this.period = period;
        this.dx = dx;
        this.waveFunctionType = waveFunctionType;
        this.potentialType = potentialType;
    }

    @Override
    public void main() {
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;

        // All processes calulculates their part of x, y and potential
        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            x[i] = (-period / 2) + dx * i;
            y[i] = WaveFunctions.waveFunction(x[i], waveFunctionType);
            potential[i] = Potential.potential(x[i], potentialType);

            // Put all variables on processor nr 0
            PCJ.asyncPut(y[i], 0, RunComputations.SharedRunJob.y, i);
            PCJ.asyncPut(x[i], 0, RunComputations.SharedRunJob.x, i);
            PCJ.asyncPut(potential[i], 0, RunComputations.SharedRunJob.potential, i);
        }

        PCJ.barrier();

        if (procID == 0) {
            InitializeWaveFunctionStorage initWFStorage = new InitializeWaveFunctionStorage(
                    // get y
                    PCJ.get(0, RunComputations.SharedRunJob.y),
                    // get x
                    PCJ.get(0, RunComputations.SharedRunJob.x),
                    // get potential
                    PCJ.get(0, RunComputations.SharedRunJob.potential)
            );

            PCJ.asyncBroadcast(initWFStorage, RunComputations.SharedRunJob.initWFStorage);

        }
    }
}
