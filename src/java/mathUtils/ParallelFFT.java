package mathUtils;

import compute.RunComputations;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.pcj.RegisterStorage;
import org.pcj.StartPoint;

@RegisterStorage(RunComputations.SharedRunJob.class)
public class ParallelFFT implements StartPoint {
    int procID = PCJ.myId();
    int procCount = PCJ.threadCount();

    Complex[] y;
    double[] x;

    static int runNumber;

    public ParallelFFT(Complex[] y, double[] x) {
        this.y = y;
        this.x = x;

        runNumber = 0;
    }

    public void main () {
        int lengthOfPiece = (int) y.length / procCount;
        Complex[] yFragment = new Complex[lengthOfPiece];

        int j = 0;
        for (int i = 0; i < y.length; i++) {
            // If process is in first half
            if (procID < (procCount / 2)){
                if ((i % procCount) == procID * 2) {
                    yFragment[j] = y[i];
                    j++;
                }
            }

            // If process is in the second half
            else {
                if ((i % procCount) == (procID - procCount / 2) * 2 + 1) {
                    yFragment[j] = y[i];
                    j++;
                }
            }

        }

        // Run number used in paralellization.
        runNumber++;
    }
}
