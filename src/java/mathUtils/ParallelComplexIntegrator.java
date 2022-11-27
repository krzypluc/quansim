package mathUtils;

import compute.RunComputations;
import org.apache.commons.math3.complex.Complex;
import org.pcj.PCJ;
import org.pcj.RegisterStorage;
import org.pcj.StartPoint;
import org.pcj.Storage;


@RegisterStorage(RunComputations.SharedRunJob.class)
public class ParallelComplexIntegrator implements StartPoint {
    Complex[] y;
    double[] x;
    public ParallelComplexIntegrator(Complex[] y, double[] x) {
        this.y = y;
        this.x = x;
    }

    @Override
    public void main() {
        double dx = x[1] - x[0];
        Complex sumOfValues = Complex.ZERO;
        Complex firstValue = Complex.ZERO;
        Complex lastValue = Complex.ZERO;

        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();
        int lengthOfpiece = (int) y.length / procCount;

        for (int i = procID * lengthOfpiece; i < (procID + 1) * lengthOfpiece; i++) {
            sumOfValues = sumOfValues.add(y[i]);

            // Saving first value - use in integral
            if (i == procID * lengthOfpiece) {
                firstValue = y[i];
            }

            // Saving last value - use in integral
            if (i == (procID + 1) * lengthOfpiece - 1) {
                lastValue = y[i];
            }
        }

        Complex parIntegral = sumOfValues
                .multiply(2.0)
                .subtract(lastValue)
                .subtract(firstValue)
                .multiply(dx / 2);

        PCJ.put(parIntegral, PCJ.myId(), RunComputations.SharedRunJob.integral);

        if (procID == 0) {
            parIntegral = PCJ.reduce(
                    (subtotal, element) -> subtotal.add(element),
                    RunComputations.SharedRunJob.integral
            );
            // Broadcast to all processes
            PCJ.broadcast(parIntegral, RunComputations.SharedRunJob.integral);
        }
    }
}
