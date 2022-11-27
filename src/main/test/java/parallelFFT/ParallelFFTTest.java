package parallelFFT;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.Precision;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class ParallelFFTTest {

    @Test
    void testForwardTransform() {

        Complex[] y = new Complex[] {
                Complex.valueOf(1, 0),
                Complex.valueOf(1, 0),
                Complex.valueOf(0, 1),
                Complex.valueOf(0, 1),
                Complex.valueOf(1, 0),
                Complex.valueOf(1, 0),
                Complex.valueOf(0, 1),
                Complex.valueOf(0, 1),
        };

        Complex[] yExpected = new Complex[] {
                Complex.valueOf(4.0, 4.0),
                Complex.valueOf(0, 0),
                Complex.valueOf(0, -4),
                Complex.valueOf(0, 0),
                Complex.valueOf(0, 0),
                Complex.valueOf(0, 0),
                Complex.valueOf(4, 0),
                Complex.valueOf(0, 0)

        };

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] y_transformed = fft.transform(y, TransformType.FORWARD);
        double epsilon = 1e-15;
        for (int i = 0; i < y.length; i++){
            assertTrue(Precision.equals(y_transformed[i].getReal(), yExpected[i].getReal(), epsilon));
            assertTrue(Precision.equals(y_transformed[i].getImaginary(), yExpected[i].getImaginary(), epsilon));
        }
    }
}