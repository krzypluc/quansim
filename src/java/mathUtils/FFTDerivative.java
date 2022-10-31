package mathUtils;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class FFTDerivative {

    public static Complex[] derivativeComplex(Complex[] y, double[] x, int diffDegree) {
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] y_transformed = fft.transform(y, TransformType.FORWARD);
        double[] fk = FTUtils.freq(x);

        for (int i = 0; i < y_transformed.length; i++) {
            y_transformed[i] = y_transformed[i]
                    .multiply(Complex.I.multiply(fk[i] * 2 * Math.PI).pow(diffDegree));

        }

        Complex[] y_derr = fft.transform(y_transformed, TransformType.INVERSE);
        return y_derr;
    }

    public static Complex[] derivativeComplex(Complex[] y, double[] x) {
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] y_transformed = fft.transform(y, TransformType.FORWARD);
        double[] fk = FTUtils.freq(x);

        for (int i = 0; i < y_transformed.length; i++) {
            y_transformed[i] = y_transformed[i]
                    .multiply(Complex.I.multiply(fk[i] * 2 * Math.PI).pow(2));

            if (y_transformed[i].isNaN()) {
                y_transformed[i] = Complex.ZERO;
            }
        }

        Complex[] y_derr = new Complex[y.length];
        y_derr = fft.transform(y_transformed, TransformType.INVERSE);

        for (Complex nr : y_derr) {
            System.out.println(nr);
        }

        return y_derr;
    }
}
