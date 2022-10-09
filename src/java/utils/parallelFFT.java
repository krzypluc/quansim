package utils;


import compute.RunJob;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.PCJ;
import org.pcj.PcjFuture;

import java.util.ArrayList;
import java.util.List;

public class parallelFFT {
    public static Complex[] FFTComplex1D(Complex[] y){
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        int lengthOfpiece = (int) y.length / procCount;

        int numberOfLoops = (int) (Math.log(procCount) / Math.log(2));

        Complex[] yFragment = new Complex[lengthOfpiece];

        int j = 0;
        for (int i = 0; i < y.length; i++){
            if (i % procCount == procID){
                yFragment[j] = y[i];
                j++;
            }
        }

        PcjFuture<ArrayList<Complex>> future = null;
        future = PCJ.asyncGet(procID + 1, RunJob.SharedRunJob.nextYFragment);

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] transformed = fft.transform(yFragment, TransformType.FORWARD);

        ArrayList<Complex> nextYFragment = future.get();

        return yFragment;
    }
}
