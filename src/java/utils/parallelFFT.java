package utils;


import compute.RunJob;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.pcj.PCJ;
import org.pcj.PcjFuture;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static compute.RunJob.SharedRunJob.nextYFragment;

public class parallelFFT {
    public static Complex[] FFTComplex1D(Complex[] y) {
        int procID = PCJ.myId();
        int procCount = PCJ.threadCount();

        int lengthOfpiece = (int) y.length / procCount;

        int numberOfLoops = (int) (Math.log(procCount) / Math.log(2));

        Complex[] yFragment = new Complex[lengthOfpiece];

        int j = 0;
        for (int i = 0; i < y.length; i++) {
            if (i % procCount == procID) {
                yFragment[j] = y[i];
                j++;
            }
        }

        HashSet<Integer> processesNumbers = new HashSet<Integer>();

        for (int i = 0; i < procCount; i++) {
            processesNumbers.add(i);
        }

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] transformed = fft.transform(yFragment, TransformType.FORWARD);

        Complex[] yFragmentBuffer = new Complex[];

        while (processesNumbers.size() > 1) {
            for (int i = 0; i < processesNumbers.size(); i++) {
                if (i % 2 == 0) {
                    processesNumbers.remove(i);
                }

                if (processesNumbers.contains(procID)) {
                    PCJ.waitFor(nextYFragment);
                    yFragmentBuffer = new Complex[yFragment.length * 2];

                    for (int i = 0; i < yFragment.length; i++) {
                        yFragmentBuffer[2 * i] = yFragment[i];
                        yFragmentBuffer[2 * i + 1] = nextYFragment[i];
                    }

                    yFragment = yFragmentBuffer;


                }

                if (!processesNumbers.contains(procID)) {
                    PCJ.asyncPut(transformed, procID - 1, nextYFragment)
                }

            }
        }

        ArrayList<Complex> nextYFragment = future.get();


        return yFragment;
    }
}
