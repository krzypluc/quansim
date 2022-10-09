package utils;

import compute.Complex;

import java.util.function.Function;


public class mapDistribute {
    public static Complex[] mapComplex1D(Complex[] array, Function<Complex, Complex> func, int processId, int numberOfProcesses) {
        Complex[] ret = new Complex[array.length];

        int lengthOfpiece = (int) array.length / numberOfProcesses;

        for (int i = processId * lengthOfpiece; i < (processId + 1) * lengthOfpiece; i++){
            ret[i] = func.apply(array[i]);
        }

        return ret;
    }
}
