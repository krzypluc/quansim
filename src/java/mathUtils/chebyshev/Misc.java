package mathUtils.chebyshev;

import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Misc {
    public static double getR(double dt, double[] potential, double mass){
        Double[] pot = ArrayUtils.toObject(potential);
        List list = Arrays.asList(pot);
        int min = (int) Collections.min(list);
        int max = (int) Collections.max(list);

        return min;
    }
}
