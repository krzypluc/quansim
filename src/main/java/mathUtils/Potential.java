package mathUtils;

public class Potential {
    public enum PotentialType {
        HARMONIC,
        NONE
    }

    public static double potential(double x, PotentialType potential){

        switch (potential) {
            case HARMONIC:
                return Math.pow(x, 2) / 2;
            case NONE:
                return 0.0;

            default:
                throw new RuntimeException("There is no such implemented potential.");
        }
    }
}
