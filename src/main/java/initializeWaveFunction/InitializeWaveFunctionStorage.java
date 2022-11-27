package initializeWaveFunction;

import org.apache.commons.math3.complex.Complex;

import java.io.Serializable;

public class InitializeWaveFunctionStorage implements Serializable {
    Complex[] y;
    double[] x;
    double[] potential;

    public InitializeWaveFunctionStorage(Complex[] y, double[] x, double[] potential) {
        this.y = y;
        this.x = x;
        this.potential = potential;
    }

    public Complex[] getY() {
        return y;
    }

    public void setY(Complex[] y) {
        this.y = y;
    }

    public double[] getX() {
        return x;
    }

    public void setX(double[] x) {
        this.x = x;
    }

    public double[] getPotential() {
        return potential;
    }

    public void setPotential(double[] potential) {
        this.potential = potential;
    }
}
