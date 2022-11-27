package parallelIntegrator;

import org.apache.commons.math3.complex.Complex;

import java.io.Serializable;

public class ParalellComplexIntegratorStorage implements Serializable {
    Complex integralValue;

    public ParalellComplexIntegratorStorage(Complex integralValue) {
        this.integralValue = integralValue;
    }

    public ParalellComplexIntegratorStorage() {
        this.integralValue = Complex.ZERO;
    }

    public ParalellComplexIntegratorStorage add(ParalellComplexIntegratorStorage element){
        return new ParalellComplexIntegratorStorage(
                this.integralValue.add(element.integralValue)
        );
    }

    public Complex getIntegralValue() {
        return integralValue;
    }

    public void setIntegralValue(Complex integralValue) {
        this.integralValue = integralValue;
    }
}