package compute;

import java.io.Serializable;

public class Complex implements Serializable {
    double real;
    double imag;

    public Complex(double real, double imag) {
        this.real = real;
        this.imag = imag;
    }

    public double getReal() {
        return real;
    }

    public double getImag() {
        return imag;
    }

    public Complex add(Complex b){
        double real = this.getReal() + b.getReal();
        double imag = this.getImag() + b.getImag();

        return new Complex(real, imag);
    }

    public Complex subtract(Complex b){
        double real = this.getReal() - b.getReal();
        double imag = this.getImag() - b.getImag();

        return new Complex(real, imag);
    }

    public Complex multiply(Complex b){
        double real = this.getReal() * b.getReal() - this.getImag() * b.getImag();
        double imag = this.getReal() * b.getImag() - this.getImag() * b.getReal();

        return new Complex(real, imag);
    }

    public Complex multiply(double b){
        return new Complex(this.getReal() * b, this.getImag());
    }

    public Complex multiply(int b){
        return new Complex(this.getReal() * b, this.getImag());
    }

    public Complex conjugation(){
        return new Complex(this.getImag(), this.getReal());
    }

    public double modulus(){
        return Math.sqrt(this.getReal() * this.getReal() + this.getImag() * this.getImag());
    }

    @Override
    public String toString() {
        return real + " + " + imag + 'i';
    }
}
