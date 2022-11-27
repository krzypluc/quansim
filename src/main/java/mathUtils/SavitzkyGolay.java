package mathUtils;

import org.apache.commons.math3.complex.Complex;

public class SavitzkyGolay {
    public static Complex smoothWindow7(Complex[] function, int index) {
        int n = function.length;

        if (index == 0){
            Complex elem0, elem1, elem2, elem3;
            elem0 = function[0].multiply(7);
            elem1 = function[1].multiply(6);
            elem2 = function[2].multiply(3);
            elem3 = function[3].multiply(-2);

            return elem0.add(elem1.add(elem2).add(elem3)).divide(21);
        }
        if (index == 1){
            Complex elemp1, elem0, elem1, elem2, elem3;
            elemp1 = function[0].multiply(6);
            elem0 = function[1].multiply(7);
            elem1 = function[2].multiply(6);
            elem2 = function[3].multiply(3);
            elem3 = function[3].multiply(-2);

            return elemp1.add(elem0.add(elem1).add(elem2).add(elem3)).divide(21);
        }
        if (index == 2){
            Complex elemp2, elemp1, elem0, elem1, elem2, elem3;
            elemp2 = function[0].multiply(3);
            elemp1 = function[1].multiply(6);
            elem0 = function[2].multiply(7);
            elem1 = function[3].multiply(6);
            elem2 = function[4].multiply(3);
            elem3 = function[5].multiply(-2);

            return elemp1.add(elem0.add(elem1).add(elem2).add(elem3)).divide(21);
        }

        if (index == n - 1){
            Complex elem0, elem1, elem2, elem3;
            elem0 = function[n - 1].multiply(7);
            elem1 = function[n - 2].multiply(6);
            elem2 = function[n - 3].multiply(3);
            elem3 = function[n - 4].multiply(-2);

            return elem0.add(elem1.add(elem2).add(elem3)).divide(21);
        }
        if (index == n - 2){
            Complex elemp1, elem0, elem1, elem2, elem3;
            elemp1 = function[n - 1].multiply(6);
            elem0 = function[n - 2].multiply(7);
            elem1 = function[n - 3].multiply(6);
            elem2 = function[n - 4].multiply(3);
            elem3 = function[n - 5].multiply(-2);

            return elemp1.add(elem0.add(elem1).add(elem2).add(elem3)).divide(21);
        }
        if (index == n - 3){
            Complex elemp2, elemp1, elem0, elem1, elem2, elem3;
            elemp2 = function[n - 1].multiply(3);
            elemp1 = function[n - 2].multiply(6);
            elem0 = function[n - 3].multiply(7);
            elem1 = function[n - 4].multiply(6);
            elem2 = function[n - 5].multiply(3);
            elem3 = function[n - 6].multiply(-2);

            return elemp1.add(elem0.add(elem1).add(elem2).add(elem3)).divide(21);
        }

        Complex elemp3, elemp2, elemp1, elem0, elem1, elem2, elem3;
        elemp3 = function[index - 3].multiply(-2);
        elemp2 = function[index - 2].multiply(3);
        elemp1 = function[index - 1].multiply(6);
        elem0 = function[index].multiply(7);
        elem1 = function[index + 1].multiply(6);
        elem2 = function[index + 2].multiply(3);
        elem3 = function[index + 3].multiply(-2);

        return elemp3.add(elemp2.add(elemp1).add(elem0).add(elem1).add(elem2).add(elem3)).divide(21);
    }
}
