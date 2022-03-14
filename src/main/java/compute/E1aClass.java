package compute;

import java.util.Random;
import java.io.File;
import org.pcj.PCJ;
import org.pcj.StartPoint;
import org.pcj.Storage;
import org.pcj.RegisterStorage;

/*
author: Krzysztof Płuciennik
*/

@RegisterStorage(E1aClass.Shared.class)
public class E1aClass implements StartPoint {

    @Storage(E1aClass.class)
    enum Shared { nAll, A, B, C, startTime}
    int nAll = 1048576;
    double[] A = new double[nAll];
    double[] B = new double[nAll];
    double[] C = new double[nAll];
    long startTime;

    @Override
    public void main() {
        PCJ.barrier();

        //Obliczenie długości gragmentów dla poszczególnych procesorów.
        int n = nAll / PCJ.threadCount();

        if (PCJ.myId() == 0) {
            //Rozpoczecie odliczania czasu pracy
            startTime = System.nanoTime();

            //Wypełnianie liczbami losowymi
            Random r = new Random();
            for (int i = 0; i < nAll; i++){
                B[i] = r.nextDouble();
                C[i] = r.nextDouble();
            }
            //Broadcast wylosowanych liczb do procesorów
            PCJ.broadcast(B, Shared.B);
            PCJ.broadcast(C, Shared.C);
        }

        PCJ.barrier();

        //Każdy procesor wykonuje obliczenia dla przydzielonego fragmentu tablicy
        for (int i = PCJ.myId() * n ; i < n * (PCJ.myId() + 1); i++){
            A[i] = B[i] + C[i];
        }

        PCJ.barrier();

        //Procesor o numerze 0 zbiera wszystkie wyniki.
        if (PCJ.myId() == 0) {
            for (int procNum = 1; procNum < PCJ.threadCount(); procNum++){
                for (int i = procNum * n; i < (procNum + 1) * n; i++){
                    A[i] = PCJ.get(procNum, Shared.A, i);
                }
            }

            //Procesor wypisuje na ekran 1 i ostatni element wektora.
            System.out.println(A[0]);
            System.out.println(A[nAll-1]);

            //Obliczenie i wypisanie czasu pracy.
            long elapsedTime = System.nanoTime() - startTime;
            System.out.println("Liczba procesorow: " + PCJ.threadCount() + " Czas: " + elapsedTime * 1.0E-9);
        }
    }

    public static void main(String[] args) throws Throwable {
        String nodesFile  = "nodes.txt";
        PCJ.executionBuilder (E1aClass.class)
                .addNodes(new File(nodesFile))
                .start();
    }
}