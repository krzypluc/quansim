package miscellaneous;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import org.apache.commons.math3.complex.Complex;

import java.util.Map;

public class hdf5Handler {
    public static void saveYandDerivative(
            double[] x,
            Complex[][] yHistory,
            Complex[][] yDerHistory,
            int timesteps, Map<String,
            String> filePaths) {
        int yLength = x.length;

        String groupName = GroupName.getGroupName();
        String hdf5FileName = (String) filePaths.get("hdf5JavaFile");

        IHDF5Writer writer = HDF5Factory.open(hdf5FileName);
        writer.writeDoubleArray(groupName + "/x", x);

        double[][][] yTransformedDouble = new double[timesteps][yLength][2];
        double[][][] yDerTransformedDouble = new double[timesteps][yLength][2];
        for (int i = 0; i < timesteps; i++) {
            for (int j = 0; j < yLength; j++) {
                yTransformedDouble[i][j][0] = yHistory[i][j].getReal();
                yTransformedDouble[i][j][1] = yHistory[i][j].getImaginary();

                yDerTransformedDouble[i][j][0] = yDerHistory[i][j].getReal();
                yDerTransformedDouble[i][j][1] = yDerHistory[i][j].getImaginary();
            }
        }

        for (int i = 0; i < timesteps; i++) {
            writer.writeDoubleMatrix(groupName + "/y/" + i, yTransformedDouble[i]);
            writer.writeDoubleMatrix(groupName + "/yDer/" + i, yDerTransformedDouble[i]);
        }

        writer.close();
    }
}
