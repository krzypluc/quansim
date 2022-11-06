import h5py
import numpy as np
import yaml
from yaml.loader import SafeLoader

from utils import diagrams

if __name__ == "__main__":
    # Load config
    with open("../config/config.yml") as yamlFile:
        data = yaml.load(yamlFile, Loader=SafeLoader)

    h5FilePath = data["filePaths"]["hdf5PythonFile"]

    with h5py.File(h5FilePath, 'r') as h5File:
        calculationFolders = list(h5File.keys())
        calculationFolders.remove("__DATA_TYPES__")

        calculationFolders = sorted(calculationFolders, key=int, reverse=True)
        datasetName = calculationFolders[0]

        y_derr = h5File[datasetName + "/y_double_derr"][:]
        y = h5File[datasetName + "/y"][:]
        x = h5File[datasetName + "/x"][:]

    y_derr = np.array([complex(x, y) for x, y in y_derr])

    diagrams.ComplexTwoPlots(x, y, y_derr, "Wave function", "Derivative of WF")
