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
    startTime = data["time"]["startTime"]
    endTime = data["time"]["endTime"]
    dt = data["time"]["dt"]

    with h5py.File(h5FilePath, 'r') as h5File:
        calculationFolders = list(h5File.keys())
        calculationFolders.remove("__DATA_TYPES__")

        calculationFolders = sorted(calculationFolders, key=int, reverse=True)
        datasetName = calculationFolders[0]

        timesteps = int((float(endTime) - float(startTime)) / dt) + 1
        y = []
        yDer = []
        for i in range(timesteps):
            y.append(h5File[datasetName + "/y/" + str(i)][:])
            yDer.append(h5File[datasetName + "/yDer/" + str(i)][:])
        x = h5File[datasetName + "/x"][:]
        potential = h5File[datasetName + "/potential"][:]


    for func, funcDer in zip(y, yDer):
        func = np.array([complex(x, y) for x, y in func])
        funcDer = np.array([complex(x, y) for x, y in funcDer])

    for func, funcDer in zip(y, yDer):
        diagrams.ComplexThreePlots(x, func, potential, funcDer, title1="Wave Function", title2="Potential", title3="Wave Function Derivative")
