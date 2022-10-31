import h5py
import yaml
from yaml.loader import SafeLoader

if __name__=="__main__":

    # Load config
    with open("../config/config.yml") as yamlFile:
        data = yaml.load(yamlFile, Loader=SafeLoader)

    h5FilePath = data.get("hdf5PythonFolderFilePath")

    with h5py.File(h5FilePath, 'r') as h5File:
        calculationFolders = list(h5File.keys())
        calculationFolders.remove("__DATA_TYPES__")

        calculationFolders = sorted(calculationFolders, key=int, reverse=True)
        datasetName = calculationFolders[0]

        dataset = h5File[datasetName + "/y_double_derr"][:]

    print(dataset)