import sys
import os

import gmlConfiguration as gmlCM
import gmlRun as gmlRM
from gml import Helper as gmlH


def main():
    if len(sys.argv) != 2:
        print("Usage: python geomodelator.py [<PATH>/]<yaml configuration file>\n" \
            "  - <PATH> is optional; if not provided, the script will search in " \
            "the current directory.")
        sys.exit(1)

    configurationFile = sys.argv[1]
    gmlConfigurationJson = False

    if not os.path.exists(configurationFile):
        print(f"Error: File '{configurationFile}' not found.")
        sys.exit(1)

    if configurationFile.split('.')[-1] == 'yaml' or \
        configurationFile.split('.')[-1] == 'yml':
        gmlConfigurationJson = gmlH.yamlToJson(configurationFile)
    elif configurationFile.split('.')[-1] == 'json':
        gmlConfigurationJson = gmlH.loadJson(configurationFile)
    elif configurationFile.split('.')[-1] == 'py':
        gmlConfigurationJson = gmlH.pyConfToJson(configurationFile)

        b = 3

    if gmlConfigurationJson:
        gmlConfiguration = gmlCM.gmlConfiguration(gmlConfigurationJson)
        gmlConfiguration.configureModel()

        gmlRun = gmlRM.gmlRun(gmlConfiguration)
        gmlRun.generateModel()

        a = 2


if __name__ == '__main__':
    main()
