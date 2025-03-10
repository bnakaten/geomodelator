
[![DOI:10.5880/GFZ.3.4.2025.001](https://img.shields.io/badge/DOI-10.5880%2FGFZ.3.4.2025.001-blue.svg)](https://doi.org/10.5880/GFZ.3.4.2025.001) [![DOI:10.5880/GFZ.3.4.2024.003](https://img.shields.io/badge/DOI-10.5880%2FGFZ.3.4.2024.003-blue.svg)](https://doi.org/10.5880/GFZ.3.4.2024.003)

# GEOMODELATOR

A python library to generate simple structured 2.5D and 3D cell-based VTK models/files with layer, fault and seams as corresponding element/zone groups for further processing or integrating into a simulator framework (TRANSPORTSE).

[Visit Wiki](https://git.gfz.de/bnakaten/geomodelator/-/wikis)

![GEO MODEL generaTOR to build up 3D cell-based model](geomodelator.png "GEO MODEL generATOR to build up 3d cell-based model")


# Requirements

**GEOMODELATOR** requires Python = 3.13 and uses the following packages:
  - python=3.13
  - gdal
  - geos
  - matplotlib
  - numpy
  - pyvista
  - scipy
  - yaml

# I) Installation

1. Download **GEOMODELATOR**.

    ```
    git clone https://git.gfz.de/bnakaten/geomodelator

    #git -c http.sslVerify=false  clone https://git.gfz.de/bnakaten/geomodelator

    cd geomodelator

    git checkout tags/v2.0
    ```


2. Setup conda envrionment conda.

    ```
    mamba env create -f environment.yml

    mamba activate gml
    ```

3. Create a new model project directory or duplicate the existing demo directory and rename it.

4. Copy the point cloud files in CSV format (x, y, z) for the desired horizon and fault surfaces into the input directory.


# II) Use GEOMODELATOR

1. Set up the model details by customizing the configuration file in accordance with the manual.

    ```
    vi configuration.yaml
    ```

2. Generate the model.

    ```
    python geomodelator.py [path to the new project directory]/configuration.yaml
    ```

    **e.g. run demo:**
    ```
    python geomodelator.py examples/demo/configuration.yaml
    ```

    -> Output:

    - VTK file with 3D model including all data: model.vts
    - VTK files for horizon and fault surfaces: ....vtu
    - Numpy array files with all 3d model data: model_all_data_numpy_array.npz

    - Compressed file with all input files: Input.zip
    - Compressed file with all output files: Output.zip

3. The further use of the model data..

    - Use in TRANSPORTSE.
    - Parametrise the model and model cells by custom python code.
    - Visualize the result files in Paraview. For reference, the example demo output_default includes all the output files that are generated when the demo model is executed.

    ```
    allModelData = numpy.load('model_all_data_numpy_array.npz')
    modelData = allModelData['array1']
    ```
