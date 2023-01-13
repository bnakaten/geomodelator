# GEOMODELATOR

A python library to generate simple structured 2d+ and 3d cell-based VTK models/files with layer, fault and seams as corresponding element/zone groups for further processing or integrating into a simulator framework (TRANSPORTSE).

[For further documentation use Wiki](https://git.gfz-potsdam.de/bnakaten/geomodelator/-/wikis)

![GEO MODEL generaTOR to build up 3d cell-based model](images/geomodelator01.png "GEO MODEL generATOR to build up 3d cell-based model")


# Requirements

**GEOMODELATOR** requires Python >= 3.6 and uses the following packages:
  - ipython
  - matplotlib
  - numpy
  - pyevtk
  - scipy

# I) Preprocessing

1. Download **GEOMODELATOR** and enter the directory.

    ```
    git clone https://git.gfz-potsdam.de/bnakaten/geomodelator

    cd geomodelator
    ```


2. Setup conda envrionment conda. 

    ```
    conda env create -f environment.yml

    conda activate geomodelator
    ```
 

3. Create an new model project directory by copying a suitable example directory.

    For 2d/3d with or without rotated model data:
    ```
    cp -r examples/3d [model name]
    ```
    
    Enter the model directory:
    ```
    cd [model name]
    ```

4. Copy csv point cloud files to the 1_input directory. Rename the files equivalent to the described name scheme.

    1. Point cloud files which represent top or bottom of a geological formation:
        - layer-01.csv
        - layer-02.csv
        - layer-03.csv
        - ...
    
    2. Point cloud files which represent fault structures:
        - fault-01.csv
        - fault-02.csv
        - fault-03.csv
        - ...
    
    3. Point cloud files which represent seam structures:
        - seam-01.csv
        - seam-02.csv
        - seam-03.csv
        - ...

    4. Optional: The discretization of the model can be provide by numpy files (numpy.save(...)) for dx, dy and dz.


# II) Generate Model

1. Setup the model configuration by following the inline documentation.

    ```
    vi config.py
    ```

2. Model generation.

    ```
    cd ..

    ipython geomodelator.py [model name]/config.py
    ```
    
    -> Output:
    
    - VTK file with 3d model including all data: [model name]/2_output/model.vtu
    - VTK files for all used layer surfaces:     [model name]/2_output/layer-xx.vtu
    - VTK files for all used fault surfaces:     [model name]/2_output/fault-xx.vtu
    - VTK files for all used seam surfaces:      [model name]/2_output/seam-xx.vtu
    
    - Numpy array files with 3d model data
        * all:            [model name]/2_output/model_all_data_numpy_array.npz
        * cell ids:       [model name]/2_output/model_cell_id_numpy_array.npz
        * layer parition: [model name]/2_output/model_layer_parition_data_numpy_array.npz
        * fault parition: [model name]/2_output/model_fault_parition_data_numpy_array.npz
        * seam parition:  [model name]/2_output/model_seam_parition_data_numpy_array.npz
    
    Note: 
    Example how to use the compressed numpy data arrays:
    
    ```
    lgrid = numpy.load('model_all_data_numpy_array.npz')
    model_data = lgrid['array1']
    ```
