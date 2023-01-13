# License: GNU General Public License, Version 3, 29 June 2007
# Copyright © 2022 Helmholtz Centre Potsdam GFZ German Research Centre for
# Geosciences, Potsdam, Germany

#!/usr/bin/env python
# coding: utf-8

'''
    See config.py for input parameter documentation
'''


import sys
import os

import numpy as np
from pyevtk.hl import gridToVTK

import gml

# check if config file name is provided
if len(sys.argv) < 2:
    print('\n\tPlease use this command to run GEOMODELATOR:\n' + \
         '\033[1m\tpython [-i] [-u] [path]geomodelator.py [path]<config file name>' + \
         '[2>&1 | tee [path]<log file name>]\033[0m' + \
         '\n\n\t\033[1mExamples:\033[0m' + \
         '\n\tpython -u geomodelator.py config.py 2>&1 | tee config.log' + \
         '\n\tpython -i geomodelator.py example_3d/config.py' + \
         '\n\tpython -i geomodelator.py ~/geomodelator/example_3d/config.py\n'
            )
    sys.exit()

# import geomodelator config.py
cfile = os.path.basename(sys.argv[1]).split('.')[0]
cdir = os.path.dirname(sys.argv[1]) + "/"
sys.path.insert(0, cdir)
cfg = __import__(cfile)
# import config as cfg


global_settings = gml.Settings(
    s3d = cfg.S3D,
    csv_columns = cfg.CSVCOLUMNS,
    model_p1 = cfg.MODELP1,
    model_p2 = cfg.MODELP2,
    partition_width = cfg.PARTITIONWIDTH,
    dx = cfg.dx,
    dy = cfg.dy,
    dz = cfg.dz
)

global_settings.path_input_data = cdir + global_settings.path_input_data
global_settings.path_output_data = cdir + global_settings.path_output_data


# generate model
model = gml.Model(global_settings)

# generate cell center model from model
cx, cy, cz = model.calculate_cell_centers()

cmodel = gml.Model(global_settings, x=cx, y=cy, z=cz)

# generate partitions from point cloud interpolated to cell center model (x,y)
layers = gml.load_structure_files(cmodel)
faults = gml.load_structure_files(cmodel, "fault")
seams = gml.load_structure_files(cmodel, "seam")

# identify partitions an get a 3d data array
partition_model_data = gml.generate_layer_partitions(cmodel, layers)
all_model_data, fault_model_data = gml.generate_special_partitions(
    cmodel, faults, old_grid_data=partition_model_data, name="fault")
all_model_data, seam_model_data = gml.generate_special_partitions(
    cmodel, seams, old_grid_data=all_model_data, name="seam")

cell_ids, i_ids, j_ids, k_ids = gml.generate_model_cell_ids(cmodel)

vtk_grid_x, vtk_grid_y, vtk_grid_z = gml.rotate_grid(model)

# combine partition model data with model as vtk
gridToVTK(
    global_settings.path_output_data + "model",
    vtk_grid_x,
    vtk_grid_y,
    vtk_grid_z,
    cellData = {
        "cell_ids": cell_ids,
        "all": all_model_data,
        "partitions": partition_model_data,
        "faults": fault_model_data,
        "seams": seam_model_data,
        "i": i_ids,
        "j": j_ids,
        "k": k_ids
    }
)

print("")
print("")
print("====================================")
print("\tGenerate Output")
print("====================================")
print("")
print("- VTK file with model including all data:    ", end="")
print(global_settings.path_output_data + "model.vtu")
print("- VTK files for all used layer surfaces:     ", end="")
print(global_settings.path_output_data + "layer-xx.vtu")
print("- VTK files for all used fault surfaces:     ", end="")
print(global_settings.path_output_data + "fault-xx.vtu")
print("- VTK files for all used seam surfaces:      ", end="")
print(global_settings.path_output_data + "seam-xx.vtu")
print("")

print("- Numpy array files with model data")

FILENAME = global_settings.path_output_data + "model_all_data_numpy_array"
np.savez_compressed(FILENAME, array1=all_model_data)
print("\t* all:             " + FILENAME + ".npz")

FILENAME = global_settings.path_output_data + "model_cell_id_numpy_array"
np.savez_compressed(FILENAME, array1=cell_ids)
print("\t* cell ids:        " + FILENAME + ".npz")

FILENAME = global_settings.path_output_data + "layer_partition_data_numpy_array"
np.savez_compressed(FILENAME, array1=partition_model_data)
print("\t* layer partition: " + FILENAME + ".npz")

FILENAME = global_settings.path_output_data + "fault_partition_data_numpy_array"
np.savez_compressed(FILENAME, array1=fault_model_data)
print("\t* fault partition: " + FILENAME + ".npz")

FILENAME = global_settings.path_output_data + "seam_partition_data_numpy_array"
np.savez_compressed(FILENAME, array1=seam_model_data)
print("\t* seam partition:  " + FILENAME + ".npz")

print("")
print("")
print("====================================")
print("\tAll done!")
print("====================================")
print("")
