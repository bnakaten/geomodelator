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
import glob

import numpy as np
from pyevtk.hl import gridToVTK

import gml

# check if config file name is provided
if len(sys.argv) < 2:
    print('\033[1mError: Missing config file!\033[0m')
    gml.print_help()

# import geomodelator config.py
cfile = os.path.basename(sys.argv[1]).split('.')[0]
cdir = os.path.dirname(sys.argv[1])
if os.path.dirname(sys.argv[1]) != "":
    cdir += "/"

sys.path.insert(0, cdir)

filenamelist = sorted(glob.glob(cdir + cfile + '.py'))
if not filenamelist:
    print('\033[1mError: Config file does not exist!\033[0m')
    gml.print_help()
    sys.exit(1)

cfg = __import__(cfile)

global_settings = gml.Settings(
    s3d = cfg.S3D,
    csv_columns = cfg.CSVCOLUMNS,
    model_p1 = cfg.MODELP1,
    model_p2 = cfg.MODELP2,
    model_dimension = cfg.MODELDIMENSION,
    # model_p8 = cfg.MODELP2,
    partition_width = cfg.PARTITIONWIDTH,
    dx = cfg.dx,
    dy = cfg.dy,
    dz = cfg.dz
)

global_settings.path_input_data = cdir + global_settings.path_input_data
global_settings.path_output_data = cdir + global_settings.path_output_data
global_settings.path_shape_data = cdir + global_settings.path_shape_data


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
partition_model_data = gml.generate_layer_partitions(
    cmodel, layers
)
all_model_data, fault_model_data=gml.generate_special_partitions(
    cmodel, faults, old_grid_data=partition_model_data, name="fault"
)
all_model_data, seam_model_data=gml.generate_special_partitions(
    cmodel, seams, old_grid_data=all_model_data, name="seam"
)

# mask the model into active and inactive zones by using a shapefile
all_model_data, active_model_data=gml.generate_active_partitions(
    cmodel, old_grid_data=all_model_data, name="active"
)

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
        "active": active_model_data,
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
print("- VTK surface files for all used structures:     ", end="")
print(global_settings.path_output_data + "....vtk")
print("- VTK point files for all used sturctures:      ", end="")
print(global_settings.path_output_data + "....vtk")
print("")

print("- Numpy array files with model data")

FILENAME = global_settings.path_output_data + "model_all_data_numpy_array"
np.savez_compressed(FILENAME, array1=all_model_data)
print("\t* all:             " + FILENAME + ".npz")

FILENAME = global_settings.path_output_data + "model_detail_data_numpy_array"
np.savez_compressed(FILENAME, cellids=cell_ids, layers=partition_model_data, \
faults=fault_model_data, seams=seam_model_data, active=active_model_data)
print("\t* details:        " + FILENAME + ".npz")

print("")
print("")
print("====================================")
print("\tAll done!")
print("====================================")
print("")
