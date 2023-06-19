# License: GNU General Public License, Version 3, 29 June 2007
# Copyright © 2022 Helmholtz Centre Potsdam GFZ German Research Centre for
# Geosciences, Potsdam, Germany

#!/usr/bin/env python
# coding: utf-8

'''
    GEOMODELATOR configuration file
'''
import os
import numpy as np

## enter the path to the project directory relative to the 
## calling python script e.g. geomodelator.py
CALLDIR = 'examples/3d/' 

## type of csv data file 
## in case of 2d data: False
## in case of 3d data: True
S3D = True

# number and order of coordinate columns in csv data files e.g.
## ['x', 'z'] := two columns
## ['y', 'x'] := two columns
## ['y', 'z', 'x'] := three columns
CSVCOLUMNS = ['x', 'y', 'z']
## delimiter in csv file
CSVDELIMITER = ','
## header in csv file (1) or not (0)
CSVHEADER = 1

#    ====================================
#        Model
#    ====================================
#
#                  +-------------------------------------+
#                 /|                                    /|
#                / |                                   / |
#               /  |                                  /  |
#              /   |                                 /   |
#             /    |                                /    |
#            /     |                               /     |
#           +-------------------------------------+      |
#           |      |                              |      |
#           |      |                              |      |
#           |      |                              |      |
#           |      +------------------------------|------+
#           |     /                               |     /
#           |    /                                |    /
#           |   /                                 |   /
#           |  /                                  |  /
#           | /                                   | /
#           |/                                    |/
#           +-------------------------------------+
#       MODELP1                               MODELP2
#
## left lower front corner of the model
MODELP1 = [3000, 6500, -2000]

## right lower front corner of the model (MODELP2 and MODELP1 will be
## both used in case of a rotated model in contrary to a model with 
## boundaries parallel to the x, y and z Cartesian axes.
MODELP2 = [10000, 6500, -2000]

## x,y and z distance of the model
MODELDIMENSION = [3900, 11900, 2000]

################################################################################
##                  Skip the following lines until line 107 and do not change!!!
################################################################################

# PATH AND FILE DETAILS ########################################################

## path details

### input path should contain for each layer file (.csv tab-delimiter and 1 line
### header, .txt or .xyz delimiter ' ' without header) and min. 3 lines
### of tab-delimited coordinates (x y z).
IPATH = "1_input/"

### output directory will be created if not exists. The content files will
### be overwritten if necessary.
OPATH = "2_output/"

################################################################################
##                  Modify the following lines if needed !!!
################################################################################

# MODEL DIMENSION AND SPATIAL DISCRETIZATION ##################################

## Example 1:
## Define by numpy arrays

# dx = np.array([0.5, 0.5, 0.5,  0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5])
# dy = np.array([0.5, 0.5, 0.5, 0.5, 0.3, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5])
# dz = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5])

## Example 2:
## Define consistent discretization by using nx, ny and nz

nx = 40
ny = 80
nz = 40
dx = np.asarray([MODELDIMENSION[0]/nx]*nx)

# dy = np.asarray([MODELDIMENSION[1]/ny]*ny)
r = int(MODELDIMENSION[1]/ny)
dy = np.ones(5)*(r*4)
dy = np.append(dy, np.ones(4)*(r*2))
dy = np.append(dy, np.ones(34)*r)
dy = np.append(dy, np.ones(4)*(r*2))
dy = np.append(dy, np.ones(2)*(r*4))

# dz = np.asarray([MODELDIMENSION[2]/nz]*nz)
r = int(MODELDIMENSION[2]/nz)
dz = np.ones(6)*(r*2)
dz = np.append(dz, np.ones(20)*r)
dz = np.append(dz, np.ones(4)*(r*2))

## Example 3:
## Define by numpy array files

# dx = np.load(IPATH + 'dx.npy')
# dy = np.load(IPATH + 'dy.npy')
# dz = np.load(IPATH + 'dz.npy')

## Example 4:
## Define by numpy arrays with inconsistent distances (tartan grid)

## starting with 30x 10 m element discretization in x direction
# dx = np.ones(30)*10
## followed by 160x 5 m element discretization in x direction
# dx = np.append(dx, np.ones(160)*5)
## ending with 90x 10 m element discretization in x direction
# dx = np.append(dx, np.ones(90)*10)

## in y direction only 1 element
# dy = np.asarray([10]*1)

## starting with 15x 10 m element discretization in z direction
# dz = np.ones(15)*10
## followed by 60x 5 m element discretization in z direction
# dz = np.append(dz, np.ones(60)*5)
## ending with 90x 10 m element discretization in z direction
# dz = np.append(dz, np.ones(90)*10)

# dx = np.load(CALLDIR + IPATH + 'dx.npy')
# dy = np.load(CALLDIR + IPATH + 'dy.npy')
# dz = np.load(CALLDIR + IPATH + 'dz.npy')


# MODEL FAULT AND SEAM PARAMETER ##################################

PARTITIONWIDTH = {}

## Extend this lines for defined fault or seam width
## by using the filename without file extension as specifier.
## If not specified the width will be one element on each
## side of the fault or seam surface.
## The width should be a positive value and is the distance on each side of the
## layer. Actually the width value is the have of the real fault or seam width.
PARTITIONWIDTH["fault-01"] = 10
PARTITIONWIDTH["fault-09"] = 10
PARTITIONWIDTH["fault-12"] = 10

# MODEL LAYER state  ##################################
LAYERSTATE = {}

# LAYERSTATE["layer-01"] = 0
