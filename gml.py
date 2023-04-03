# License: GNU General Public License, Version 3, 29 June 2007
# Copyright © 2022 Helmholtz Centre Potsdam GFZ German Research Centre for
# Geosciences, Potsdam, Germany

#!/usr/bin/env python
# coding: utf-8

'''
GEOMODELATOR library file
'''
import os
import sys
import glob
import math

import numpy as np

from scipy.interpolate import griddata
from scipy.spatial.qhull import Delaunay

from osgeo import ogr, gdal

from pyevtk.hl import unstructuredGridToVTK, pointsToVTK

import pyvista as pv

# CLASSES AND FUNCTIONS ########################################################
class Settings():
    '''
        Data structure for global settings
    '''
    # global static variable
    s3d = True
    path_input_data = ''
    path_output_data = ''
    csv_columns = None
    csv_header = 1
    csv_delimiter = ','
    model_p1 = None
    model_p2 = None
    model_b1 = None
    model_b2 = None
    model_p8 = None
    model_dimension = None
    partition_width = None
    dx = None
    dy = None
    dz = None

    # global cache variable
    rotation_angle = 0
    partition_number = 1

    def __init__(
        self,
        s3d,
        csv_columns,
        model_p1,
        model_p2,
        model_dimension,
        partition_width,
        dx,
        dy,
        dz,
        path_input_data="1_input/",
        path_output_data="2_output/",
        path_shape_data="1_input/shapefiles/",
        mask_file='model_mask_raster.tif'
    ):
        self.s3d = s3d
        self.csv_columns = csv_columns
        self.model_p1 = model_p1
        self.model_p2 = model_p2
        self.model_dimension = model_dimension
        self.partition_width = partition_width
        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.path_input_data = path_input_data
        self.path_shape_data = path_shape_data
        self.path_output_data = path_output_data

        self.mask_file = mask_file


class Model():
    '''
        Model holds all important data
    '''
    settings = None
    x = None
    y = None
    z = None

    def __init__(self, settings, x=None, y=None, z=None):
        '''
            Constructor for Model class.
        '''
        self.settings = settings

        if x is None:
            x = np.zeros(len(self.settings.dx)+1)
            y = np.zeros(len(self.settings.dy)+1)
            z = np.zeros(len(self.settings.dz)+1)

            x[0] = 0
            y[0] = 0
            z[0] = 0

            for i, value in enumerate(self.settings.dx):
                x[i+1] = value + x[i]

            for i, value in enumerate(self.settings.dy):
                y[i+1] = value + y[i]

            for i, value in enumerate(self.settings.dz):
                z[i+1] = value + z[i]

            self.x = x
            self.y = y
            self.z = z

            print("\n====================================")
            print("\tModel")
            print("====================================")
            print("                          nx                   ")
            print("                  P8                P7         ")
            print("                 +-----------------+           ")
            print("           ny   /|                /|           ")
            print("               / |               / |           ")
            print("           P5 /  |           P6 /  |   nz      ")
            print("             +-----------------+   |           ")
            print("             |   |             |   |           ")
            print("             |   |             |   |           ")
            print("             |   |P4           |   |P3         ")
            print("             |   +-------------|---+           ")
            print("             |  /              |  /            ")
            print("             | /               | /             ")
            print("             |/                |/              ")
            print("             +-----------------+               ")
            print("           P1                 P2               ")
            print("                                               \n")
            print("Origin coordinate (= P1 original) (x, y, z): (", end="")
            print(f"{self.settings.model_p1[0]},", end="")
            print(f"{self.settings.model_p1[1]},", end="")
            print(f"{self.settings.model_p1[2]})")
            print("Discretization (nx, ny, nz): (", end="")
            print(f"{len(self.settings.dx)}, ", end="")
            print(f"{len(self.settings.dy)}, ", end="")
            print(f"{len(self.settings.dz)})\n")
            print("!Normalized model corner coordinates (x, y, z):")
            print(f"\tP1 ({self.x[0]}, {self.y[0]}, {self.z[0]})")
            print(f"\tP2 ({self.x[-1]}, {self.y[0]}, {self.z[0]})")
            print(f"\tP3 ({self.x[-1]}, {self.y[-1]}, {self.z[0]})")
            print(f"\tP4 ({self.x[0]}, {self.y[-1]}, {self.z[0]})")
            print(f"\tP5 ({self.x[0]}, {self.y[0]}, {self.z[-1]})")
            print(f"\tP6 ({self.x[-1]}, {self.y[0]}, {self.z[-1]})")
            print(f"\tP7 ({self.x[-1]}, {self.y[-1]}, {self.z[-1]})")
            print(f"\tP8 ({self.x[0]}, {self.y[-1]}, {self.z[-1]})\n")
            self.settings.model_p8 = [self.x[0], self.y[-1], self.z[-1]]
            self.settings.model_b1 = [self.x[0], self.y[0], self.z[0]]
            self.settings.model_b2 = [self.x[-1], self.y[-1], self.z[-1]]

            number_of_cells = len(self.settings.dx)*len(self.settings.dy)*\
                len(self.settings.dz)
            print(f"\n\n\tNumber off cells/elements: {number_of_cells}\n")
        else:
            if x is not None:
                self.x = x
                self.nx = len(x)

            if y is not None:
                self.y = y
                self.ny = len(y)

            if z is not None:
                self.z = z
                self.nz = len(z)

                if self.settings.rotation_angle == 0:
                    self.settings.rotation_angle = calculate_full_angle(
                        self.settings.model_p2,
                        self.settings.model_p1
                    )

                    print("Rotation angle: " + \
                        f"{np.degrees(self.settings.rotation_angle).round()}\n")


    def calculate_cell_centers(self):
        '''
            Method to calculate cell/element center coordinates.
        '''
        return (self.x[1:] + self.x[:-1])/2, (self.y[1:] + self.y[:-1])/2, \
                (self.z[1:] + self.z[:-1])/2

class Layer_Set:
    '''
        Layer_Set class
    '''
    def __init__(self):
        self._A = None
        self._B = None
        self._orientation = None

    @property
    def A(self):
        ''' getter surface A '''
        return self._A

    @A.setter
    def A(self, value):
        ''' setter surface A '''
        self._A = value

    @property
    def B(self):
        ''' getter surface B '''
        return self._B

    @B.setter
    def B(self, value):
        ''' setter surface B '''
        self._B = value

    @property
    def orientation(self):
        ''' getter surface orientation '''
        return self._orientation

    @orientation.setter
    def orientation(self, value):
        ''' setter surface orientation '''
        self._orientation = value


def interpolate(data, x_vector, y_vector, indexing='xy', method=0):
    '''
        Method to interpolate a surface to plane (e.g. xy or yz or xz).
    '''
    interp_method = ['nearest', 'linear', 'cubic']

    # coords = np.array([[data[i,0], data[i,1]] for i in range(len(data))])
    # values = np.array(data[:,2])
    coords = (data[:,0], data[:,1])
    values = data[:,2]

    x_grid, y_grid = np.meshgrid(x_vector, y_vector, indexing=indexing)

    try:
        grid = griddata(
            coords,
            values,
            (x_grid, y_grid),
            method=interp_method[method]
        )
    except ValueError:
        grid = None

    return grid


def load_and_organize_point_cloud_data(cmodel, filename):
    '''
        Load point cloud from text file
    '''
    data = None

    extensions = ['.csv', '.txt', '.tmp', '.xyz']

    if filename.lower().endswith(tuple(extensions)):
        columns = []

        if 'x' in cmodel.settings.csv_columns:
            x_id = list(cmodel.settings.csv_columns).index('x')
            columns.append(x_id)

        if 'y' in cmodel.settings.csv_columns:
            y_id = list(cmodel.settings.csv_columns).index('y')
            columns.append(y_id)

        if 'z' in cmodel.settings.csv_columns:
            z_id = list(cmodel.settings.csv_columns).index('z')
            columns.append(z_id)

        try:
            print(f"-> {filename}")
            data = np.genfromtxt(
                filename,
                skip_header = cmodel.settings.csv_header,
                delimiter = cmodel.settings.csv_delimiter,
                usecols = tuple(columns)
            )

        except ValueError:
            print('Please check if the csv file column number and the csv \
                column oder array fit together.')
            sys.exit()


    return data

def load_point_cloud_files(filename):
    '''
        Load point cloud from text file
    '''
    data = None
    if filename.lower().endswith('.tmp'):
        data = np.genfromtxt(filename, skip_header=1,  delimiter=',')

    return data



def calculate_full_angle(P, origin):
    '''
        Calculate full angle (360 degrees) between two points
    '''
    return math.radians(
        math.atan2(P[1] - origin[1],
        P[0] - origin[0]) * (180.0 / math.pi)
    )


def rotate_coordinate(cmodel, P, clockwise=-1):
    '''
        Rotate all points in list by rotation angle in origin
    '''
    # rotate x and y coordinates clockwise
    x = (P[0] * math.cos(clockwise*cmodel.settings.rotation_angle)) +\
            (P[1] * math.sin(clockwise*cmodel.settings.rotation_angle))
    y = -(P[0] * math.sin(clockwise*cmodel.settings.rotation_angle)) +\
            (P[1] * math.cos(clockwise*cmodel.settings.rotation_angle))
    P = [x, y, P[2]]

    return P

def shift_coordinate(P, origin):
    '''
        Shift list of points into origin
    '''
    return [P[i] - origin[i] for i in range(3)]


def normalized_grid(model):
    '''
        Rotate model by rotation angle
    '''
    z_matrix, y_matrix, x_matrix = np.meshgrid(
        model.z,
        model.y,
        model.x,
        indexing='ij'
    )

    z_list = [item for layer in z_matrix for row in layer for item in row]
    y_list = [item for layer in y_matrix for row in layer for item in row]
    x_list = [item for layer in x_matrix for row in layer for item in row]

    coord_list = np.array(
        [
            [x_list[i],y_list[i],z_list[i]] for i in range(len(x_list))
        ]
    )

    vtk_grid_x = np.array([coord[0] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))
    vtk_grid_y = np.array([coord[1] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))
    vtk_grid_z = np.array([coord[2] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))

    return vtk_grid_x, vtk_grid_y, vtk_grid_z


def rotate_grid(model):
    '''
        Rotate model by rotation angle
    '''
    z_matrix, y_matrix, x_matrix = np.meshgrid(
        model.z,
        model.y,
        model.x,
        indexing='ij'
    )

    z_list = [item for layer in z_matrix for row in layer for item in row]
    y_list = [item for layer in y_matrix for row in layer for item in row]
    x_list = [item for layer in x_matrix for row in layer for item in row]

    coord_list = np.array(
        [
            [x_list[i],y_list[i],z_list[i]] for i in range(len(x_list))
        ]
    )

    if model.settings.rotation_angle != 0:
        coord_list = np.array(
            [
                rotate_coordinate(model, coord_list[i,:], clockwise=-1)\
                                for i in range(len(coord_list))
            ]
        )

    if model.settings.model_p1 != [0,0,0]:
        coord_list =  np.array(
            [
                shift_coordinate(coord_list[i,:],\
                    [-value for value in model.settings.model_p1])\
                    for i in range(len(coord_list))
            ]
        )

    vtk_grid_x = np.array([coord[0] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))
    vtk_grid_y = np.array([coord[1] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))
    vtk_grid_z = np.array([coord[2] for coord in coord_list]).reshape((\
                            len(model.z),len(model.y),len(model.x)))

    return vtk_grid_x, vtk_grid_y, vtk_grid_z


def convert_2d_data_to_3d_data(cmodel, points):
    '''
        Convert 2d csv data to 2d+/3d data
    '''
    ## convert 2d data to interpolatable 3d data
    # extend data by 3rd dimension coordinate shifted by -50 (default)
    # model_extension_width = math.dist(
        # cmodel.settings.model_p2,
        # cmodel.settings.model_p1
    # )/2. + 50.
    model_extension_width = cmodel.settings.model_dimension[1]/2+50
    print(model_extension_width)


    third_column = np.full(len(points[:,0]), -model_extension_width)
    print(third_column)
    missing_column = 0

    if 'x' not in cmodel.settings.csv_columns:
        points = np.column_stack((third_column,points))
    elif 'y' not in cmodel.settings.csv_columns:
        column_x = points[:,0]
        column_z = points[:,1]
        points = np.column_stack((column_x,third_column))
        points = np.column_stack((points,column_z))
        missing_column = 1
    elif 'z' not in cmodel.settings.csv_columns:
        points = np.column_stack((points,third_column))
        missing_column = 2

    print(missing_column)
    # add to data a copy of the existing data with shifted 3rd coordinate to 100
    tmp_points = np.copy(points)
    tmp_points[:,missing_column] += model_extension_width
    print(points)
    print(tmp_points)
    return np.vstack((points,tmp_points))


def structure_in_model_boundaries(cmodel, points):
    '''
        Check if structure is in model boundaries
    '''
    # print(cmodel.settings.model_b1)
    # print(cmodel.settings.model_b2)
    # print(points)

    # np.save("test", points)
    # sys.exit()

    return np.any(
        (cmodel.settings.model_b1 <= points).all(axis=1) &
        (points <= cmodel.settings.model_b2).all(axis=1)
    )


def load_structure_files(cmodel, structure_type='layer', vtk=True):
    '''
        Load text files containing layer, fault or seam point cloud data
    '''
    layers = Layer_Set()
    layers.A = {}
    layers.B = {}
    layers.orientation = {}
    i = 0

    print(f"\nProcess point {structure_type} cloud files:\n")

    filenamelist = sorted(
        glob.glob(cmodel.settings.path_input_data + structure_type + '-*.csv')
    )


    for filename in filenamelist:
        filename = os.path.basename(filename)

        #### load xyz raster file
        points = load_and_organize_point_cloud_data(
                cmodel,
                cmodel.settings.path_input_data + filename
        )

        pointsToVTK(
            cmodel.settings.path_output_data + "/" + \
            os.path.splitext(filename)[0] + "_points",
            np.ascontiguousarray(points[:,0]),
            np.ascontiguousarray(points[:,1]),
            np.ascontiguousarray(points[:,2])
        )

        # for 2d csv data
        if len(cmodel.settings.csv_columns) == 2:
            points = convert_2d_data_to_3d_data(cmodel, points)

        # shift csv data to origin
        if cmodel.settings.model_p1 != [0,0,0]:
            points =  np.array(
                [
                    shift_coordinate(points[i,:], cmodel.settings.model_p1)\
                        for i in range(len(points))
                ]
            )

        # rotate csv data so that the data are parallel to x axis
        if cmodel.settings.rotation_angle != 0:
            points = np.array(
                [
                    rotate_coordinate(cmodel, points[i,:], clockwise=1)\
                        for i in range(len(points))
                ]
            )

        # for 2d+ model extend 2d csv, from a point line to a surface
        # the existing points are duplicated and
        # the first copy is shifted on y axis to -50*abs(model width/2) and
        # the other copy is shifted on y asis to 50*abs(model width/2)
        if not cmodel.settings.s3d:
            model_extension_width = cmodel.settings.model_dimension[1]/2+50
            # print(model_extension_width)

            points[:,1] -= model_extension_width
            tmp_points = np.copy(points)
            tmp_points[:,1] += model_extension_width*2
            points = np.vstack((points,tmp_points))

        # if structure_in_model_boundaries(cmodel, points):
        np.savetxt(
            cmodel.settings.path_input_data + os.path.splitext(
                filename
            )[0] + ".tmp",
            points,
            fmt='%.2f',
            header='x,y,z',
            comments="",
            delimiter=","
        )
        # else:
            # print(f"> {filename} is ignored because coordinates out of model boundaries")


    filenamelist = sorted(
        glob.glob(cmodel.settings.path_input_data + structure_type + '-*.tmp')
    )

    if filenamelist:
        print(f"\nBest plane for interpolation of {structure_type}" + \
            " values to the grid:\n")

    for filename in filenamelist:

        filename = os.path.basename(filename)
        #### load xyz raster file
        points = load_point_cloud_files(
            cmodel.settings.path_input_data + filename
        )

        ## set NAN to 0
        points[np.isnan(points)] = 0

        ex = np.array(cmodel.x)
        if not cmodel.settings.s3d:
            ey = np.asarray([-model_extension_width*3/4, model_extension_width*3/4])
        else:
            ey = np.array(cmodel.y)
        ez = np.array(cmodel.z)


        layers.orientation[i]= "xy"
        elevationXY = interpolate(points, ex, ey, method=2)

        elevationXZ = None
        if cmodel.ny > 1:
            elevationXZ = interpolate(
                np.array([points[:,0], points[:,2], points[:,1]]).T,
                ex, ez, method=2
            )

        elevationYZ = interpolate(
            np.array([points[:,1], points[:,2], points[:,0]]).T,
            ey, ez, method=2
        )


        if elevationXY is None:
            layers.orientation[i]= "yz"

            if elevationYZ is None:
                layers.orientation[i]= "xz"

                if cmodel.ny > 1 and elevationXZ is None:
                    print(' coordinates of ' + structure_type +\
                            ' could not be interpolated and will be ignored!\n')
                    continue

        elevationXY_quality = 0
        elevationXZ_quality = 0
        elevationYZ_quality = 0

        if elevationXY is not None:
            elevationXY_quality = np.sum(~np.isnan(elevationXY))/\
                (cmodel.nx*cmodel.ny)
        if cmodel.ny > 1 and elevationXZ is not None:
            elevationXZ_quality = np.sum(~np.isnan(elevationXZ))/\
                (cmodel.nx*cmodel.nz)
        if elevationYZ is not None:
            elevationYZ_quality = np.sum(~np.isnan(elevationYZ))/\
                (cmodel.ny*cmodel.nz)

        interpolation_quality = elevationXY_quality

        if cmodel.ny > 1 and interpolation_quality < elevationXZ_quality:
            layers.orientation[i]= "xz"
            interpolation_quality = elevationXZ_quality
        if interpolation_quality < elevationYZ_quality:
            layers.orientation[i]= "yz"
            interpolation_quality = elevationXZ_quality

        width = 0
        if filename.lower().split(".")[0] in cmodel.settings.partition_width:
            width = cmodel.settings.partition_width[filename.split(".")[0]]

        # np.set_printoptions(threshold=sys.maxsize)
        if layers.orientation[i] == "xy":
            # print("---------XY------------")
            # print(elevationXY)
            a, b = np.meshgrid(ex, ey)
            c = elevationXY
            layers.A[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationXY
            }
            layers.B[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationYZ
            }

        if layers.orientation[i] == "xz":
            # print("---------XZ------------")
            # print(elevationXZ)
            a, c = np.meshgrid(ex, ez)
            b = elevationXZ
            layers.A[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationXZ
            }
            layers.B[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationYZ
            }

        if layers.orientation[i] == "yz":
            # print("---------YZ------------")
            # print(elevationYZ)
            b, c = np.meshgrid(ey, ez)
            a = elevationYZ
            layers.A[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationYZ
            }
            layers.B[i] = {
                "name" : filename.split(".")[0],
                "width": width,
                "elevation": elevationXY
            }

        if vtk:
            print(f"\tinterpolation plane {layers.orientation[i]}\t->\t",end="")

            vtk_filename = cmodel.settings.path_output_data + \
                os.path.splitext(filename)[0]

            status = points_to_pyvista_surface(
                cmodel,
                vtk_filename,
                a.flatten(),
                b.flatten(),
                c.flatten()
            )

            # status = points_to_vtk_surface(
                # cmodel,
                # vtk_filename,
                # a.flatten(),
                # b.flatten(),
                # c.flatten()
            # )

            if status:
                print(f"{vtk_filename}.vtk")
            else:
                print("\n\tWARNING: Interpolation of surface not successful." +\
                     f"\n\tBuilding of {vtk_filename}.vtu not possible!")

        os.remove(cmodel.settings.path_input_data + \
            os.path.splitext(filename)[0] + '.tmp')

        i += 1

    return layers


def generate_layer_partitions(cmodel, layers):
    '''
        Generate model cell clusters regarding layers
        and add a partition number to this cells
    '''
    grid_data = np.ones((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    tmodel = None
    flag = 0

    print("\n====================================")
    print("\tGenerate layer partitions")
    print("====================================\n")
    print("Total number of partitions: ", end='')

    for l in range(0, len(layers.A)):

        if layers.orientation[l] == "xy":
            ijk = 2
        elif layers.orientation[l] == "xz":
            ijk = 1
        else:
            ijk = 0

        tmodel = np.meshgrid(cmodel.x, cmodel.y, cmodel.z, indexing='ij')[ijk]

        for i in range(0, cmodel.nx):
            for j in range(0, cmodel.ny):
                for k in range(0, cmodel.nz):
                    if layers.orientation[l] == "xy":
                        if layers.A[l]['elevation'][j, i] > tmodel[i, j, k]:
                            grid_data[k, j, i] += l+1
                            cmodel.settings.partition_number = max(
                                [
                                    cmodel.settings.partition_number,
                                    grid_data[k,j,i]
                                ]
                            )
                            flag = 1

                    elif layers.orientation[l] == "xz":
                        if layers.A[l]['elevation'][k, i] > tmodel[i, j, k]:
                            grid_data[k, j, i] += l+1
                            cmodel.settings.partition_number = max(
                                [
                                    cmodel.settings.partition_number,
                                    grid_data[k,j,i]
                                ]
                            )
                            flag = 1

                    elif layers.orientation[l] == "yz":
                        if layers.A[l]['elevation'][k, j] > tmodel[i, j, k]:
                            grid_data[k, j, i] += l+1
                            cmodel.settings.partition_number = max(
                                [
                                    cmodel.settings.partition_number,
                                    grid_data[k,j,i]
                                ]
                            )
                            flag = 1

    print(str(cmodel.settings.partition_number+1), end='')
    if flag:
        print(' (0-' + str(cmodel.settings.partition_number) + ')', end='')
    print('')

    return grid_data


def generate_model_cell_ids(cmodel):
    '''
        Generate model cell ids for model post processing
    '''
    grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    grid_i_ids = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    grid_j_ids = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    grid_k_ids = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)

    for k in range(0, cmodel.nz):
        for j in range(0, cmodel.ny):
            for i in range(0, cmodel.nx):
                grid_data[k, j, i] = i + j*cmodel.nx + k*cmodel.nx*cmodel.ny
                grid_i_ids[k, j, i] = i
                grid_j_ids[k, j, i] = j
                grid_k_ids[k, j, i] = k

    return grid_data, grid_i_ids, grid_j_ids, grid_k_ids


def generate_active_partitions(cmodel,old_grid_data=None,name="active"):
    '''
        Generate model cell clusters regarding the shapefile mask 
        and add a patiation number to this cells
    '''

    filenamelist = sorted(
        glob.glob(cmodel.settings.path_shape_data + '*.shp')
    )

    special_grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    tiff_file = cmodel.settings.path_output_data + cmodel.settings.mask_file

    if len(filenamelist):
        print("\n====================================")
        print("\tGenerate " + name + " partitions")
        print("====================================\n")
        print("-> " + filenamelist[-1])

        shape = ogr.Open(filenamelist[-1])

        shape_layer = shape.GetLayer()
        x_min, x_max, y_min, y_max = shape_layer.GetExtent()

        x_pixel_size = int(np.min(cmodel.settings.dx)/2)
        y_pixel_size = int(np.min(cmodel.settings.dy)/2)

        x_pixel_res = int((x_max-x_min)/x_pixel_size)
        y_pixel_res = int((y_max-y_min)/y_pixel_size)

        target_ds = gdal.GetDriverByName('GTiff').Create(
            tiff_file, x_pixel_res, y_pixel_res, 1, gdal.GDT_Byte
        )

        target_ds.SetGeoTransform(
            (x_min, x_pixel_size, 0, y_min, 0, y_pixel_size)
        )

        band = target_ds.GetRasterBand(1)
        NoData_value = 0
        band.SetNoDataValue(NoData_value)
        band.FlushCache()
        gdal.RasterizeLayer(target_ds, [1], shape_layer)
        target_ds = None
        mask = gdal.Open(tiff_file).ReadAsArray()

        x_mask_vector = np.array(
            [
                x_min+((x_max-x_min)/x_pixel_res)*i
                for i in range(0, x_pixel_res)
            ]
        )
        y_mask_vector = np.array(
            [
                y_min+((y_max-y_min)/y_pixel_res)*i
                for i in range(0, y_pixel_res)
            ]
        )

        x_mask, y_mask = np.meshgrid(x_mask_vector, y_mask_vector)
        mask_data = np.array(
            (x_mask.flatten(), y_mask.flatten(), mask.flatten())
        )

        x_vector_model = cmodel.x[:] + x_min
        y_vector_model = cmodel.y[:] + y_min
        model_mask = interpolate(
            mask_data.T, x_vector_model, y_vector_model[::-1]
        )

        old_grid_data[:,model_mask[::-1]< 90] = 0
        special_grid_data[:,model_mask[::-1]>=90] = 1

        cmodel.settings.partition_number += 1

    return old_grid_data, special_grid_data



def identify_partition_elements(
    cmodel,
    grid_data,
    special_grid_data,
    a,
    layer_orientation,
    layer,
    partition_number
):
    '''
        identify elements that are part of partition
    '''
    model = np.meshgrid(cmodel.x, cmodel.y, cmodel.z, indexing='ij')[a.vd]

    for i in range(0, cmodel.nx - a.vi):
        for j in range(0, cmodel.ny - a.vj):
            for k in range(0, cmodel.nz - a.vk):

                layer_elevation = 0

                if layer_orientation == "xy":
                    try:
                        layer_elevation = layer["elevation"][j, i]
                    except:
                        layer_elevation = layer["elevation"][k, j]
                elif layer_orientation == "xz":
                    try:
                        layer_elevation = layer["elevation"][k, i]
                    except:
                        layer_elevation = layer["elevation"][k, j]
                else:
                    try:
                        layer_elevation = layer["elevation"][k, j]
                    except:
                        layer_elevation = layer["elevation"][j, i]

                nk = k + a.vk
                nj = j + a.vj
                ni = i + a.vi

                if model[i, j, k] <= layer_elevation < model[ni,nj,nk] and \
                    layer["elevation"] is not None:

                    grid_data[k,j,i] = cmodel.settings.partition_number+1
                    grid_data[nk,nj,ni] = cmodel.settings.partition_number+1

                    if special_grid_data[k, j, i] == 0:
                        special_grid_data[k, j, i] = partition_number
                    if special_grid_data[nk,nj,ni] == 0:
                        special_grid_data[nk,nj,ni] = partition_number

                    t = 1
                    ptk = k - t*a.vk
                    ptj = j - t*a.vj
                    pti = i - t*a.vi

                    valid_offset = 0 <= (i - t)*a.vi and \
                        0 <= (j - t)*a.vj and 0 <= (k - t)*a.vk

                    while valid_offset and\
                        layer_elevation - layer["width"] < model[pti,ptj,ptk]:

                        grid_data[ptk,ptj,pti] = \
                            cmodel.settings.partition_number + 1

                        if special_grid_data[ptk,ptj,pti] == 0:
                            special_grid_data[ptk,ptj,pti] = partition_number

                        t += 1
                        ptk = k - t*a.vk
                        ptj = j - t*a.vj
                        pti = i - t*a.vi

                    t = 1
                    ntk = k + t*a.vk
                    ntj = j + t*a.vj
                    nti = i + t*a.vi

                    valid_offset = (i + t)*a.vi + (j + t)*a.vj + (k + t)*a.vk <\
                        (cmodel.nx-1)*a.vi+(cmodel.ny-1)*a.vj+(cmodel.nz-1)*a.vk

                    while valid_offset and\
                        layer_elevation + layer["width"] > model[nti,ntj,ntk]:

                        grid_data[ntk,ntj,nti] = \
                            cmodel.settings.partition_number + 1

                        if special_grid_data[ntk,ntj,nti] == 0:
                            special_grid_data[ntk,ntj,nti] = partition_number

                        t += 1
                        ntk = k + t*a.vk
                        ntj = j + t*a.vj
                        nti = i + t*a.vi

    return grid_data, special_grid_data


class Check_Control():
    '''
        Data structure to hold the moving vector for the parition check
    '''
    def __init__(self, vd, vi, vj, vk):
        ## vector direction
        self.vd = vd
        ## vector i switch
        self.vi = vi
        ## vector i switch
        self.vj = vj
        ## vector i switch
        self.vk = vk


def generate_special_partitions(cmodel,layers,old_grid_data=None,name="fault"):
    '''
        Generate model cell clusters regarding faults or seams
        and add a patiation number to this cells
    '''
    SPECIALpartition_number = 1

    grid_data = np.copy(old_grid_data)
    if grid_data is None:
        grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)

    special_grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)

    print("\n====================================")
    print("\tGenerate " + name + " partitions")
    print("====================================\n")
    print(name + "   \t| partition number")
    print("------------------------------------")

    for l in range(0, len(layers.A)):

        ## set direction for partion checks
        # for every layer two checks are performed "a" and "b"
        if layers.orientation[l] == "xy":
            # in case of structure (layer/fautl/seam) surface is orientated
            # in xy orientation
            ## a do check in z (k = 2) direction
            a = Check_Control(2, 0, 0, 1)
            ## b do check in x (i = 0) direction
            b = Check_Control(0, 1, 0, 0)
        elif layers.orientation[l] == "xz":
            # in case of structure (layer/fautl/seam) surface is orientated
            # in xz orientation
            ## a do check in y (j = 1) direction
            a = Check_Control(1, 0, 1, 0)
            ## b do check in x (i = 0) direction
            b = Check_Control(2, 1, 0, 0)
        else:
            # in case of structure (layer/fautl/seam) surface is orientated
            # in yz orientation
            ## a do check in x (i = 0) direction
            a = Check_Control(0, 1, 0, 0)
            ## b do check in x (i = 0) direction
            b = Check_Control(2, 1, 0, 0)


        grid_data, special_grid_data = identify_partition_elements(
            cmodel,
            grid_data,
            special_grid_data,
            a,
            layers.orientation[l],
            layers.A[l],
            SPECIALpartition_number
        )

        grid_data, special_grid_data = identify_partition_elements(
            cmodel,
            grid_data,
            special_grid_data,
            b,
            layers.orientation[l],
            layers.B[l],
            SPECIALpartition_number
        )

        print(layers.A[l]["name"] + "\t| " + \
            str(cmodel.settings.partition_number + 1))
        SPECIALpartition_number += 1
        cmodel.settings.partition_number += 1

    return grid_data, special_grid_data


def points_to_pyvista_surface(cmodel, filename, x, y, z):
    '''
        Create vtk output of layers, faults or seams
    '''

    status = True
    points = np.zeros( (len(x), 3) )
    for i, x_value in enumerate(x):
        points[i,0] = x_value
        points[i,1] = y[i]
        points[i,2] = z[i]
    points = points[~np.isnan(points).any(axis=1),:]

    if cmodel.settings.rotation_angle != 0:
        points = np.array(
            [rotate_coordinate(cmodel, points[i,:]) for i in range(len(points))]
        )

    if cmodel.settings.model_p1 != [0,0,0]:
        points =  np.array(
            [
                shift_coordinate(
                    points[i,:], [-value for value in cmodel.settings.model_p1]
                ) for i in range(len(points))
            ]
        )

    points = np.array(
        [
            i
                if i[2]<cmodel.settings.model_p8[2]
                else [np.nan, np.nan, np.nan]
            for i in points
        ]
    )

    points = np.array(
        [
            i
                if i[2]>cmodel.settings.model_p1[2]
                else [np.nan, np.nan, np.nan]
            for i in points
        ]
    )
    points = points[~np.isnan(points).any(axis=1),:]

    vtk_filename = filename + ".vtk"
    pcl = pv.PolyData(points)
    mesh = pcl.delaunay_2d(tol=1e-06)
    mesh.save(vtk_filename)

    return status



def points_to_vtk_surface(cmodel, filename, x, y, z):
    '''
        Create vtk output of layers, faults or seams
    '''

    points = np.zeros( (len(x), 3) )
    for i, x_value in enumerate(x):
        points[i,0] = x_value
        points[i,1] = y[i]
        points[i,2] = z[i]
    points = points[~np.isnan(points).any(axis=1),:]

    if cmodel.settings.rotation_angle != 0:
        points = np.array(
            [rotate_coordinate(cmodel, points[i,:]) for i in range(len(points))]
        )

    if cmodel.settings.model_p1 != [0,0,0]:
        points =  np.array(
            [
                shift_coordinate(
                    points[i,:], [-value for value in cmodel.settings.model_p1]
                ) for i in range(len(points))
            ]
        )

    try:
        # 2d surface
        if not cmodel.settings.s3d:
            triangles = Delaunay(points,qhull_options="QJ Pp Qw")
        # 3d surface
        else:
            triangles = Delaunay(points,qhull_options="Qbb Qt Qc Qz Q12")

        # list of triangles that form the tesselation
        ncells = triangles.simplices.shape[0]

        connections =  np.zeros(ncells * 3)
        for i in range(ncells):
            ii = i * 3
            connections[ii]     = triangles.simplices[i,0]
            connections[ii + 1] = triangles.simplices[i,1]
            connections[ii + 2] = triangles.simplices[i,2]

        offset = np.zeros(ncells)
        for i in range(ncells):
            offset[i] = (i + 1) * 3

        # vtk triangle cell type id = 5
        cell_type = np.ones(ncells) * 5

        for i in range(len(points)):
            x[i] = points[i,0]
            y[i] = points[i,1]
            z[i] = points[i,2]


        unstructuredGridToVTK(
            filename + "_qhull" , x, y, z,
            connectivity = connections,
            offsets = offset,
            cell_types = cell_type,
            cellData = None,
            pointData = {"Elevation" : z}
        )
        status = True
    except ValueError:
        status = False

    return status


def print_help():
    '''
        Help text
    '''
    print('\n\tPlease use this command to run GEOMODELATOR:\n' + \
         '\033[1m\tpython [-i] [-u] [path]geomodelator.py ' + \
         '[path]<config file name> ' + \
         '[2>&1 | tee [path]<log file name>]\033[0m' + \
         '\n\n\t\033[1mExamples:\033[0m' + \
         '\n\tpython -u geomodelator.py examples/2d/config.py 2>&1 ' + \
         '| tee config.log' + \
         '\n\tpython -i geomodelator.py examples/3d/config.py')
    sys.exit()
