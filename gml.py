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

from pyevtk.hl import unstructuredGridToVTK


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
    partition_width = None
    dx = None
    dy = None
    dz = None

    # global cache variable
    rotation_angle = 0
    partition_number = 0

    def __init__(
        self,
        s3d,
        csv_columns,
        model_p1,
        model_p2,
        partition_width,
        dx,
        dy,
        dz
    ):
        self.s3d = s3d
        self.csv_columns = csv_columns
        self.model_p1 = model_p1
        self.model_p2 = model_p2
        self.partition_width = partition_width
        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.path_input_data = "1_input/"
        self.path_output_data = "2_output/"


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

            print("")
            print("====================================")
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
            print("                                               ")
            print("")
            print("Origin coordinate (x, y, z): (", end="")
            print(f"{self.settings.model_p1[0]},", end="")
            print(f"{self.settings.model_p1[1]},", end="")
            print(f"{self.settings.model_p1[2]})", end="")
            print("")
            print("Discretization (nx, ny, nz): (", end="")
            print(f"{len(self.settings.dx)}, ", end="")
            print(f"{len(self.settings.dy)}, ", end="")
            print(f"{len(self.settings.dz)})")
            print("")
            print("Normalized model corner coordinates (x, y, z):")
            print(f"\tP1 ({self.x[0]}, {self.y[0]}, {self.z[0]})")
            print(f"\tP2 ({self.x[-1]}, {self.y[0]}, {self.z[0]})")
            print(f"\tP3 ({self.x[-1]}, {self.y[-1]}, {self.z[0]})")
            print(f"\tP4 ({self.x[0]}, {self.y[-1]}, {self.z[0]})")
            print(f"\tP5 ({self.x[0]}, {self.y[0]}, {self.z[-1]})")
            print(f"\tP6 ({self.x[-1]}, {self.y[0]}, {self.z[-1]})")
            print(f"\tP7 ({self.x[-1]}, {self.y[-1]}, {self.z[-1]})")
            print(f"\tP8 ({self.x[0]}, {self.y[-1]}, {self.z[-1]})")
            print("")
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

    def calculate_cell_centers(self):
        '''
            Method to calculate cell/element center coordinates.
        '''
        return (self.x[1:] + self.x[:-1])/2, (self.y[1:] + self.y[:-1])/2, \
                (self.z[1:] + self.z[:-1])/2


class Scattered_Surface:
    '''
        Scattered_Surface class
    '''
    coords = None
    elevation = None
    name = ""

    #### initialize
    def __init__(self, coords=None, elevation=None, data=None, name="default"):
        '''
            Constructor to the scattered surface class.
        '''

        if data is not None:
            self.coords = (data[:,0], data[:,1])
            self.elevation = data[:,2]
        else:
            self.coords = coords
            self.elevation = elevation

        self.name = name

    def interpolate(self, x, y, method=1):
        '''
            Method to interpolate a surface to plane (e.g. xy or yz or xz).
        '''

        interp_method = ['nearest', 'linear', 'cubic']
        gx, gy = np.meshgrid(x, y)

        try:
            grid = griddata(
                self.coords,
                self.elevation,
                (gx, gy),
                method=interp_method[method]
            )
        except ValueError:
            grid = None

        return grid

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
            print(filename)
            data = np.genfromtxt(filename, skip_header=cmodel.settings.csv_header,\
                delimiter=cmodel.settings.csv_delimiter, usecols=tuple(columns))

        except ValueError:
            print('Please check if the csv file column number and the csv\
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
    return math.radians(math.atan2(P[1] - origin[1], P[0] - origin[0]) *\
        (180.0 / math.pi))


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


def rotate_grid(model):
    '''
        Rotate model by rotation angle
    '''
    z_matrix, y_matrix, x_matrix = np.meshgrid(model.z, model.y, model.x,\
                                indexing='ij')

    z_list = [item for layer in z_matrix for row in layer for item in row]
    y_list = [item for layer in y_matrix for row in layer for item in row]
    x_list = [item for layer in x_matrix for row in layer for item in row]

    coord_list = np.array([[x_list[i],y_list[i],z_list[i]]\
                        for i in range(len(x_list))])

    if model.settings.rotation_angle != 0:
        coord_list = np.array([rotate_coordinate(model, coord_list[i,:], clockwise=-1)\
                                for i in range(len(coord_list))])

    if model.settings.model_p1 != [0,0,0]:
        coord_list =  np.array([shift_coordinate(coord_list[i,:],\
                                [-value for value in model.settings.model_p1])\
                                for i in range(len(coord_list))])


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
    model_extension_width = math.dist(cmodel.settings.model_p2, cmodel.settings.model_p1)/2+50

    third_column = np.full(len(points[:,0]), -model_extension_width)
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

    # add to data a copy of the existing data with shifted 3rd coordinate to 100
    tmp_points = np.copy(points)
    tmp_points[:,missing_column] += model_extension_width*2
    return np.vstack((points,tmp_points))


def load_structure_files(cmodel, structure_type='layer', vtk=True):
    '''
        Load text files containing layer, fault or seam point cloud data
    '''
    layers = Layer_Set()
    layers.A = {}
    layers.B = {}
    layers.orientation = {}
    i = 0

    if cmodel.settings.rotation_angle == 0:
        cmodel.settings.rotation_angle = calculate_full_angle(
            cmodel.settings.model_p2,
            cmodel.settings.model_p1
        )

        print("")
        print(f"Rotation angle: {np.degrees(cmodel.settings.rotation_angle).round()}")
        print("")

    filenamelist = sorted(glob.glob(cmodel.settings.path_input_data + structure_type + '-*.csv'))

    for filename in filenamelist:
        filename = os.path.basename(filename)

        #### load xyz raster file
        points = load_and_organize_point_cloud_data(
                cmodel,
                cmodel.settings.path_input_data + filename
        )

        # for 2d csv data
        if len(cmodel.settings.csv_columns) == 2:
            points = convert_2d_data_to_3d_data(cmodel, points)

        # shift csv data to origin
        if cmodel.settings.model_p1 != [0,0,0] :
            points =  np.array([shift_coordinate(points[i,:], cmodel.settings.model_p1)\
                            for i in range(len(points))])

        # rotate csv data so that the data are parallel to x axis
        if cmodel.settings.rotation_angle != 0:
            points = np.array([rotate_coordinate(cmodel, points[i,:], clockwise=1)\
                            for i in range(len(points))])

        # for 2d+ model extend 2d csv, from a point line to a surface
        # the existing points are duplicated and
        # the first copy is shifted on y axis to -50*abs(model width/2) and
        # the other copy is shifted on y asis to 50*abs(model width/2)
        if cmodel.settings.s3d is False:
            model_extension_width = cmodel.settings.dy.sum()/2+50

            points[:,1] -= model_extension_width
            tmp_points = np.copy(points)
            tmp_points[:,1] += model_extension_width*2
            points = np.vstack((points,tmp_points))


        np.savetxt(
            cmodel.settings.path_input_data + os.path.splitext(filename)[0] + ".tmp",
            points,
            fmt='%.2f',
            header='x,y,z',
            comments="",
            delimiter=","
        )


    filenamelist = sorted(glob.glob(cmodel.settings.path_input_data + structure_type + '-*.tmp'))

    if filenamelist:
        print("")
        print(f"Best plane for interpolation of {structure_type} values to the grid:")
        print("")

    for filename in filenamelist:

        filename = os.path.basename(filename)
        #### load xyz raster file
        points = load_point_cloud_files(cmodel.settings.path_input_data + filename)

        ## set NAN to 0
        points[np.isnan(points)] = 0

        surfaceXY = Scattered_Surface(data=points)
        surfaceXZ = Scattered_Surface(
            data=np.array([points[:,0], points[:,2], points[:,1]]).T
        )
        surfaceYZ = Scattered_Surface(
            data=np.array([points[:,1], points[:,2], points[:,0]]).T
        )

        ex = 0
        if len(cmodel.x) < 2:
            ex = np.array([cmodel.x-10, cmodel.x + 10])
        else:
            ex = np.array([cmodel.x, cmodel.x])

        ey = 0
        if len(cmodel.y) < 2:
            ey = np.array([cmodel.y-10, cmodel.y + 10])
        else:
            ey = np.array([cmodel.y, cmodel.y])

        ez = 0
        if len(cmodel.z) < 2:
            ez = np.array([cmodel.z-10, cmodel.z + 10])
        else:
            ez = np.array([cmodel.z, cmodel.z])

        layers.orientation[i]= "xy"
        elevationXY = surfaceXY.interpolate(ex, ey, 1)

        if cmodel.ny > 1:
            elevationXZ = surfaceXZ.interpolate(ex, ez, 1)

        elevationYZ = surfaceYZ.interpolate(ey, ez, 1)

        if elevationXY is None:
            layers.orientation[i]= "yz"

            if elevationYZ is None:
                layers.orientation[i]= "xz"

                if cmodel.ny > 1 and elevationXZ is None:
                    print(' the point cloud of this ' + structure_type +\
                                ' could not be interpolated!\n')
                    continue

        elevationXY_quality = 0
        elevationXZ_quality = 0
        elevationYZ_quality = 0

        if elevationXY is not None:
            elevationXY_quality = np.sum(~np.isnan(elevationXY))
        if cmodel.ny > 1 and elevationXZ is not None:
            elevationXZ_quality = np.sum(~np.isnan(elevationXZ))
        if elevationYZ is not None:
            elevationYZ_quality = np.sum(~np.isnan(elevationYZ))

        quality = elevationXY_quality

        if cmodel.ny > 1 and quality < elevationXZ_quality:
            layers.orientation[i]= "xz"
            quality = elevationXZ_quality

        if quality < elevationYZ_quality:
            layers.orientation[i]= "yz"
            quality = elevationYZ_quality


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
            print(f"\tinterpolation plane {layers.orientation[i]}\t->\t", end="")

            vtk_filename = cmodel.settings.path_output_data + os.path.splitext(filename)[0]

            status = points_to_vtk_surface(cmodel, vtk_filename, a.flatten(), b.flatten(),\
                                    c.flatten())
            if status:
                print(f"{vtk_filename}.vtu")
                os.remove(cmodel.settings.path_input_data + os.path.splitext(filename)[0] + '.tmp')
            else:
                print("\n\tWARNING: Interpolation of surface not successful.\n\t" +
                        f"Building of {vtk_filename}.vtu not posible!")
        i += 1

    # print("\n---------------------\nnumber of layers: %s\n" % i)

    return layers


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


def generate_layer_partitions(cmodel, layers):
    '''
        Generate model cell clusters regarding layers
        and add a partition number to this cells
    '''
    grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)
    tmodel = None
    flag = 0

    print("")
    print("====================================")
    print("\tGenerate layer partitions")
    print("====================================")
    print("")
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
                                [cmodel.settings.partition_number,grid_data[k,j,i]]
                            )
                            flag = 1

                    elif layers.orientation[l] == "xz":
                        if layers.A[l]['elevation'][k, i] > tmodel[i, j, k]:
                            grid_data[k, j, i] += l+1
                            cmodel.settings.partition_number = max(
                                [cmodel.settings.partition_number, grid_data[k,j,i]]
                            )
                            flag = 1

                    elif layers.orientation[l] == "yz":
                        if layers.A[l]['elevation'][k, j] > tmodel[i, j, k]:
                            grid_data[k, j, i] += l+1
                            cmodel.settings.partition_number = max(
                                [cmodel.settings.partition_number, grid_data[k,j,i]]
                            )
                            flag = 1

    print(str(cmodel.settings.partition_number+1), end='')
    if flag:
        print(' (0-' + str(cmodel.settings.partition_number) + ')', end='')
    print('')

    return grid_data


def generate_special_partitions(cmodel,layers,old_grid_data=None,name="fault"):
    '''
        Generate model cell clusters regarding faults or seams
        and add a patiation number to this cells
    '''
    SPECIALpartition_number = 1


    ## debug switch: high light elements left and right of a fault/seam in
    ## different colors
    hl = 0
    hl_offset = 8

    grid_data = np.copy(old_grid_data)
    if grid_data is None:
        grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)

    special_grid_data = np.zeros((cmodel.nz, cmodel.ny, cmodel.nx), np.int32)

    print("")
    print("")
    print("====================================")
    print("\tGenerate " + name + " partitions")
    print("====================================")
    print("")
    print(name + "   \t| partition number")
    print("------------------------------------")

    for l in range(0, len(layers.A)):

        if layers.orientation[l] == "xy":
            aijk = 2
            ai = 0
            aj = 0
            ak = 1
            bijk = 0
            bi = 1
            bj = 0
            bk = 0
        elif layers.orientation[l] == "xz":
            aijk = 1
            ai = 0
            aj = 1
            ak = 0
            bijk = 2
            bi = 1
            bj = 0
            bk = 0
        else:
            aijk = 0
            ai = 1
            aj = 0
            ak = 0
            bijk = 2
            bi = 0
            bj = 0
            bk = 1


        atmodel = np.meshgrid(cmodel.x, cmodel.y, cmodel.z, indexing='ij')[aijk]
        btmodel = np.meshgrid(cmodel.x, cmodel.y, cmodel.z, indexing='ij')[bijk]

        for i in range(0, cmodel.nx-ai):
            for j in range(0, cmodel.ny-aj):
                for k in range(0, cmodel.nz-ak):

                    layerA_elevation = 0
                    layerB_elevation = 0

                    if layers.orientation[l] == "xy":
                        layerA_elevation = layers.A[l]["elevation"][j, i]
                        if layers.B[l]["elevation"] is not None:
                            layerB_elevation = layers.B[l]["elevation"][k, j]
                    elif layers.orientation[l] == "xz":
                        layerA_elevation = layers.A[l]["elevation"][k, i]
                        if layers.B[l]["elevation"] is not None:
                            layerB_elevation = layers.B[l]["elevation"][k, j]
                    else:
                        layerA_elevation = layers.A[l]["elevation"][k, j]
                        if layers.B[l]["elevation"] is not None:
                            layerB_elevation = layers.B[l]["elevation"][j, i]


                    if atmodel[i, j, k] <= layerA_elevation < atmodel[i+ai, j+aj, k+ak]:

                        grid_data[k, j, i] = cmodel.settings.partition_number + 1
                        grid_data[k+ak, j+aj, i+ai] = cmodel.settings.partition_number + 1

                        if special_grid_data[k, j, i] == 0:
                            special_grid_data[k, j, i] = SPECIALpartition_number + hl*0
                        if special_grid_data[k+ak, j+aj, i+ai] == 0:
                            special_grid_data[k+ak, j+aj, i+ai] = SPECIALpartition_number + hl*1

                        a = 1
                        while 0 <= (i-a)*ai + (j-a)*aj + (k-a)*ak and\
                            layerA_elevation - layers.A[l]["width"] <\
                            atmodel[i-a*ai, j-a*aj, k-a*ak]:

                            grid_data[k-a*ak, j-a*aj, i-a*ai] = cmodel.settings.partition_number + 1

                            if special_grid_data[k-a*ak, j-a*aj, i-a*ai] == 0:
                                special_grid_data[k-a*ak, j-a*aj, i-a*ai] = SPECIALpartition_number\
                                + hl*2

                            a += 1

                        # try:
                        a = 1
                        while (i+a)*ai + (j+a)*aj + (k+a)*ak <\
                            (cmodel.nx-1)*ai + (cmodel.ny-1)*aj + (cmodel.nz-1)*ak and\
                            layerA_elevation + layers.A[l]["width"] >\
                            atmodel[i+a*ai, j+a*aj, k+a*ak]:

                            grid_data[k+a*ak, j+a*aj, i+a*ai] = cmodel.settings.partition_number + 1

                            if special_grid_data[k+a*ak, j+a*aj, i+a*ai] == 0:
                                special_grid_data[k+a*ak, j+a*aj, i+a*ai] = SPECIALpartition_number\
                                + hl*2

                            a += 1

                    if layers.B[l]["elevation"] is not None and \
                        i*bi + j*bj + k*bk < (cmodel.nx-1)*bi + (cmodel.ny-1)*bj + (cmodel.nz-1)*bk\
                        and btmodel[i, j, k] <= layerB_elevation < btmodel[i+1*bi, j+1*bj, k+1*bk]:

                        grid_data[k, j, i] = cmodel.settings.partition_number + 1
                        grid_data[k+bk, j+bj, i+bi] = cmodel.settings.partition_number + 1

                        if special_grid_data[k, j, i] == 0:
                            special_grid_data[k, j, i] = SPECIALpartition_number + hl*4
                        if special_grid_data[k+bk, j+bj, i+bi] == 0:
                            special_grid_data[k+bk, j+bj, i+bi] = SPECIALpartition_number + hl*1

                            a = 1
                            while (i+a)*bi + (j+a)*bj + (k+a)*bk <\
                                (cmodel.nx-1)*bi + (cmodel.ny-1)*bj + (cmodel.nz-1)*bk and\
                                layerB_elevation + layers.B[l]["width"] >\
                                btmodel[i+a*bi, j+a*bj, k+a*bk]:

                                grid_data[k+a*bk,j+a*bj,i+a*bi] = cmodel.settings.partition_number+1

                                if special_grid_data[k+a*bk, j+a*bj, i+a*bi] == 0:
                                    special_grid_data[k+a*bk,j+a*bj,i+a*bi] =\
                                        SPECIALpartition_number + hl*6

                                a += 1

                            a = 1
                            while 0 <= (i-a)*bi + (j-a)*bj + (k-a)*bk and\
                                layerB_elevation - layers.B[l]["width"]\
                                < btmodel[i-a*bi, j-a*bj, k-a*bk]:

                                grid_data[k-a*bk,j-a*bj,i-a*bi] = cmodel.settings.partition_number+1

                                if special_grid_data[k-a*bk, j-a*bj, i-a*bi] == 0:
                                    special_grid_data[k-a*bk,j-a*bj,i-a*bi] =\
                                        SPECIALpartition_number + hl*7

                                a += 1

        print(layers.A[l]["name"] + "\t| " + str(cmodel.settings.partition_number + 1))
        SPECIALpartition_number += 1 + hl*hl_offset
        cmodel.settings.partition_number += 1

    return grid_data, special_grid_data


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
        if cmodel.settings.s3d is False:
            triangles = Delaunay(points,qhull_options="QJ Pp Qw")
        # 3d surface
        else:
            triangles = Delaunay(points)

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
            filename, x, y, z,
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
