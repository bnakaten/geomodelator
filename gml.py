'''
GEOMODELATOR library file
'''

import csv
import json
import os
import re
import sys
import glob
import math
import numpy as np

from scipy.interpolate import griddata

from scipy.interpolate import interpn

from osgeo import ogr, gdal

import pyvista as pv

import logging

import yaml

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
#logging.basicConfig(stream=sys.stderr, level=logging.WARNING)


class Configuration:
    def __init__(self):
        self.modelGrid()
        self.dataLocation()
        self.csvFormat()
        self.partitionDetails()
        self.feature()

    def modelGrid(
        self,
        modelTyp3d=True,
        cornerPoint1=[.0,.0,.0],
        cornerPoint2=[.2,.0,.0],
        modelDimension=[.5,.0,.0],
        discretisationX=None,
        discretisationY=None,
        discretisationZ=None,
        cellNumberX=1,
        cellNumberY=1,
        cellNumberZ=1,
    ):
        self.modelType3d = modelTyp3d
        self.cornerPoint1 = cornerPoint1
        self.cornerPoint2 = cornerPoint2
        self.modelDimension = modelDimension
        self.discretisationX = discretisationX
        self.discretisationY = discretisationY
        self.discretisationZ = discretisationZ
        self.cellNumberX = cellNumberX
        self.cellNumberY = cellNumberY
        self.cellNumberZ = cellNumberZ

    def dataLocation(
        self,
        inputPath="",
        outputPath="",
        shapeFilePath="",
        modelFileVtk='model.vts',
    ):
        self.inputPath = inputPath
        self.outputPath = outputPath
        self.shapeFilePath = shapeFilePath
        self.modelFileVtk = modelFileVtk

    def csvFormat(
        self,
        columns=['x','y','z'],
        header=1,
        delimiter=',',
    ):
        self.csvColumns = columns
        self.csvHeader = int(header)
        self.csvDelimiter = delimiter

    def partitionDetails(
        self,
        structureWidth={},
        structureRank={},
        partitionCounter=0,
    ):
        self.structureWidth = structureWidth
        self.structureRank = structureRank
        self.partitionCounter = partitionCounter

    def feature(
        self,
        structureFileList=[],
        structureTypeList=[],
        structureExtensionList=[],
        maskFileName='mask.tiff',
        part={},
        rotationAngle=.0,
        modelExtensionWith=15.,
    ):
        self.maskFileName = maskFileName
        self.partitionList = ParitionList()
        self.structureFileList = structureFileList
        self.structureExtensionList = structureExtensionList

        self.structureTypeList = structureTypeList
        self.part = part
        self.rotationAngle = rotationAngle
        self.modelExtensionWith = modelExtensionWith


class Model:
    def __init__(
        self,
        configuration: Configuration,
        ):
        self.configuration = configuration

        if self.configuration.rotationAngle == 0:
            self.configuration.rotationAngle = Helper.calculate_full_angle(
                self.configuration.cornerPoint2,
                self.configuration.cornerPoint1,
            )


    def generateGridPoints(
            self
        ):
        x = np.zeros(len(self.configuration.discretisationX)+1)
        y = np.zeros(len(self.configuration.discretisationY)+1)
        z = np.zeros(len(self.configuration.discretisationZ)+1)

        x[0] = 0

        for i, value in enumerate(self.configuration.discretisationX):
            x[i+1] = value + x[i]

        self.x = x

        y[0] = 0

        for i, value in enumerate(self.configuration.discretisationY):
            y[i+1] = value + y[i]

        self.y = y

        z[0] = 0

        for i, value in enumerate(self.configuration.discretisationZ):
            z[i+1] = value + z[i]

        self.z = z

        logging.info("====================================")
        logging.info("\tModel")
        logging.info("====================================")
        logging.info("                          nx                   ")
        logging.info("                  P8                P7         ")
        logging.info("                 +-----------------+           ")
        logging.info("           ny   /|                /|           ")
        logging.info("               / |               / |           ")
        logging.info("           P5 /  |           P6 /  |   nz      ")
        logging.info("             +-----------------+   |           ")
        logging.info("             |   |             |   |           ")
        logging.info("             |   |             |   |           ")
        logging.info("             |   |P4           |   |P3         ")
        logging.info("             |   +-------------|---+           ")
        logging.info("             |  /              |  /            ")
        logging.info("             | /               | /             ")
        logging.info("             |/                |/              ")
        logging.info("             +-----------------+               ")
        logging.info("           P1                 P2               ")

        logMessage = f"({self.configuration.cornerPoint1[0]}," + \
            f"{self.configuration.cornerPoint1[1]}," + \
            f"{self.configuration.cornerPoint1[2]})"
        logging.info("Origin coordinate (= P1 original) (x, y, z): " + logMessage)

        logMessage = f"({len(self.configuration.discretisationX)}, "+ \
            f"{len(self.configuration.discretisationY)}, " + \
            f"{len(self.configuration.discretisationZ)})"
        logging.info("Discretisation (nx, ny, nz): " + logMessage)

        logging.info("Normalized model corner coordinates (x, y, z):")
        logging.info(f"\tP1 ({self.x[0]}, {self.y[0]}, {self.z[0]})")
        logging.info(f"\tP2 ({self.x[-1]}, {self.y[0]}, {self.z[0]})")
        logging.info(f"\tP3 ({self.x[-1]}, {self.y[-1]}, {self.z[0]})")
        logging.info(f"\tP4 ({self.x[0]}, {self.y[-1]}, {self.z[0]})")
        logging.info(f"\tP5 ({self.x[0]}, {self.y[0]}, {self.z[-1]})")
        logging.info(f"\tP6 ({self.x[-1]}, {self.y[0]}, {self.z[-1]})")
        logging.info(f"\tP7 ({self.x[-1]}, {self.y[-1]}, {self.z[-1]})")
        logging.info(f"\tP8 ({self.x[0]}, {self.y[-1]}, {self.z[-1]})")

        numberOfCells = len(self.configuration.discretisationX)*\
            len(self.configuration.discretisationY)*\
            len(self.configuration.discretisationZ)
        logging.info(f"\tNumber off cells/elements: {numberOfCells}")


        logging.info("Rotation angle: " + \
            f"{np.degrees(self.configuration.rotationAngle).round()}")

        self.complete()


    def complete(self):
        self.nx = len(self.x)
        self.ny = len(self.y)
        self.nz = len(self.z)


class CenterPointModel(Model):
    def generateGridPoints(
            self,
            cornerPointModel: Model,
        ):
        self.x = np.zeros(len(cornerPointModel.x)-1)
        self.y = np.zeros(len(cornerPointModel.y)-1)
        self.z = np.zeros(len(cornerPointModel.z)-1)

        self.x[:] = (cornerPointModel.x[1:] + cornerPointModel.x[:-1])/2
        self.y[:] = (cornerPointModel.y[1:] + cornerPointModel.y[:-1])/2
        self.z[:] = (cornerPointModel.z[1:] + cornerPointModel.z[:-1])/2

        self.complete()


class StructureList:
    def __init__(
        self,
        a = None,
        orientation = None,
        option = None,
        type = 'layer'
        ):
        self.A = a
        self.orientation = orientation
        self.option = option
        self.type = type
        self.id = 0


class ParitionList:
    def __init__(
        self,
        id=None,
        name=None,
        modelA=None,
        modelB=None,
        ):
        self.id = id
        self.name = name
        self.modelA = modelA
        self.modelB = modelB


class IndexStepControl:
    '''
        Data structure to hold the moving vector for the parition check
    '''
    def __init__(
            self,
            vectorDirection,
            vectorI,
            vectorJ,
            vectorK
        ):
        ## vector direction
        self.vd = vectorDirection
        ## vector i switch
        self.vi = vectorI
        ## vector j switch
        self.vj = vectorJ
        ## vector k switch
        self.vk = vectorK


class Helper:

    # load_and_organize_point_cloud_data
    @staticmethod
    def loadAndOrganizeStructureFileData(
        model: CenterPointModel,
        filename: str
        ):
        '''
            Load point cloud from text file
        '''
        data = None

        extensions = ['.csv', '.asc', '.txt', '.tmp', '.xyz']

        if filename.lower().endswith(tuple(extensions)):
            columns = []

            if filename.split(".")[-1] == 'asc':
                model.configuration.csvHeader = 1
                model.configuration.csvDelimiter = " "

                with open(filename, "rt") as f:
                    csvReader = csv.reader(f)
                    header = next(csvReader)
                    model.configuration.csvColumns = header[-1].split(
                        model.configuration.csvDelimiter
                    )[1:]


            if 'x' in model.configuration.csvColumns:
                x_id = list(model.configuration.csvColumns).index('x')
                columns.append(x_id)

            if 'y' in model.configuration.csvColumns:
                y_id = list(model.configuration.csvColumns).index('y')
                columns.append(y_id)

            if 'z' in model.configuration.csvColumns:
                z_id = list(model.configuration.csvColumns).index('z')
                columns.append(z_id)

            try:
                logging.info(f"-> {filename}")
                data = np.genfromtxt(
                    filename,
                    skip_header = model.configuration.csvHeader,
                    delimiter = model.configuration.csvDelimiter,
                    usecols = tuple(columns)
                )

            except ValueError:
                logging.info(
                    'Please check if the csv file column number and the csv \
                    column oder array fit together.'
                )
                raise

        return data


    @staticmethod
    def loadPointCloudFiles(filename):
        '''
            Load point cloud from text file
        '''

        data = None
        if filename.lower().endswith('.tmp'):
            data = np.genfromtxt(filename, skip_header=1,  delimiter=',')
        elif filename.lower().endswith('.asc'):
            data = np.genfromtxt(filename, skip_header=1,  delimiter=' ')

        return data


    @staticmethod
    def calculate_full_angle(P, origin):
        '''
            Calculate full angle (360 degrees) between two points
        '''
        return math.radians(
            math.atan2(P[1] - origin[1],
            P[0] - origin[0]) * (180.0 / math.pi)
        )


    @staticmethod
    def rotateCoordianteAroundOrigin(model, P, clockwise=-3):
        '''
            Rotate all points in list by rotation angle in origin
        '''
        # rotate x and y coordinates clockwise
        x = (P[0] * math.cos(clockwise*model.configuration.rotationAngle)) +\
                (P[1] * math.sin(clockwise*model.configuration.rotationAngle))
        y = -(P[0] * math.sin(clockwise*model.configuration.rotationAngle)) +\
                (P[1] * math.cos(clockwise*model.configuration.rotationAngle))
        P = [x, y, P[2]]

        return P


    @staticmethod
    def shiftCoordinateToOrigin(P, origin):
        '''
            Shift list of points into origin
        '''
        return [P[i] - origin[i] for i in range(3)]


    @staticmethod
    def interpolatePlane(
        xyGrid,
        vectorX,
        vectorY,
        indexing='xy',
        method=2,
        ):
        '''
            Method to interpolate a surface to plane (e.g. xy or yz or xz).
        '''
        interp_method = ['nearest', 'linear', 'cubic']

        coords = (xyGrid[:,0], xyGrid[:,1])
        values = xyGrid[:,2]

        x_grid, y_grid = np.meshgrid(vectorX, vectorY, indexing=indexing)


        try:
            grid = griddata(
                coords,
                values,
                (x_grid, y_grid),
                method=interp_method[method]
            )

            if np.isnan(grid).all():
                grid = None

        except ValueError:
            grid = None

        return grid


    def rotateModel(
        self,
        model: Model,
        a=None,
        b=None,
        c=None,
        indexing='ij',
        ):
        '''
            Rotate model by rotation angle
        '''

        if a is None:
            a = model.x

        if b is None:
            b = model.y

        if c is None:
            c = model.z

        z_matrix, y_matrix, x_matrix = np.meshgrid(
            c,
            b,
            a,
            indexing=indexing
        )

        z_list = [item for layer in z_matrix for row in layer for item in row]
        y_list = [item for layer in y_matrix for row in layer for item in row]
        x_list = [item for layer in x_matrix for row in layer for item in row]

        coord_list = np.array(
            [
                [x_list[i],y_list[i],z_list[i]] for i in range(len(x_list))
            ]
        )

        if model.configuration.rotationAngle != 0:
            coord_list = np.array(
                [
                    self.rotateCoordianteAroundOrigin(
                        model, coord_list[i,:], clockwise=-1
                    ) for i in range(len(coord_list))
                ]
            )

        if model.configuration.cornerPoint1 != [0,0,0]:
            coord_list =  np.array(
                [
                    self.shiftCoordinateToOrigin(
                        coord_list[i,:],
                        [-value for value in model.configuration.cornerPoint1]
                    ) for i in range(len(coord_list))
                ]
            )

        vtk_grid_x = np.array([coord[0] for coord in coord_list]).reshape((\
                                len(c),len(b),len(a)))
        vtk_grid_y = np.array([coord[1] for coord in coord_list]).reshape((\
                                len(c),len(b),len(a)))
        vtk_grid_z = np.array([coord[2] for coord in coord_list]).reshape((\
                                len(c),len(b),len(a)))

        return vtk_grid_x, vtk_grid_y, vtk_grid_z


    @staticmethod
    def normalizeGrid(self, model):
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


    @staticmethod
    def convert2dTo3dData(model, points):
        '''
            Convert 2d csv data to 2d+/3d data
        '''
        ## convert 2d data to interpolatable 3d data
        third_column = np.full(
            len(points[:,0]),
            -(model.configuration.modelDimension[1]/2) \
                -model.configuration.modelExtensionWith
        )
        missing_column = 0

        if 'x' not in model.configuration.csvColumns:
            points = np.column_stack((third_column,points))
        elif 'y' not in model.configuration.csvColumns:
            column_x = points[:,0]
            column_z = points[:,1]
            points = np.column_stack((column_x,third_column))
            points = np.column_stack((points,column_z))
            missing_column = 1
        elif 'z' not in model.configuration.csvColumns:
            points = np.column_stack((points,third_column))
            missing_column = 2

        # add to data a copy of the existing data with shifted 3rd coordinate to 100
        tmp_points = np.copy(points)
        tmp_points[:,missing_column] += (
            model.configuration.modelDimension[1] +model.configuration.modelExtensionWith
        )
        return np.vstack((points,tmp_points))


    # @staticmethod
    def pointsToPyvistaSurface(self, model, filename, x, y, z):
        '''
            Create vtk output of layers, faults or seams
        '''

        status = True
        points = np.zeros( (len(x), 3) )
        for i, x_value in enumerate(x):
            points[i,0] = x_value
            points[i,1] = y[i]
            try:
                points[i,2] = z[i]
            except:
                points[i,2] = None
                pass

        points = points[~np.isnan(points).any(axis=1),:]

        if model.configuration.rotationAngle != 0:
            points = np.array(
                [
                    self.rotateCoordianteAroundOrigin(model, points[i,:], clockwise=-1)
                    for i in range(len(points))
                ]
            )

        if model.configuration.cornerPoint1 != [0,0,0]:
            points =  np.array(
                [
                    self.shiftCoordinateToOrigin(
                        points[i,:],
                        [-value for value in model.configuration.cornerPoint1]
                    )
                    for i in range(len(points))
                ]
            )

        try:
            points = points[~np.isnan(points).any(axis=1),:]
            vtk_filename = filename + ".vtk"
            pcl = pv.PolyData(points)
            mesh = pcl.delaunay_2d(tol=1e-02)
            #mesh = pcl.delaunay_3d().extract_surface()
            mesh.save(vtk_filename)
        except:
            pass


        return status


    @staticmethod
    def pyConfToJson(confFile):
        variables = {}

        # Regular expression to match variable assignments
        pattern = r'(\w+)\s*=\s*(.+)'

        with open(confFile, 'r') as file:
            for line in file:
                # Match variable assignments
                match = re.match(pattern, line.strip())
                if match:
                    variable_name, value = match.groups()
                    # Try to evaluate the value and store it in the dictionary
                    try:
                        value = eval(value)
                        variables[variable_name] = value
                    except:
                        # If evaluation fails, store the value as a string
                        variables[variable_name] = value

        configuration = {}
        for key, value in variables.items():
            if key == 'CALLDIR':
                configuration['config'] = {}
                configuration['config']['base'] = value
                configuration['config']['shapePath'] = ''

            elif key == 'S3D':
                configuration['model'] = {}
                configuration['model']['3d'] = value

            elif key == 'IPATH':
                configuration['config']['inputPath'] = value

            elif key == 'OPATH':
                configuration['config']['outputPath'] = value

            elif key == 'SPATH':
                configuration['config']['shapePath'] = value

            elif key == 'CSVCOLUMNS':
                configuration['config']['csv'] = {}
                configuration['config']['csv']['columns'] = value

            elif key == 'CSVDELIMITER':
                configuration['config']['csv']['delimiter'] = value

            elif key == 'CSVHEADER':
                configuration['config']['csv']['header'] = value

            elif key == 'MODELP1':
                configuration['model']['cornerPoint1'] = value

            elif key == 'MODELP2':
                configuration['model']['cornerPoint2'] = value

            elif key == 'MODELDIMENSION':
                configuration['model']['dimension'] = value

            elif key == 'nx':
                configuration['model']['cellNumberX'] = value
                configuration['model']['discretizationType'] = "n"

            elif key == 'ny':
                configuration['model']['cellNumberY'] = value

            elif key == 'nz':
                configuration['model']['cellNumberZ'] = value

            elif key == 'dx':
                print("The conversion of dx is not possible automaticaly! Please enter the real dx values as numerical values in the configuration file and name them \"DX\". e.g. \"DX = [10, 20, 10, 5, 5]\". The original variable dx must be removed or commented out.")
                sys.exit(1)
            elif key == 'DX':
                configuration['model']['discretizationX'] = value
                configuration['model']['discretizationType'] = "d"

            elif key == 'dy':
                print("The conversion of dy is not possible automaticaly! Please enter the real dx values as numerical values in the configuration file and name them \"DY\". e.g. \"DY = [10, 20, 10, 5, 5]\". The original variable dy must be removed or commented out.")
                sys.exit(1)
            elif key == 'DY':
                configuration['model']['discretizationY'] = value

            elif key == 'dz':
                print("The conversion of dz is not possible automaticaly! Please enter the real dx values as numerical values in the configuration file and name them \"DZ\". e.g. \"DZ = [10, 20, 10, 5, 5]\".  The original variable dz must be removed or commented out.")
                sys.exit(1)
            elif key == 'DZ':
                configuration['model']['discretizationZ'] = value

        fileList = Helper.loadFileList(
            configuration['config']['base'] + configuration['config']['inputPath']
        )

        configuration['structure'] = {}
        configuration['structure']['file'] = fileList


        return json.dumps(configuration, indent=4);


    @staticmethod
    def loadFileList(directory):
        fileList = []
        for filename in os.listdir(directory):
            if os.path.isfile(os.path.join(directory, filename)):
                fileList.append(filename)
        return fileList


    @staticmethod
    def yamlToJson(yamlFile):
        with open(yamlFile, 'r') as yamlFile:
            yamlData = yaml.safe_load(yamlFile)

        return json.dumps(yamlData, indent=4)

    @staticmethod
    def loadJson(jsonFile):
        with open(jsonFile, 'r') as jsonFile:
            jsonData = json.load(jsonFile)

        return json.dumps(jsonData, indent=4)


class Partitioning(Helper):

    def __init__(
        self,
        model: Model
    ):
        self.model = model


    def partitionateModel(
        self,
        model: CenterPointModel,
        structureType='layer',
        vtk=True,
        filenameList=None,
        filenameType=None
    ):
        '''
            Load text files containing layer, fault or seam point cloud data
        '''
        layers = StructureList(type=structureType)
        layers.A = {}
        layers.orientation = {}
        layers.fileType = {}
        layers.id = 0
        i = 0

        filenameMappingList = {}

        all_elevs = None

        logging.info(f"Process point {structureType} cloud files ...")


        if filenameList:
            filenameList = [(index, item) for index, item in enumerate(filenameList)]
            filenameList = dict(filenameList)
        else:
            return layers, all_elevs

        for key, filename in filenameList.items():

            if key >= 0 and filenameType[key] != structureType:
                continue

            logging.info(model.configuration.inputPath + filename)
            filename = os.path.basename(filename)

            if not os.path.exists(model.configuration.inputPath + filename):
                continue

            #### load xyz raster file
            points = self.loadAndOrganizeStructureFileData(
                    model,
                    model.configuration.inputPath + filename
            )

            # for 2d csv data
            if str(0) in model.configuration.csvColumns:
                points = self.convert2dTo3dData(model, points)

            # shift csv data to origin
            if model.configuration.cornerPoint1 != [0,0,0]:
                points =  np.array(
                    [
                        self.shiftCoordinateToOrigin(
                            points[i,:],
                            model.configuration.cornerPoint1
                        )
                        for i in range(len(points))
                    ]
                )

            # rotate csv data so that the data are parallel to x axis
            if model.configuration.rotationAngle != 0:
                points = np.array(
                    [
                        self.rotateCoordianteAroundOrigin(
                            model,
                            points[i,:],
                            clockwise=1
                        )
                        for i in range(len(points))
                    ]
                )

            # for 2d+ model extend 2d csv, from a point line to a surface
            if not model.configuration.modelType3d:
                points[:,1] += model.configuration.modelDimension[1]/2

                tmp_points_left = np.copy(points)
                left_shift = -model.configuration.modelDimension[1]/2 \
                    -model.configuration.modelExtensionWith
                tmp_points_left[:,1] += left_shift

                tmp_points_right = np.copy(points)
                right_shift = model.configuration.modelDimension[1]/2 \
                    +model.configuration.modelExtensionWith
                tmp_points_right[:,1] += right_shift
                points = np.vstack((tmp_points_left, tmp_points_right))

            if  "seam" in filename:
                tmp_points_left = np.copy(points)
                left_shift_x = -model.configuration.modelDimension[1]/(
                    model.configuration.cellNumberY*2
                )
                tmp_points_left[:,0] += left_shift_x
                left_shift_y = -model.configuration.modelDimension[1]/(
                    model.configuration.cellNumberY*2
                )
                tmp_points_left[:,1] += left_shift_y

                tmp_points_right = np.copy(points)
                right_shift_x = \
                    model.configuration.modelDimension[0]/model.configuration.cellNumberX
                tmp_points_right[:,0] += right_shift_x
                right_shift_y = \
                    model.configuration.modelDimension[1]/model.configuration.cellNumberY
                tmp_points_right[:,1] += right_shift_y
                points = np.vstack((tmp_points_left, tmp_points_right))


            # if structure_in_model_boundaries(model, points):
            np.savetxt(
                model.configuration.inputPath + os.path.splitext(
                    filename
                )[0] + ".tmp",
                points,
                fmt='%.2f',
                header='x,y,z',
                comments="",
                delimiter=","
            )

            filenameMappingList[os.path.splitext(filename)[0] + '.tmp'] = filename

        if filenameMappingList:
            logging.info(f"Best plane for interpolation of {structureType}" + \
                " values to the grid:")


        for filename in filenameMappingList:
            glob.glob(filename)
            filename = os.path.basename(filename)
            #### load xyz raster file
            points = self.loadPointCloudFiles(
                model.configuration.inputPath + filename
            )

            ## set NAN to 0
            points[np.isnan(points)] = 0

            ex = np.array(model.x)
            ey = np.array(model.y)

            if not model.configuration.modelType3d:
                eyLength = ey.size
                if(eyLength == 1):
                    eyLength = 2

                ey = np.linspace(
                    -model.configuration.modelExtensionWith,
                    model.configuration.modelDimension[1] \
                        +model.configuration.modelExtensionWith,
                    eyLength
                )

            ez = np.array(model.z)

            interpolation_quality = 0
            layers.orientation[i]= "xy"
            elevationXY = self.interpolatePlane(points, ex, ey)
            elevationXY = self.fixInterpolationVariationInZDirection(model, elevationXY)

            elevationXZ = None
            if model.configuration.modelType3d:
                newPoints = np.array([points[:,0], points[:,2], points[:,1]]).T
                elevationXZ = self.interpolatePlane(
                    newPoints,
                    ex, ez
                )

            newPoints = np.array([points[:,1], points[:,2], points[:,0]]).T
            elevationYZ = self.interpolatePlane(newPoints, ey, ez)
            elevationYZ = self.fixInterpolationVariationInXDirection(model, elevationYZ)

            all_elevs = [elevationXY, elevationXZ, elevationYZ]


            if elevationXY is None:
                layers.orientation[i]= "yz"

                if elevationYZ is None:
                    layers.orientation[i]= "xz"

                    if model.ny > 1 and elevationXZ is None:
                        logging.info(' coordinates of ' + structureType +\
                                ' could not be interpolated and will be ignored!\n')
                        continue

            elevationXY_quality = 0
            elevationXZ_quality = 0
            elevationYZ_quality = 0

            if elevationXY is not None:
                elevationXY_quality = np.sum(~np.isnan(elevationXY))/\
                    (model.nx*model.ny)

            if model.ny > 1 and elevationXZ is not None and structureType != "layer":
                elevationXZ_quality = np.sum(~np.isnan(elevationXZ))/\
                    (model.nx*model.nz)

            if elevationYZ is not None and structureType != "layer":
                elevationYZ_quality = np.sum(~np.isnan(elevationYZ))/\
                    (model.ny*model.nz)

            if interpolation_quality < elevationXY_quality:
                layers.orientation[i] = "xy"
                interpolation_quality = elevationXY_quality

            if interpolation_quality < elevationXZ_quality and structureType != "layer":
                layers.orientation[i] = "xz"
                interpolation_quality = elevationXZ_quality

            if interpolation_quality < elevationYZ_quality and structureType != "layer":
                layers.orientation[i] = "yz"
                interpolation_quality = elevationYZ_quality

            width = 0
            if filename.lower().split(".")[0] in model.configuration.structureWidth:
                width = model.configuration.structureWidth[filename.split(".")[0]]

            if filename.split(".")[0] in model.configuration.structureRank:
                layer_option = model.configuration.structureRank[filename.split(".")[0]]
            else:
                layer_option = 0


            positionInFileList = model.configuration.structureFileList.index(
                filenameMappingList[filename]
            )
            layers.fileType[i] =  model.configuration.structureTypeList[
                positionInFileList
            ]

            if layers.orientation[i] == "xy":
                layers.A[i] = {
                    "name" : filename.split(".")[0],
                    "width": width/2,
                    "elevation": elevationXY,
                    "option" : layer_option,
                    "fileType" : model.configuration.structureTypeList[
                        positionInFileList
                    ],
                }
                layers.id = 1

            if layers.orientation[i] == "xz":
                layers.A[i] = {
                    "name" : filename.split(".")[0],
                    "width": width/2,
                    "elevation": elevationXZ,
                    "option" : layer_option,
                    "fileType" : model.configuration.structureTypeList[positionInFileList],
                }
                layers.id = 2

            if layers.orientation[i] == "yz":
                layers.A[i] = {
                    "name" : filename.split(".")[0],
                    "width": width/2,
                    "elevation": elevationYZ,
                    "option" : layer_option,
                    "fileType" : model.configuration.structureTypeList[positionInFileList],
                }
                layers.id = 3

            logging.info(
                f"\tinterpolation plane {layers.orientation[i]}\t-> {layers.id} \t"
            )

            if vtk:
                xmesh, ymesh = np.meshgrid(ex, ey)
                zmesh = elevationXY

                vtk_filename = model.configuration.outputPath + \
                    os.path.splitext(filename)[0]

                status = self.pointsToPyvistaSurface(
                    model,
                    vtk_filename,
                    xmesh.flatten(),
                    ymesh.flatten(),
                    zmesh.flatten()
                )

                if status:
                    logging.info(f"\t\t{vtk_filename}.vtk")
                else:
                    logging.warning(
                        "\tWARNING: Interpolation of surface not successful." +\
                        f"\tBuilding of {vtk_filename}.vtu not possible!"
                    )

            os.remove(model.configuration.inputPath + \
                os.path.splitext(filename)[0] + '.tmp')

            i += 1

        return layers, all_elevs


    def fixInterpolationVariationInZDirection(self, model, elevationMatrix):
        if not model.configuration.modelType3d and elevationMatrix is not None:

            elevationMatrixRowMean = np.nanmean(elevationMatrix, axis=0)
            elevationMatrix = np.tile(
                elevationMatrixRowMean, (elevationMatrix.shape[0], 1)
            )
        return elevationMatrix


    def fixInterpolationVariationInXDirection(self, model, elevationMatrix):
        if not model.configuration.modelType3d and elevationMatrix is not None:
            elevationMatrixRowMean = np.nanmean(elevationMatrix, axis=1)
            elevationMatrix = np.tile(
                elevationMatrixRowMean, (elevationMatrix.shape[1], 1)
            )
        return elevationMatrix




    def setNewParitionId(
        self,
        layers,
        l,
        lx,
        ly,
        current_cell_center_elevation,
        pt
    ):
        '''
            Set new parition ids
        '''
        new_pid = False

        current_layer_elevation = layers.A[l]['elevation'][ly, lx]

        # ist die das mittelpunkt des aktullen elementes unterhalb der aktuellen
        # horizontoberfläche (layer surface) dann gehört das element zur parition
        # des neuen layers
        if current_cell_center_elevation < current_layer_elevation:
            current_layer_priority = layers.A[l]['option']
            new_pid = l == 0

            # new_pid = True

            # falls aber bereits eine andere partition existiert die eine höhere
            # priorität (rank) hat dann eben nicht, also wir hier noch mal gegen alle
            # anderen horizontoberflächen die priorität geprüft
            for pi in (range(l)):
                try:
                    old_layer_elevation = layers.A[pi]['elevation'][ly, lx]
                    old_layer_priority = pt.mode1[pi]
                except:
                    pass

                if old_layer_priority == current_layer_priority:
                    new_pid = True
                # an already processed layer surface with higher priority cuts the
                # current layer surface
                elif old_layer_priority > current_layer_priority and \
                    old_layer_elevation > current_cell_center_elevation:
                    new_pid = True

                # if old_layer_elevation >= current_layer_elevation and \
                #    old_layer_priority >= current_layer_priority:
                #     new_pid = True

                # if current_layer_elevation < old_layer_elevation and \
                #     current_layer_priority == old_layer_priority:
                #     new_pid = True

                # if current_cell_center_elevation < old_layer_elevation < \
                #     current_layer_elevation and current_layer_priority < \
                #     old_layer_priority:
                #     new_pid = True

                # if current_layer_elevation < old_layer_elevation and \
                #     current_layer_priority < old_layer_priority:
                #     new_pid = True

                # if old_layer_elevation < current_cell_center_elevation < \
                #     current_layer_elevation and current_layer_priority < \
                #     old_layer_priority:
                #     new_pid = False

        return new_pid


    def generateModelCellIds(self, model):
        '''
            Generate model cell ids for model post processing
        '''
        grid_data = np.zeros((model.nz, model.ny, model.nx), np.int32)
        grid_i_ids = np.zeros_like(grid_data)
        grid_j_ids = np.zeros_like(grid_data)
        grid_k_ids = np.zeros_like(grid_data)

        for k in range(0, model.nz):
            for j in range(0, model.ny):
                for i in range(0, model.nx):
                    grid_data[k, j, i] = i + j*model.nx + k*model.nx*model.ny
                    grid_i_ids[k, j, i] = i
                    grid_j_ids[k, j, i] = j
                    grid_k_ids[k, j, i] = k

        return grid_data, grid_i_ids, grid_j_ids, grid_k_ids


    def generateActivePartition(
        self,
        model: CenterPointModel,
        old_grid_data=None,
        name="active",
        filenameList=None,
        filenameType=None
        ):
        '''
            Generate model cell clusters regarding the shapefile mask
            and add a patiation number to this cells
        '''


        if filenameList:
            if not old_grid_data:
                old_grid_data = np.zeros((model.nz, model.ny, model.nx), np.int32)

            special_grid_data = np.zeros_like(old_grid_data)

            filenameList = [(index, item) for index, item in enumerate(filenameList)]
            filenameList = dict(filenameList)
        else:
            if not old_grid_data:
                old_grid_data = np.ones((model.nz, model.ny, model.nx), np.int32)

            special_grid_data = np.ones_like(old_grid_data)

            return old_grid_data, special_grid_data

        maskFilename = False
        for key, filename in filenameList.items():

            if key >= 0 and filenameType[key] == 'mask':
                maskFilename = model.configuration.inputPath + filename
                maskFileKey = key
                break

        if maskFilename:

            if maskFilename.split('.')[-1] == 'asc':

                logging.info("====================================")
                logging.info("\tGenerate " + name + " partitions")
                logging.info("====================================")
                logging.info("-> " + maskFilename)

                filename = os.path.basename(maskFilename)
                #### load xyz raster file
                points = self.loadPointCloudFiles(
                    model.configuration.inputPath + filename
                )


                ## set NAN to 0
                points[np.isnan(points)] = -1

                if len(model.x) < 2 and len(model.y) < 2:
                    x = model.x[0]
                    ex = np.array([x-abs(x), x, x+abs(x)])
                    y = model.y[0]
                    ey = np.array([y-abs(y), y, y+abs(y)])
                else:
                    ex = np.array(model.x)
                    ey = np.array(model.y)

                elevationXY = self.interpolatePlane(points, ex, ey)

                if len(model.x) < 2 and len(model.y) < 2:
                    if elevationXY != None and elevationXY[1,1] == 0:
                        elevationXY = np.array([elevationXY[1,1]])
                    else:
                        elevationXY = None


                if elevationXY.any() != None:
                    elevationXY_ = np.flipud(elevationXY)
                    elevationXY_[np.isnan(elevationXY_)] = -1

                    old_grid_data[:,elevationXY_[::-1]==0] = 1
                    old_grid_data[:,elevationXY_[::-1]==-1] = 0
                    special_grid_data[:,elevationXY_[::-1]==0] = 1
                    special_grid_data[:,elevationXY_[::-1]==-1] = 0
                else:
                    old_grid_data[:,:] = 1
                    special_grid_data[:,:] = 1

                pointCloud = pv.PolyData(points)

                vtk_filename = model.configuration.outputPath + \
                    filename.split('.')[0] + '.vtk'
                pcl = pointCloud.delaunay_3d(alpha=3.0)
                surface = pcl.extract_geometry()
                surface.save(vtk_filename)

            else:

                tiff_file = \
                    model.configuration.outputPath + model.configuration.maskFileName

                if not old_grid_data:
                    old_grid_data = np.zeros((model.nz, model.ny, model.nx), np.int32)

                logging.info("====================================")
                logging.info("\tGenerate " + name + " partitions")
                logging.info("====================================")
                logging.info("-> " + maskFilename)

                shape = ogr.Open(maskFilename)

                shape_layer = shape.GetLayer()
                x_min, x_max, y_min, y_max = shape_layer.GetExtent()

                x_pixel_size = int(np.min(model.configuration.discretisationX)/2)
                y_pixel_size = int(np.min(model.configuration.discretisationY)/2)

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

                x_vector_model = model.x[:] + x_min
                y_vector_model = model.y[:] + y_min
                model_mask = Helper.interpolatePlane(
                    mask_data.T, x_vector_model, y_vector_model[::-1]
                )

                old_grid_data[:,model_mask[::-1]>110] = 1
                special_grid_data[:,model_mask[::-1]>110] = 1

        else:
            if not old_grid_data.all():
                old_grid_data = np.ones((model.nz, model.ny, model.nx), np.int32)

            special_grid_data = np.ones_like(old_grid_data)

        return old_grid_data, special_grid_data


    def GenerateLayerPartitions(
        self,
        model,
        layers,
        old_grid_data
    ):
        '''
            Generate model cell clusters regarding layers
            and add a partition number to this cells
        '''
        grid_data = np.ones((model.nz, model.ny, model.nx), np.int32)

        if old_grid_data is None:
            old_grid_data = np.ones_like(grid_data)

        partition_data = np.ones_like(grid_data)
        special_grid_data = np.zeros_like(grid_data)

        if len(layers.A) == 0:
            grid_data = old_grid_data

        SPECIALpartitionCounter = 1
        tmodel = None

        logging.info("====================================")
        logging.info("\tGenerate structure partitions")
        logging.info("====================================")

        pt = model.configuration.partitionList
        pt.pname = []
        pt.mode1 = []
        pt.pid = []

        rankLayers = []

        model.configuration.part['all'] = 0
        model.configuration.part['default'] = 1

        logging.info("inaktiv  \t-> 0")
        model.configuration.partitionCounter += 1
        logging.info("default  \t-> 1")

        for l in range(0, len(layers.A)):
            model.configuration.partitionCounter += 1
            pt.pname.append(layers.A[l]['name'])
            pt.mode1.append(layers.A[l]['option'])

            if layers.A[l]['option'] > 0:
                rankLayers.append(
                    [model.configuration.partitionCounter, layers.A[l]['option']]
                )
                # print(rankLayers[-1][0],rankLayers[-1][1])
                print(layers.orientation[l])

            pt.pid.append(model.configuration.partitionCounter)
            model.configuration.part[layers.A[l]['name']] = \
                model.configuration.partitionCounter
            logging.info(f"{pt.pname[l]}\t-> {pt.pid[l]}")

            if layers.orientation[l] == "xy":
                ijk = 2
                indexSteps = IndexStepControl(2, 0, 0, 1)
            elif layers.orientation[l] == "xz":
                ijk = 1
                indexSteps = IndexStepControl(1, 0, 1, 0)
            else:
                ijk = 0
                indexSteps = IndexStepControl(0, 1, 0, 0)

            tmodel = np.meshgrid(model.x, model.y, model.z, indexing='ij')[ijk]


            if layers.A[l]['fileType'] == "layer":
                for i in range(0, model.nx):
                    for j in range(0, model.ny):
                        for k in range(0, model.nz):
                            if layers.orientation[l] == "xy":
                                lx = i
                                ly = j
                            elif layers.orientation[l] == "xz":
                                lx = i
                                ly = k
                            else:
                                lx = j
                                ly = k

                            new_pid = self.setNewParitionId(
                                layers, l, lx, ly, tmodel[i, j, k], pt
                            )

                            if old_grid_data[k,j,i] != 0:
                                if new_pid:
                                    grid_data[k, j, i] = pt.pid[l]
                            else:
                                grid_data[k, j, i] = 0


                            if partition_data[k,j,i] != 0:
                                if new_pid:
                                    partition_data[k, j, i] = pt.pid[l]
                            else:
                                partition_data[k, j, i] = 0


            elif layers.A[l]['fileType'] == "fault":
                grid_data, special_grid_data = self.identifyParitionElements(
                    model,
                    grid_data,
                    special_grid_data,
                    indexSteps,
                    layers.orientation[l],
                    layers.A[l],
                    layers,
                    SPECIALpartitionCounter,
                    rankLayers
                )

                SPECIALpartitionCounter += 1

        logging.info("Total number of partitions: " + str(model.configuration.partitionCounter+1))

        return grid_data, partition_data, special_grid_data


    def identifyParitionElements(
        self,
        model,
        grid_data,
        special_grid_data,
        a,
        layer_orientation,
        layer,
        layers,
        partitionCounter,
        rankLayers
    ):
        '''
            identify elements that are part of partition
        '''
        tmodel = np.meshgrid(model.x, model.y, model.z, indexing='ij')[a.vd]
        flag = 0
        flag2 = 0
        flag3 = 0

        for i in range(0, model.nx - a.vi):
            for j in range(0, model.ny - a.vj):
                for k in range(0, model.nz - a.vk):


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
                            layer_elevation = layer["elevation"][k, j]


                    nk = k + a.vk
                    nj = j + a.vj
                    ni = i + a.vi
                    # value = np.float64('nan')
                    
                    ## if element is inactive do nothing on this element and go to next element
                    # if  grid_data[k,j,i] == 0:
                    #     continue

                    # ## if the layer elevation value is not set for this element jump to next element
                    # if layer["elevation"] is None:
                    #     special_grid_data[k, j, i] = 0
                    #     continue
                    
                    ## ???? if all layer elevation values are not set for this element jump to next element
                    if np.isnan(np.float64(layer_elevation)):
                        continue


                    ## ???? compare current layer ranking with all other layer rankings
                    if np.isnan(np.float64(layer_elevation)) and any(
                        value[0] > grid_data[k,j,i] and \
                            layer["option"] < value[1] for value in rankLayers
                    ):
                        continue


                    # if flag == 0:
                    #     print(layer["width"])
                    #     flag = 1
                    
                    # if flag == 0:
                    #     idx, idxRank  = next(((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[k,j,i]), None)
                    #     print("> ",idx, idxRank)
                    #     flag = 1
                    
                    # if flag2 == 0:
                    #     idx, idxRank  = next(((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[k,j,i]), None)
                    #     if tmodel[i, j, k] <= layer_elevation < tmodel[ni,nj,nk] and idxRank < layer["option"]:
                    #         print(">> ",layer["option"] ,  idxRank)
                    #         flag2 = 1

                    


                    # if np.isnan(np.float64(layer_elevation)) or any(
                    #     value[0] > grid_data[k,j,i] and \
                    #         layer["option"] < value[1] for value in rankLayers
                    # ):
                    #     continue
                    # gridPID = grid_data[k,j,i]
                    # tmodelElevationCC = tmodel[i, j, k]
                    # tmodelElevationNC = tmodel[ni,nj,nk]
                    # if tmodel[i, j, k] <= layer_elevation < tmodel[ni,nj,nk] and \
                    #     layer["elevation"] is not None:
                    
                    ## if the current element elevation value is between 
                    ## the current model element elevation value and the next element elevation value
                    ## this indicates that the layer is between the tow elements
                    idx, idxRank = next(
                        ((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[k, j, i]),
                        (None, None)
                    )
                    if tmodel[i, j, k] <= layer_elevation < tmodel[ni,nj,nk]:
                        # and idxRank is not None and  idxRank < layer["option"]:
                        # (idxRank is None or  idxRank < layer["option"]):
                    #if tmodel[i, j, k] <= layer_elevation < tmodel[ni,nj,nk]:



                        # for value in rankLayers:
                        #     print("a", tmodel[i, j, k], layer_elevation, tmodel[ni,nj,nk]) 
                        #     print("b", value[0], grid_data[k,j,i], layer["option"], value[1]) 


                        if idxRank is not None and  idxRank < layer["option"]:
                            grid_data[k,j,i] = model.configuration.partitionCounter
                            grid_data[nk,nj,ni] = model.configuration.partitionCounter

                            if special_grid_data[k, j, i] == 0:
                                special_grid_data[k, j, i] = partitionCounter
                            if special_grid_data[nk,nj,ni] == 0:
                                special_grid_data[nk,nj,ni] = partitionCounter

                        t = 1
                        ptk = k - t*a.vk
                        ptj = j - t*a.vj
                        pti = i - t*a.vi


                        valid_offset = False
                        if 0 <= pti and 0 <= ptj and 0 <= ptk:
                            valid_offset = True

                            idx, idxRank = next(
                                ((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[ptk,ptj,pti]),
                                (None, None)
                            )

                        # print("> k: " + str(k) + " - " +str(t) + " * " +str(a.vk) + " = " +str(ptk) + " ; ", end="")
                        # print("  j: " + str(j) + " - " +str(t) + " * " +str(a.vj) + " = " +str(ptj) + " ; ", end="")
                        # print("  i: " + str(i) + " - " +str(t) + " * " +str(a.vi) + " = " +str(pti) + " ; ")                                
                        # print("* " + str(t) + " : " + str(valid_offset) + " : "  + str(layer_elevation - layer["width"]) + " : " + str(tmodel[pti,ptj,ptk]) + " : " + str(idxRank) + " : " +str(layer["option"]) + " - ")


                        while valid_offset and\
                            layer_elevation - layer["width"] < tmodel[pti,ptj,ptk]:
                            # and idxRank is not None and  idxRank < layer["option"]:
                            # (idxRank is None or  idxRank < layer["option"]):

                            if idxRank is not None and  idxRank < layer["option"]:
                                grid_data[ptk,ptj,pti] = \
                                    model.configuration.partitionCounter

                                if special_grid_data[ptk,ptj,pti] == 0:
                                    special_grid_data[ptk,ptj,pti] = partitionCounter
                            
                            t += 1
                            ptk = k - t*a.vk
                            ptj = j - t*a.vj
                            pti = i - t*a.vi

                            valid_offset = False
                            if 0 <= pti and 0 <= ptj and 0 <= ptk:                 
                                # print("= " + str(t) + " : "  + str(layer_elevation - layer["width"]) + " : " + str(tmodel[pti,ptj,ptk]) + " : " + str(idxRank) + " : " +str(layer["option"]) + " - ")

                                valid_offset = True

                                idx, idxRank = next(
                                    ((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[ptk,ptj,pti]),
                                    (None, None)
                                )


                        t = 1
                        ntk = k + t*a.vk
                        ntj = j + t*a.vj
                        nti = i + t*a.vi


                        valid_offset = False
                        if nti < model.nx and ntj < model.ny and ntk < model.nz:
                            valid_offset = True

                            idx, idxRank = next(
                                ((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[ntk,ntj,nti]),
                                (None, None)
                            )

                        while valid_offset and\
                            layer_elevation + layer["width"] > tmodel[nti,ntj,ntk]:
                            # and idxRank is not None and idxRank <= layer["option"]:

                            if idxRank is not None and  idxRank < layer["option"]:
                                grid_data[ntk,ntj,nti] = \
                                    model.configuration.partitionCounter

                                if special_grid_data[ntk,ntj,nti] == 0:
                                    special_grid_data[ntk,ntj,nti] = partitionCounter

                            t += 1
                            ntk = k + t*a.vk
                            ntj = j + t*a.vj
                            nti = i + t*a.vi


                            valid_offset = False
                            if nti < model.nx and ntj < model.ny and ntk < model.nz:
                                valid_offset = True

                            
                                idx, idxRank = next(
                                    ((sub[0], sub[1]) for sub in rankLayers if sub[0] == grid_data[ntk,ntj,nti]),
                                    (None, None)
                                )
                        # print(str(t) + " : " + str(layer_elevation + layer["width"]) + " : " + str(tmodel[nti,ntj,ntk]) + " : " + str(idxRank) + " : " +str(layer["option"]))

                        # print(str(i) + " ; " +str(t) + " ; " +str(a.vi) + " ; " +str(model.nx) + " ; ")
                    else:
                        special_grid_data[k, j, i] = 0

        return grid_data, special_grid_data