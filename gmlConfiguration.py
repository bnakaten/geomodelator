import sys
import numpy as np
import json

import logging

import gml as gml

# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.basicConfig(stream=sys.stderr, level=logging.WARNING)

class gmlConfiguration:

    def __init__(self, runSettings):
        if type(runSettings) is dict:
            self.rs = runSettings
        else:
            self.rs = json.loads(runSettings)

        self.rs['run'] = False

    def configureModel(self, baseDir = None):
        if baseDir:
            self.rs['config']['base'] = baseDir

        if 'compress' not in self.rs['config']['csv']:
            self.rs['config']['csv']['compress'] = False

        def isNotNanOrNone(value):
            return value != "NaN" and value != "None"

        if self.rs['model']['discretisationType'] == 'd':
            discretisationX = [
                v for v in self.rs['model']['discretisationX'] if isNotNanOrNone(v)
            ]
            discretisationY = [
                v for v in self.rs['model']['discretisationY'] if isNotNanOrNone(v)
            ]
            discretisationZ = [
                v for v in self.rs['model']['discretisationZ'] if isNotNanOrNone(v)
            ]
        else:
            discretisationX = np.asarray(
                [self.rs['model']['dimension'][0]/self.rs['model']['cellNumberX']]\
                *self.rs['model']['cellNumberX']
            )
            discretisationY = np.asarray(
                [self.rs['model']['dimension'][1]/self.rs['model']['cellNumberY']]\
                *self.rs['model']['cellNumberY']
            )
            discretisationZ = np.asarray(
                [self.rs['model']['dimension'][2]/self.rs['model']['cellNumberZ']]\
                *self.rs['model']['cellNumberZ']
            )

        if 'extension' not in self.rs['structure'] or \
            self.rs['structure']['extension'] == None:
            self.rs['structure']['extension'] = []
            for file in self.rs['structure']['file']:
                self.rs['structure']['extension'].append(file.split('.')[-1])

        gml_gs = gml.Configuration()

        gml_gs.csvFormat(
            columns=self.rs['config']['csv']['columns'],
            header=int(self.rs['config']['csv']['header']),
            delimiter=self.rs['config']['csv']['delimiter'],
        )

        gml_gs.modelGrid(
            modelTyp3d=self.rs['model']['3d'],
            cornerPoint1=self.rs['model']['cornerPoint1'],
            cornerPoint2=self.rs['model']['cornerPoint2'],
            modelDimension=self.rs['model']['dimension'],
            discretisationX=discretisationX,
            discretisationY=discretisationY,
            discretisationZ=discretisationZ,
            cellNumberX=self.rs['model']['cellNumberX'],
            cellNumberY=self.rs['model']['cellNumberY'],
            cellNumberZ=self.rs['model']['cellNumberZ'],
        )

        if 'inputPath' not in self.rs['config']:
            self.rs['config']['inputPath'] = ''

        if 'outputPath' not in self.rs['config']:
            self.rs['config']['outputPath'] = ''

        if 'shapePath' not in self.rs['config']:
            self.rs['config']['shapePath'] = ''

        if 'modelFileVtk' not in self.rs['config']:
            self.rs['config']['modelFileVtk'] = 'model.vts'


        baseDir = \
            self.rs['config']['base'] + '/' if self.rs['config']['base'] != "" else ''
        inputDir = \
            self.rs['config']['inputPath'] + '/' \
                if self.rs['config']['inputPath']  != "" else ''
        outputDir = \
            self.rs['config']['outputPath'] + '/' \
                if self.rs['config']['outputPath'] != "" else ''
        shapeDir = \
            self.rs['config']['shapePath'] + '/' \
                if self.rs['config']['shapePath'] != "" else ''

        gml_gs.dataLocation(
            inputPath = baseDir + inputDir,
            outputPath = baseDir + outputDir,
            shapeFilePath = baseDir + shapeDir,
            modelFileVtk = baseDir + outputDir + self.rs['config']['modelFileVtk'],
        )

        if 'structureWidth' not in self.rs['model']:
            self.rs['model']['structureWidth'] = []

        if 'structureRank' not in self.rs['model']:
            self.rs['model']['structureRank'] = []

        gml_gs.partitionDetails(
            structureWidth=self.rs['model']['structureWidth'],
            structureRank=self.rs['model']['structureRank'],
        )

        if 'structure' in self.rs:
            gml_gs.feature(
               maskFileName='mask.tiff',
               structureFileList=self.rs['structure']['file'],
               structureExtensionList=self.rs['structure']['extension'],
               structureTypeList=self.rs['structure']['type'],
            )

            gml_gs.part =  {}
            gml_gs.partitionCounter = 0

        self.gml_rs = gml_gs