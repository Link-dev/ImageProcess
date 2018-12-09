# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 12:43:57 2018

@author: ZJX
"""

from osgeo import gdal
from osgeo.gdalconst import *
import numpy as np
import matplotlib.pyplot as plt

FileName = ''
Dataset = gdal.Open(FileName, GA_ReadOnly)

ImgXSize = Dataset.RasterXSize  
ImgYSize = Dataset.RasterYSize  
BandCount = Dataset.RasterCount

geotrans = Dataset.GetGeoTransform()  #Affine matrix
proj = Dataset.GetProjection() #projection

Data = []
for i in range(0, BandCount):
    band = Dataset.GetRasterBand(i + 1)
    data = band.ReadAsArray(0, 0, ImgXSize, ImgYSize)
    Data.append(data)
    

Driver = gdal.GetDriverByName("GTiff")
OutFile = 'Save.tif'
dataSave = Driver.Create(OutFile, ImgXSize, ImgYSize, BandCount, GDT_Float64)
dataSave.SetGeoTransform(geotrans)              #
dataSave.SetProjection(proj)                    #

for i in range(0, BandCount):
    dataSave.GetRasterBand(i+1).WriteArray(Data[i])
    
    