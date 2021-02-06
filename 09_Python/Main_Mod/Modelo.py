# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 06:24:14 2020

@author: Oscar
"""
#inspirado en Gidahatari fully geospatial model
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy as fp
from flopy.utils.gridgen import Gridgen
from flopy.utils.reference import SpatialReference
import  shapefile as sf #tomado de github flopy --pip install pyshp
from osgeo import gdal # tomado de gidahatari
import pandas as pd
from scipy.interpolate import griddata

import mplleaflet

#inspirado en grid intersection demo
import matplotlib as mpl
import flopy.discretization as fgrid
import flopy.plot as fplot
import shapely 
from shapely.geometry import Polygon, Point, LineString, MultiLineString, MultiPoint, MultiPolygon, shape#lastone form StackO.F.
from shapely.strtree import STRtree
from flopy.utils.gridintersect import GridIntersect



model_name= "modelo_Norte"
# Creation of workspace folder
workspace = os.path.join("data", model_name)
if not os.path.exists(workspace):
    os.makedirs(workspace)

#creating general simulation
    
sim =fp.mf6.MFSimulation(sim_name=model_name, version="mf6",
                         exe_name=r"C:\WRDAPP\mf6.2.0\bin\mf6",
                         sim_ws=workspace)

#setting modflow time
# [perlen, nstp, tsmult]
time_disc=[(1.0, 1, 1.0)]
tdis= fp.mf6.ModflowTdis(sim, pname="tdis",
                         time_units="SECONDS", 
                         nper=1, perioddata=time_disc)

# Create gwf model
model_nam_file=f"{model_name}.nam"
gwf = fp.mf6.ModflowGwf(sim, modelname=model_name,
                        model_nam_file=model_nam_file,
                        newtonoptions=None)


# Setting the solver
ims=  fp.mf6.modflow.mfims.ModflowIms(sim, pname="ims",
                                      complexity= "Simple")


# open shapefiles of limits and refinement
path_sh= r"D:\OneDrive - UNIVERSIDAD INDUSTRIAL DE SANTANDER\Maestria\06_Tesis\01_Tesis_Dev\02_Shp_Vect\Input_Model"
ModelLimitShp=sf.Reader(path_sh+"/Domin_Mod.shp" )
GalleryShp=sf.Reader(path_sh+"/Alineamiento_Galeria.shp" )


limitArray = np.array(ModelLimitShp.shapeRecords()[0].shape.points)
galleryArray = np.array(GalleryShp.shapeRecords()[0].shape.points)
# galleryArray=np.array([point.shape.points[0] for point in GalleryShp.shapeRecords()])

LBB=GloRefBox=ModelLimitShp.bbox
GBB=LocRefBox=GalleryShp.bbox

# To plot Area model
fig = plt.figure()
plt.plot(limitArray[:,0],limitArray[:,1])
plt.arrow(galleryArray[0,0],galleryArray[0,1],-galleryArray[0,0]+galleryArray[1,0],-galleryArray[0,1]+galleryArray[1,1], shape="full")
plt.plot([GBB[0],GBB[0],GBB[2],GBB[2],GBB[0]],[GBB[1],GBB[3],GBB[3],GBB[1],GBB[1]])
plt.plot([LBB[0],LBB[0],LBB[2],LBB[2],LBB[0]],[LBB[1],LBB[3],LBB[3],LBB[1],LBB[1]])
plt.show()

# if you want to see it on a map
# crs={'init' : 'epsg:3116'}#deprecated
mplleaflet.show(fig, epsg=3116)#muy lento por eso lo comento

print(GloRefBox)
print(LocRefBox)

#Calculating Global Model (Glo) and Local Refinement (Loc) dimensions
GloLx = GloRefBox[2] - GloRefBox[0] #x_max - x_min
GloLy = GloRefBox[3] - GloRefBox[1]
print('Global Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (GloLx,GloLy))

LocLx = LocRefBox[2] - LocRefBox[0] #x_max - x_min
LocLy = LocRefBox[3] - LocRefBox[1]
print('Local Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (LocLx,LocLy))


#Defining Global and Local Refinements, for purpose of simplicity cell x and y dimension will be the same
celGlo = 50
celRef = 20

def arrayGeneratorCol(gloRef, locRef, gloSize, locSize):

    cellArray = np.array([])

    while cellArray.sum() + gloRef[0] < locRef[0] - celGlo:
        cellArray = np.append(cellArray,[gloSize])
    while cellArray.sum() + gloRef[0] > locRef[0] - celGlo and cellArray.sum() + gloRef[0] < locRef[2] + celGlo:
        cellArray = np.append(cellArray,[locSize])
    while cellArray.sum() + gloRef[0] > locRef[2] + celGlo and cellArray.sum() + gloRef[0] < gloRef[2]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray
def arrayGeneratorRow(gloRef, locRef, gloSize, locSize):

    cellArray = np.array([])
    accumCoordinate =  gloRef[3] - cellArray.sum()

    while gloRef[3] - cellArray.sum() > locRef[3] + celGlo:
        cellArray = np.append(cellArray,[gloSize])
    while gloRef[3] - cellArray.sum() < locRef[3] + celGlo and gloRef[3] - cellArray.sum() > locRef[1] - celGlo:
        cellArray = np.append(cellArray,[locSize])
    while gloRef[3] - cellArray.sum() < locRef[1] - celGlo and gloRef[3] - cellArray.sum() > gloRef[1]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray


#Remember that DELR is space between rows, so it is the column dimension array
delRArray = arrayGeneratorCol(GloRefBox, LocRefBox, celGlo, celRef)
print("delRows= \n"+str(delRArray))

#And DELC is the space between columns, so its the row dimension array
delCArray = arrayGeneratorRow(GloRefBox, LocRefBox, celGlo, celRef)
print("delColds= \n"+str(delCArray))


#Calculating number or rows and cols since they are dependant from the discretization
nrows = delCArray.shape[0]
ncols = delRArray.shape[0]
print('Number of rows: %d and number of cols: %d' % (nrows,ncols))

#Define some parameters and values for the spatial and temporal discretization (DIS package)

#Number of layers and layer elevations
nlay = 20
mtop = 0
botm = np.array([])
for i in range(nlay):
    botm=np.append(botm, mtop-20)#Corregir/Revisar
print(botm)

# Apply the spatial and temporal discretization parameters to the DIS package
dis = fp.mf6.ModflowGwfdis(gwf, pname= "dis", nlay=nlay,
                           nrow=nrows, ncol=ncols,  delr=delRArray,
                           delc=delCArray, top=mtop, botm=botm,
                           filename= f"{model_name}.dis",xorigin=GloRefBox[0],
                           yorigin=GloRefBox[1])


#  procedure to include Recharge in shapely

recarga=sf.Reader(path_sh+"/Zonas_Rec3")
recarga1=recarga.shapeRecords()[0]
first=recarga1.shape.__geo_interface__
print(first)
shp_geom=shape(first)
# print(shp_geom)
print(type(shp_geom))


# now we use shapely to instersect
# plt.figure()
gwf.modelgrid.plot()
ix = GridIntersect(gwf.modelgrid, method="structured", rtree=True)
# %timeit ix.intersect(shp_geom) #it works!
result=ix.intersect(shp_geom)
print(result)

rch_spd=[]
for i in range(result.shape[0]):
    rch_spd.append([0,*result["cellids"][i],0.11/86400])#hay que revisar si la tupla quedó mal
    
rch=fp.mf6.ModflowGwfrch(gwf, stress_period_data=rch_spd,
                            filename=f"{model_name}.rch",pname="RCH",
                            print_input=True,print_flows=True,save_flows=True)

# we are including the idomain by using shapes
dominio=sf.Reader(path_sh+"/Domin_Mod.shp")
dominio1=dominio.shapeRecords()[0]
firstd=dominio1.shape.__geo_interface__
print(firstd)
shp_geomd=shape(firstd)


resultd= ix.intersect(shp_geomd)
domain_cells=[]
for i in range(result.shape[0]):
    domain_cells.append([0,*resultd["cellids"][i]])#hay que revisar si la tupla quedó mal

# creating the idomain matrix
idom=np.zeros((nlay,nrows,ncols))
print(nlay,nrows,ncols)
for i, value in enumerate(domain_cells):
    # print(i)
    # print(value)
    idom[tuple(value)]=1
    
    
# cloning in all layers
idom[:]=idom[0]

# assigning domain to the model
dis.idomain=idom


# cretaing drains from shapefiles

queb=sf.Reader(os.path.join(path_sh,"Queb_corr.shp"))


queb1=queb.shapeRecords()[0]
firstq=queb1.shape.__geo_interface__
shp_geomq=shape(firstq)

for i in range(1,queb.numRecords):
    queb1=queb.shapeRecords()[i]
    firstq=queb1.shape.__geo_interface__
    shp_geomq=shp_geomq.union(shape(firstq))
    
    
resultq=ix.intersect(shp_geomq)
drn_spd=[]
for i in range(resultq.shape[0]):
    drn_spd.append([0,*resultq["cellids"][i]])#falta agregar valores de quebradas, I need terrain






























