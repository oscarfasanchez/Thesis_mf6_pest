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
from shapely.geometry import Polygon, Point, LineString, MultiLineString, MultiPoint, MultiPolygon, shape#lastone from StackO.F.
from shapely.strtree import STRtree
from flopy.utils.gridintersect import GridIntersect

from flopy.utils import Raster

import copy


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
# 2/21/2021 real rain ideam until there

nsper=int(365*4.5) #number of  stress periods, 365 per year
#there are 86400 seconds
time_disc=[(86400, 1, 1.0) for _ in range(nsper-1)]#[(1,1,1.0)]
time_disc.insert(0,(1,1,1.0))#inserting the steady stress period at the beginning of list
tdis= fp.mf6.ModflowTdis(sim, pname="tdis",
                         time_units="SECONDS", 
                         nper=nsper, perioddata=time_disc,)
# Create gwf model
model_nam_file=f"{model_name}.nam"
gwf = fp.mf6.ModflowGwf(sim, modelname=model_name,
                        model_nam_file=model_nam_file,
                        newtonoptions=None, save_flows=True)


# Setting the solver
ims=  fp.mf6.modflow.mfims.ModflowIms(sim, pname="ims",
                                      complexity= "Complex")


# open shapefiles of limits and refinement
path_sh= "../02_Shp_Vect/Input_Model"
ruta=os.path.join(path_sh,"Domin_Mod.shp")
print(ruta)
ModelLimitShp=sf.Reader(ruta)
GalleryShp=sf.Reader(path_sh+"/Alineamiento_Galeria.shp" )


limitArray = np.array(ModelLimitShp.shapeRecords()[0].shape.points)
galleryArray = np.array(GalleryShp.shapeRecords()[0].shape.points)
# galleryArray=np.array([point.shape.points[0] for point in GalleryShp.shapeRecords()])

LBB=GloRefBox=ModelLimitShp.bbox
GBB=LocRefBox=GalleryShp.bbox

GloRefBox=np.array(GloRefBox)
LocRefBox=np.array(LocRefBox)
# To plot Area model
fig = plt.figure()
plt.plot(limitArray[:,0],limitArray[:,1])
plt.arrow(galleryArray[0,0],galleryArray[0,1],-galleryArray[0,0]+galleryArray[1,0],-galleryArray[0,1]+galleryArray[1,1], shape="full")
plt.plot([GBB[0],GBB[0],GBB[2],GBB[2],GBB[0]],[GBB[1],GBB[3],GBB[3],GBB[1],GBB[1]])
plt.plot([LBB[0],LBB[0],LBB[2],LBB[2],LBB[0]],[LBB[1],LBB[3],LBB[3],LBB[1],LBB[1]])
plt.show()

# if you want to see it on a map
# crs={'init' : 'epsg:3116'}#deprecated
# mplleaflet.show(fig, epsg=3116)#muy lento por eso lo comento

print(GloRefBox)
print(LocRefBox)

#Calculating Global Model (Glo) and Local Refinement (Loc) dimensions
GloLx = GloRefBox[2] - GloRefBox[0] #x_max - x_min
GloLy = GloRefBox[3] - GloRefBox[1] #y_max - y_min
print('Global Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (GloLx,GloLy))

LocLx = LocRefBox[2] - LocRefBox[0] #x_max - x_min
LocLy = LocRefBox[3] - LocRefBox[1] #y_max - y_min
print('Local Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (LocLx,LocLy))

# y si quiero dejar diagonal la malla?t
#Defining Global and Local Refinements, for purpose of simplicity cell x and y dimension will be the same
celGlo = 80
celRef = 3
"""
def arrayGeneratorCol(gloRef, locRef, gloSize, locSize):

    cellArray = np.array([])

    while cellArray.sum() + gloRef[0] < locRef[0] - gloSize:
        cellArray = np.append(cellArray,[gloSize])
    while cellArray.sum() + gloRef[0] > locRef[0] - gloSize and cellArray.sum() + gloRef[0] < locRef[2] + gloSize:
        cellArray = np.append(cellArray,[locSize])
    while cellArray.sum() + gloRef[0] > locRef[2] + gloSize and cellArray.sum() + gloRef[0] < gloRef[2]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray

def arrayGeneratorRow(gloRef, locRef, gloSize, locSize):

    cellArray = np.array([])

    while gloRef[3] - cellArray.sum() > locRef[3] + gloSize:
        cellArray = np.append(cellArray,[gloSize])
    while gloRef[3] - cellArray.sum() < locRef[3] + gloSize and gloRef[3] - cellArray.sum() > locRef[1] - gloSize:
        cellArray = np.append(cellArray,[locSize])
    while gloRef[3] - cellArray.sum() < locRef[1] - gloSize and gloRef[3] - cellArray.sum() > gloRef[1]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray
"""
import math
def arrayGeneratorRowSmooth(gloRef, locRef, gloSize, locSize):
    print("coordGlob= "+str( gloRef[3]))
    smooth=1.5#maximun cell size incremente between cells
    cellratio=gloSize/locSize
    # calculate how many cells(and lenght) do we need to reach big cell size
    cellneeds=np.log(cellratio)/np.log(smooth)#how many smooth cells are needed to change the cell size correctly
    print("smooth_cellneeds= "+str(cellneeds))
    cellcomplete=int(cellneeds)#steps needed to change cell size smoothly
    print("cellcomplete= "+str(cellcomplete))
    lenght_afected=0
    smooth_cell= np.array([])
    for i in range(cellcomplete):
        lenght_afected += locSize*smooth**(i+1)#lenght needed for transition to small size
        print("lenght_afected="+str(lenght_afected))
        smooth_cell=np.append(smooth_cell,locSize*smooth**(i+1))#smooth cells size
    # lenght_afected=(cellneeds-int(cellneeds))*gloSize
    print("lenght_afected="+str(lenght_afected))
    print("smoothcell"+str(smooth_cell))    
    cells_afected=(lenght_afected/gloSize)-lenght_afected//gloSize#partial big cell affected
    print("cells_afected="+str(cells_afected))
    # lenght_afected =+ locSize
    # correct grid position due to grid change
    delta_coords=gloSize*cells_afected#coordinate shift to adapt mesh
    print("deltaCoords:"+str(delta_coords))
    print("coordGlob= "+str( gloRef[3]))
    # cloning the list to avoid side effects
    gloRef2=list(gloRef)
    gloRef2[3] =gloRef[3] + delta_coords#coordianate corrected
    print("coordGlobCorr=",gloRef2[3])
    cellArray = np.array([])
    
    print("cellarray:", cellArray)
    
    while gloRef2[3] - cellArray.sum() > locRef[3] + gloSize + lenght_afected:
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + gloSize + lenght_afected, gloRef2[3] - cellArray.sum() > locRef[3] + gloSize + lenght_afected)
        cellArray = np.append(cellArray,[gloSize])
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + gloSize + lenght_afected, gloRef2[3] - cellArray.sum() > locRef[3] + gloSize + lenght_afected)
        
    
    i=1    
    while gloRef2[3] - cellArray.sum() > locRef[3] + smooth_cell[cellcomplete-i] and cellcomplete-i >= 0:
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + gloSize, gloRef2[3] - cellArray.sum() > locRef[3] + gloSize)
        # print(i," i/cellcom" ,cellcomplete)
        cellArray = np.append(cellArray,smooth_cell[cellcomplete-i])
        # print("smoothCell: "+str(smooth_cell[cellcomplete-i]))
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + gloSize, gloRef2[3] - cellArray.sum() > locRef[3] + gloSize)
        i += 1
        
        
    while gloRef2[3] - cellArray.sum() < locRef[3] + gloSize and gloRef2[3] - cellArray.sum() > locRef[1] - gloSize:
        cellArray = np.append(cellArray,[locSize])
        
    i=0      
    while gloRef2[3] - cellArray.sum() < locRef[1] + gloSize and gloRef2[3] - cellArray.sum() > locRef[1] - lenght_afected - smooth_cell[cellcomplete - 1]:
        cellArray = np.append(cellArray,smooth_cell[i])
        i += 1
        
    while gloRef2[3] - cellArray.sum() < locRef[1] - gloSize and gloRef2[3] - cellArray.sum() > gloRef2[1]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray


def arrayGeneratorColSmooth(gloRef, locRef, gloSize, locSize):
    print("coordGlob= "+str( gloRef[0]))
    smooth=1.5
    cellratio=gloSize/locSize
    # calculate how many cells(and lenght) do we need to reach big cell size
    cellneeds=np.log(cellratio)/np.log(smooth)
    print("smooth_cellneeds= "+str(cellneeds))
    cellcomplete=int(cellneeds)
    print("cellcomplete= "+str(cellcomplete))
    lenght_afected=0
    smooth_cell= np.array([])
    for i in range(cellcomplete):
        lenght_afected += locSize*smooth**(i+1)
        print("lenght_afected="+str(lenght_afected))
        smooth_cell=np.append(smooth_cell,locSize*smooth**(i+1))
    # lenght_afected=(cellneeds-int(cellneeds))*gloSize
    print("lenght_afected="+str(lenght_afected))
    print("smoothcell"+str(smooth_cell))    
    cells_afected=(lenght_afected/gloSize)-lenght_afected//gloSize
    print("cells_afected="+str(cells_afected))
    # lenght_afected =+ locSize
    # correct grid position due to grid change
    delta_coords=gloSize*cells_afected
    print("deltaCoords:"+str(delta_coords))
    print("coordGlob= "+str( gloRef[0]))
    # cloning the list to avoid side effects
    gloRef2=list(gloRef)
    gloRef2[0] =gloRef[0] + delta_coords
    print("coordGlobCorr=",gloRef2[0])
    cellArray = np.array([])
    
    print("cellarray:", cellArray)
    
    while gloRef2[0] + cellArray.sum() < locRef[0] - gloSize - lenght_afected:
        # print(gloRef2[0] - cellArray.sum()," comparar" ,locRef[0] + gloSize + lenght_afected, gloRef2[0] - cellArray.sum() > locRef[0] + gloSize + lenght_afected)
        cellArray = np.append(cellArray,[gloSize])
        # print(gloRef2[0] - cellArray.sum()," comparar" ,locRef[0] + gloSize + lenght_afected, gloRef2[0] - cellArray.sum() > locRef[0] + gloSize + lenght_afected)
        
    
    i=1    
    while gloRef2[0] + cellArray.sum() < locRef[0] - smooth_cell[cellcomplete-i] and cellcomplete-i >= 0:
        # print(gloRef2[0] - cellArray.sum()," comparar" ,locRef[0] + gloSize, gloRef2[0] - cellArray.sum() > locRef[0] + gloSize)
        # print(i," i/cellcom" ,cellcomplete)
        cellArray = np.append(cellArray,smooth_cell[cellcomplete-i])
        # print("smoothCell: "+str(smooth_cell[cellcomplete-i]))
        # print(gloRef2[0] - cellArray.sum()," comparar" ,locRef[0] + gloSize, gloRef2[0] - cellArray.sum() > locRef[0] + gloSize)
        i += 1
        
        
    while gloRef2[0] + cellArray.sum() > locRef[0] - gloSize and gloRef2[0] + cellArray.sum() < locRef[2] + gloSize:
        cellArray = np.append(cellArray,[locSize])
        
    i=0      
    while gloRef2[0] + cellArray.sum() > locRef[2] - gloSize and gloRef2[0] + cellArray.sum() < locRef[2] + lenght_afected + smooth_cell[cellcomplete - 1]:
        cellArray = np.append(cellArray,smooth_cell[i])
        i += 1
        
    while gloRef2[0] + cellArray.sum() > locRef[2] + gloSize and gloRef2[0] + cellArray.sum() < gloRef2[2]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray

#And DELC is the space between columns, so its the row dimension array
delCArray = arrayGeneratorRowSmooth(GloRefBox, LocRefBox, celGlo, celRef)
print("delCols= \n"+str(delCArray))


#Remember that DELR is space between rows, so it is the column dimension array
delRArray = arrayGeneratorColSmooth(GloRefBox, LocRefBox, celGlo, celRef)
print("delRows= \n"+str(delRArray))



# delCArray = arrayGeneratorRow(GloRefBox, LocRefBox, celGlo, celRef)
# print("delCols= \n"+str(delCArray))

#Calculating number or rows and cols since they are dependant from the discretization
nrows = delCArray.shape[0]
ncols = delRArray.shape[0]
print('Number of rows: %d and number of cols: %d' % (nrows,ncols))

#Define some parameters and values for the spatial and temporal discretization (DIS package)

#Number of layers and layer elevations
nlay = 7
mtop = 1000
H=500 #model bottom
botm= np.linspace(mtop-H/nlay, H, nlay)



# Apply the spatial and temporal discretization parameters to the DIS package
dis = fp.mf6.ModflowGwfdis(gwf, pname= "dis", nlay=nlay,
                           nrow=nrows, ncol=ncols,  delr=delRArray,
                           delc=delCArray, top=mtop, botm=botm,
                           filename= f"{model_name}.dis",xorigin=GloRefBox[0],
                           yorigin=GloRefBox[1])



# assigning surface raster
path_raster="../03_Raster/Input_ModelR/Superficies_R_Tiff"

surface=Raster.load(os.path.join(path_raster,"R_Topo_Union_Clip.tif"))
surface.bands

# fig2=plt.figure(figsize=(12,12))
# ax=fig2.add_subplot(1,1,1, aspect="equal")
# ax=surface.plot(ax=ax)
# plt.colorbar(ax.images[0], shrink=0.7) 
# gwf.modelgrid.plot()
# intersecting and resampling raster
dem_Matrix=surface.resample_to_grid(gwf.modelgrid, surface.bands[0], method="nearest")
dis.top=dem_Matrix
"""
botm=np.empty((nlay,nrows,ncols))
print(botm)

botm[0,:,:]=dem_Matrix-100
for i in range(1,nlay):
    botm[i,:,:]=botm[i-1,:,:]-(H/nlay)#Corregir/Revisar/provisional
    
dis.botm=botm
"""
#hydraulic conductivity definition in arrays

k_qd_qbg=1e-7*np.ones([nrows,ncols])
k_1="k_qd_qbg.txt"
k_qbo2=1e-7*np.ones([nrows,ncols])
k_2="k_qbo2.txt"
k_qbo1=1e-8*np.ones([nrows,ncols])
k_3="k_qbo1.txt"
k_roc=1e-7*np.ones([nrows,ncols])
k_4="k_roc.txt"

kv_qd_qbg=1e-1*np.ones([nrows,ncols])
kv_1="kv_qd_qbg.txt"
kv_qbo2=1e-1*np.ones([nrows,ncols])
kv_2="kv_qbo2.txt"
kv_qbo1=1e-1*np.ones([nrows,ncols])
kv_3="kv_qbo1.txt"
kv_roc=1e-1*np.ones([nrows,ncols])
kv_4="kv_roc.txt"

#Storage definition

sy_qd_qbg=1e-1*np.ones([nrows,ncols])
sy_1="sy_qd_qbg.txt"
sy_qbo2=1e-1*np.ones([nrows,ncols])
sy_2="sy_qbo2.txt"
sy_qbo1=1e-1*np.ones([nrows,ncols])
sy_3="sy_qbo1.txt"
sy_roc=1e-1*np.ones([nrows,ncols])
sy_4="sy_roc.txt"

ss_qd_qbg=1e-4*np.ones([nrows,ncols])
ss_1="ss_qd_qbg.txt"
ss_qbo2=1e-4*np.ones([nrows,ncols])
ss_2="ss_qbo2.txt"
ss_qbo1=1e-4*np.ones([nrows,ncols])
ss_3="ss_qbo1.txt"
ss_roc=1e-5*np.ones([nrows,ncols])
ss_4="ss_roc.txt"



# we are including the qbg hydraulic conductivity  by using shapes
gravoso=sf.Reader(path_sh+"/Superficies/Cont_Qbg2.shp")
gravoso1=gravoso.shapeRecords()[0]
firstkg=gravoso1.shape.__geo_interface__
# print(firstd)
shp_geomkg=shape(firstkg)


# now we use shapely to instersect
# plt.figure()
gwf.modelgrid.plot()
ix = GridIntersect(gwf.modelgrid, method="structured", rtree=False)

resultkg= ix.intersect(shp_geomkg)
gravoso_cells=[]
for i in range(resultkg.shape[0]):
    gravoso_cells.append([*resultkg["cellids"][i]])#hay que revisar si la tupla quedó mal

# creating the idomain matrix
domin_grav=np.zeros((nrows,ncols))
# print(nlay,nrows,ncols)
for i, value in enumerate(gravoso_cells):
    # print(i)
    # print(value)
    k_qd_qbg[tuple(value)]=1e-5
    ss_qd_qbg[tuple(value)]=1e-4
    sy_qd_qbg[tuple(value)]=1e-1
    domin_grav[tuple(value)]=1 #this will be necessary for pyemu 
#this should be more automatic i guess
np.savetxt(os.path.join(workspace,k_1), k_qd_qbg)
np.savetxt(os.path.join(workspace,k_2), k_qbo2)
np.savetxt(os.path.join(workspace,k_3), k_qbo1)
np.savetxt(os.path.join(workspace,k_4), k_roc)
kgeol=[k_1, k_2, k_3, k_4]

np.savetxt(os.path.join(workspace,kv_1), kv_qd_qbg)
np.savetxt(os.path.join(workspace,kv_2), kv_qbo2)
np.savetxt(os.path.join(workspace,kv_3), kv_qbo1)
np.savetxt(os.path.join(workspace,kv_4), kv_roc)
kvgeol=[kv_1, kv_2, kv_3, kv_4]

np.savetxt(os.path.join(workspace,ss_1), ss_qd_qbg)
np.savetxt(os.path.join(workspace,ss_2), ss_qbo2)
np.savetxt(os.path.join(workspace,ss_3), ss_qbo1)
np.savetxt(os.path.join(workspace,ss_4), ss_roc)
ssgeol=[ss_1, ss_2, ss_3, ss_4]

np.savetxt(os.path.join(workspace,sy_1), sy_qd_qbg)
np.savetxt(os.path.join(workspace,sy_2), sy_qbo2)
np.savetxt(os.path.join(workspace,sy_3), sy_qbo1)
np.savetxt(os.path.join(workspace,sy_4), sy_roc)
sygeol=[sy_1, sy_2, sy_3, sy_4]


def dis_layers(path_folder, name_raster, div_layers, kgeo, kvgeo, ssgeo, sygeo, bottom_model = -1000, min_thick=0):
    """
    dis_layers is made for layer discretization
    this function also helps to define idomain -1 when layer thickness is 0
    also helps to define hydraulic conductivity
    raster layers should have a possitive no data value to be corrected by minimun thickness
    
    path_folder is a string with the relative or absolute path of the folder containing raster files
    name_raster is a list of  variables containing the names of raster files, the firxt one has to be model Top
    div_layers is a list of int variables containing the number of division of each layer between raster layers
    it should have the same size as name_raster minus 1(or the same if bottom constant layer is defined[bottom <= -1000])
    bottom is the bottom height  of the model, if negative, the last raster will be used as bottom
    min_thick is a list of int Minimun raster thickness variable allowed for raster layers.
    
    Return
    botm is a matrix of each layer bottom
    demMatrix is a matrix of geological layers bottoms
    thickcells is a matrix used for identify no thick cells, later can be used to build an idomain matrix
    k, kv hydraulic condivity matrix for every layer
    
    """
    
    # deciding if there is an extra layer for model bottom
    if bottom_model == -1000: 
        botm=np.empty((div_layers.sum(),nrows, ncols))
        demMatrix=np.empty((len(name_raster),nrows, ncols))
        
    else:
        botm=np.empty((div_layers.sum(),nrows, ncols))
        demMatrix=np.empty((len(name_raster) + 1 ,nrows, ncols))
    
    # creating an empty list of k
    k=div_layers.sum()*[None]
    kv=div_layers.sum()*[None]
    # creating an empty list of Storage
    ss=div_layers.sum()*[None]
    sy=div_layers.sum()*[None]
    # k=list()
    
    print(demMatrix)
    print(demMatrix.ndim)
    print(demMatrix.shape)
    # first define model top
    Topo = Raster.load(os.path.join(path_folder, name_raster[0]))
    fig = plt.figure(figsize=(12,12))
    ax=fig.add_subplot(1, 1, 1, aspect = "equal")
    ax = Topo.plot(ax = ax)
    plt.colorbar(ax.images[0], shrink =0.7)
    gwf.modelgrid.plot()
    # intersecting and resampling raster
    demMatrix[0]=Topo.resample_to_grid(gwf.modelgrid,
                                    Topo.bands[0], method="nearest")
    print(demMatrix)
    print(type(demMatrix))
    print(demMatrix.ndim)
    print(demMatrix.shape)
    # dis.top=demMatrix[0]
    count = 0
    
    print("bottom_Shape=",botm.shape)
    thickcells=np.ones(botm.shape)
    for i in range(1, len(capas)+1):
        print("i=",i)
        # assign geologicalbottom layers to matrix
        if i  == len(capas) and bottom_model != -1000:
            demMatrix[i] = bottom_model
            print("opt1")
        elif i  == len(capas) and bottom_model == -1000:
            print("you probably didn't define a suitable bottom model!!")
            demMatrix[i] = bottom_model
            print("opt2")
        else:
            print("opt3")
            bottom = Raster.load(os.path.join(path_folder, name_raster[i]))
            fig = plt.figure(figsize=(12,12))
            ax=fig.add_subplot(1, 1, 1, aspect = "equal")
            ax = bottom.plot(ax = ax)
            plt.colorbar(ax.images[0], shrink =0.7)
            gwf.modelgrid.plot()
            # intersecting and resampling raster, the nodata values, are being used, it has to be fixed
            demMatrix[i] = bottom.resample_to_grid(gwf.modelgrid,
                                                  Topo.bands[0],
                                                  method="linear",multithread=True,
                                                  thread_pool=6, extrapolate_edges=True)
        
        
        
        demMatrix[i][demMatrix[i] + min_thick[i-1] > demMatrix[i-1]] = demMatrix[i-1][demMatrix[i] + min_thick[i-1] > demMatrix[i-1]]-min_thick[i-1]
        # demMatrix[i][demMatrix[i]  > demMatrix[i-1]] = demMatrix[i-1][demMatrix[i] > demMatrix[i-1]]
        
        print("count=",count)
        for j in range(count ,div_layers[i-1] + count):
            # assign model layer bottoms, K and thickcell to a matrix
            print("j=",j)
            if div_layers[i-1] == 1:
                botm[j]= demMatrix[i]
                thickcells[j][demMatrix[i] + 0 >= demMatrix[i-1]]=-1
                k[j]=kgeo[i-1]
                kv[j]=kvgeo[i-1]
                ss[j]=ssgeo[i-1]
                sy[j]=sygeo[i-1]
                break
            
            botm[j] = demMatrix[i-1]+(demMatrix[i]-demMatrix[i-1])*(j-count + 1)/div_layers[i-1]
            thickcells[j][demMatrix[i] + 0>= demMatrix[i-1]]=-1
            k[j]=kgeo[i-1]
            kv[j]=kvgeo[i-1]
            ss[j]=ssgeo[i-1]
            sy[j]=sygeo[i-1]
        count += div_layers[i-1]
        print("count=",count)
    # if bottom_model != -1000:
    #     botm[-1]= bottom_model
    
    
    
    return botm, demMatrix, thickcells, k, kv, ss, sy

# run dis_layer function
raster_names=["R_Topo_Union_Clip.tif", "R_Qbg_Qd.tif", "R_Qbo2.tif", "R_Qbo1.tif"]
capas=np.array([1,2,2,1])
min_thick=([1,0,0,0])
fondos, geol, thickcells, k, kv, ss, sy= dis_layers(path_raster,raster_names, capas, kgeol, kvgeol, ssgeol, sygeol, bottom_model=500, min_thick=min_thick)
nlay = fondos.shape[0]
dis.nlay = fondos.shape[0]
dis.botm=fondos
celltype=([0,0,0,1])#not used, ill defined
# define node property flow package
npf = fp.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, k33overk=True, k33=kv,
                           save_flows=True, save_specific_discharge = True)
# specifies storage
sto = fp.mf6.ModflowGwfsto(gwf, pname="sto", save_flows=True, iconvert=1, ss=ss, sy=sy, steady_state= {0:True}, transient={1:True})
#later i will need it to pyemu
np.savetxt(os.path.join(workspace,"layers"), capas)

"""       

# assigning bottom height using geology

bot_Qd = Raster.load(os.path.join(path_raster, "R_Qd.tif"))
bot_Qd.bands

figQd = plt.figure(figsize=(12,12))
ax = figQd.add_subplot(1,1,1, aspect = "equal")
ax= bot_Qd.plot(ax=ax)
plt.colorbar(ax.images[0], shrink =0.7)
gwf.modelgrid.plot()
# intersecting and resampling raster
bot_qd_matrix = bot_Qd.resample_to_grid(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters, bot_Qd.bands[0], method= "nearest")

botm[0][botm[0] < dem_Matrix] = bot_qd_matrix[botm[0] < dem_Matrix]
botm[0][botm[0] > dem_Matrix] = dem_Matrix[botm[0] > dem_Matrix]
dis.botm=botm


bot_Qbo2 = Raster.load(os.path.join(path_raster, "R_Qbo2.tif"))
bot_Qbo2.bands

figQbo2 = plt.figure(figsize=(12,12))
ax = figQbo2.add_subplot(1,1,1, aspect = "equal")
ax= bot_Qbo2.plot(ax=ax)
plt.colorbar(ax.images[0], shrink =0.7)
gwf.modelgrid.plot()
# intersecting and resampling raster
bot_qbo2_matrix = bot_Qbo2.resample_to_grid(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters, bot_Qbo2.bands[0], method= "nearest")

botm[1][botm[1] < dem_Matrix] = bot_qbo2_matrix[botm[1] < dem_Matrix]
botm[1][botm[1] > dem_Matrix] = dem_Matrix[botm[1] > dem_Matrix]
dis.botm=botm



bot_Qbo1 = Raster.load(os.path.join(path_raster, "R_Qbo1.tif"))
bot_Qbo1.bands

figQbo1 = plt.figure(figsize=(12,12))
ax = figQbo1.add_subplot(1,1,1, aspect = "equal")
ax= bot_Qbo1.plot(ax=ax)
plt.colorbar(ax.images[0], shrink =0.7)
gwf.modelgrid.plot()
# intersecting and resampling raster
bot_qbo1_matrix = bot_Qbo1.resample_to_grid(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters, bot_Qbo1.bands[0], method= "nearest")

botm[3][botm[3] < dem_Matrix] = bot_qbo1_matrix[botm[3] < dem_Matrix]
botm[3][botm[3] > dem_Matrix] = dem_Matrix[botm[3] > dem_Matrix]
dis.botm=botm


"""


#  procedure to include Recharge 
df_rain=pd.read_csv("../04_Xls/Lluvia_Ideam.csv", sep=",")#load file
df_rain["Fecha"]=pd.to_datetime(df_rain.Fecha, dayfirst=False)#format of of dates setted
df_rain["time"]=df_rain["Fecha"]-df_rain["Fecha"][0]#days since the first day
df_rain["time_s"]=df_rain["time"].astype("timedelta64[s]")#count in seconds for modflow

#shapely process 
recarga=sf.Reader(path_sh+"/Zonas_Rec3")
recarga1=recarga.shapeRecords()[0]
first=recarga1.shape.__geo_interface__
# print(first)
shp_geom=shape(first)
# print(shp_geom)
print(type(shp_geom))



# %timeit ix.intersect(shp_geom) #it works!
result=ix.intersect(shp_geom)
print(result)
# Recharge is weighted by intersected Area and defined as Recharge boundary (Neuman) condition (BC)
rch_spd=[]
#steady list for first stress period
rch_spd_st=[]
for i in range(result.shape[0]):#0.11 from water budget?REVIEW VALUE IN TRANSIENT
    rch_spd.append([0,
                    *result["cellids"][i],
                    0.1*(0.001/86400)*(result['areas'][i]/#convert mm/day to m/s , 10% rain percolated*1/1000(due to mm units in files))
                                  delCArray[result["cellids"][i][0]]/
                                  delRArray[result["cellids"][i][1]]),"rain_mult"])#it is weighted by intersected area, rain_mult timeseries multplier
    rch_spd_st.append([0,
                    *result["cellids"][i],
                    (110/1000/86400/365)*(result['areas'][i]/#convert mm/year to m/s ,10%*11000 mm --> 110 mm, steady
                                  delCArray[result["cellids"][i][0]]/
                                  delRArray[result["cellids"][i][1]]), 1])#it is multiplied by %cell area afected
for i in range(len(rch_spd)):#to correct 0 based index in python, because flopy doesn't correct in external files
    rch_spd[i][0]+=1
    rch_spd[i][1]+=1
    rch_spd[i][2]+=1
    rch_spd_st[i][0]+=1
    rch_spd_st[i][1]+=1
    rch_spd_st[i][2]+=1

df_rch_spd_st=pd.DataFrame(rch_spd_st)
df_rch_spd=pd.DataFrame(rch_spd)
df_rch_spd_st.to_csv(workspace+"/rch_0.txt", index=False, header= False, sep= " ")
df_rch_spd.to_csv(workspace+"/rch_1.txt", index=False, header= False, sep= " ")
   #we start setting "rain_mult" 
ts_data=[]
ts_data=[(df_rain["time_s"].tolist()[i], df_rain["Valor"].tolist()[i]) for i in range(df_rain.shape[0])]

rch_spd_txt={0:{"filename":"rch_0.txt"},1:{"filename":"rch_1.txt"}}

ts_dict={
    "timeseries": ts_data,
    "time_series_namerecord":"rain_mult",
    "interpolation_methodrecord":"LINEAREND"}


# debo revisar la activación de recarga en celdas secas.    
rch=fp.mf6.ModflowGwfrch(gwf, stress_period_data=rch_spd_txt,
                            filename=f"{model_name}.rch",pname="RCH",
                            auxiliary="rain_mult",
                            auxmultname="rain_mult",
                            timeseries= ts_dict,
                            print_input=True,print_flows=True,save_flows=True)

# we are including the idomain by using shapes
dominio=sf.Reader(path_sh+"/Domin_Mod.shp")
dominio1=dominio.shapeRecords()[0]
firstd=dominio1.shape.__geo_interface__
# print(firstd)
shp_geomd=shape(firstd)


resultd= ix.intersect(shp_geomd)
domain_cells=[]
for i in range(resultd.shape[0]):
    domain_cells.append([0,*resultd["cellids"][i]])#hay que revisar si la tupla quedó mal

# creating the idomain matrix
idom=np.zeros((nlay,nrows,ncols))
# print(nlay,nrows,ncols)
# Assigning the list of active cells to a matrix
for i, value in enumerate(domain_cells):
    # print(i)
    # print(value)
    idom[tuple(value)]=1
    
    
# cloning in all layers
idom[:]=idom[0]
idom=idom*thickcells #using this variable to use vertical flow cells
# assigning domain to the model
dis.idomain=idom



# creating drains from shapefiles

queb=sf.Reader(os.path.join(path_sh,"Queb_corr.shp"))


queb1=queb.shapeRecords()[0]
firstq=queb1.shape.__geo_interface__
shp_geomq=shape(firstq)

for i in range(1,queb.numRecords):
    queb1=queb.shapeRecords()[i]
    firstq=queb1.shape.__geo_interface__
    shp_geomq=shp_geomq.union(shape(firstq))
    
    #change conductance accordingly to area/lenght?
resultq=ix.intersect(shp_geomq)
drn_spd=[]
for i in range(resultq.shape[0]):
    if idom[tuple((0,*resultq["cellids"][i]))]==1:
        drn_spd.append([0,*resultq["cellids"][i], dem_Matrix[resultq["cellids"][i]]+1, 1e-5*delCArray[resultq["cellids"][i][0]]*delRArray[resultq["cellids"][i][1]]])#falta agregar valores de quebradas, I need terrain

drn=fp.mf6.ModflowGwfdrn(gwf,stress_period_data=drn_spd, filename=f"{model_name}.drn", pname="drn", print_input=True,print_flows=True,save_flows=True)


# cretaing Rivers from shapefiles as Constant heads

rio=sf.Reader(os.path.join(path_sh,"Rios.shp"))


rio1=rio.shapeRecords()[0]
firstr=rio1.shape.__geo_interface__
shp_geomr=shape(firstr)

for i in range(1,rio.numRecords):
    rio1=rio.shapeRecords()[i]
    firstr=rio1.shape.__geo_interface__
    shp_geomr=shp_geomr.union(shape(firstr))
    
    
resultr=ix.intersect(shp_geomr)
chd_spd=[]
for i in range(resultr.shape[0]):
    if idom[tuple((0,*resultr["cellids"][i]))]==1:
        chd_spd.append([0,*resultr["cellids"][i], dem_Matrix[resultr["cellids"][i]]+1 ])#falta agregar valores de quebradas, I need terrain

chd=fp.mf6.ModflowGwfchd(gwf,stress_period_data=chd_spd, filename=f"{model_name}.chd", pname="chd", print_input=True,print_flows=True,save_flows=True)
# creating the main ghc inflow boundary condition

ghb_spd=[]
ghb_spd_tr=[]
# defining conductances per thickness unit in general head boundary condition

c_1=k_qd_qbg[-2,-1]*(celGlo/2000)#*thickness_layer to get conductance

c_2=k_qbo2[0,0]*(celGlo/2000)

c_3=k_qbo1[0,0]*(celGlo/2000)

c_4=k_roc[0,0]*(celGlo/2000)

# np.savetxt(os.path.join(workspace,c_1), c_qd_qbg)

c_geol=[c_1, c_2, c_3, c_4]

# thickness for conductance
lay_thick  = np.ones(fondos.shape)
lay_thick[0]=dem_Matrix -fondos[0]
for i in range(1, fondos.shape[0]):
    lay_thick[i]=fondos[i-1]-fondos[i]
    
    

j=0
# c=[None]*capas.sum()
# c=[c_geol[0]]*capas[0]
c=[]
layer_count=[]
# create each layer with suitable conductance
for i in range(capas.size):
    for j in range(capas[i]):
        c.append(c_geol[i])#vector of suitable conductance
        layer_count.append(i) #list to remeber discretization of geological layer
        print("i= ",i ," j= ",j )
    
"""
for i in range(capas.sum()):
    if capas[:j].sum()!=i:
        # ghb_spd[1]=capas[j]
        c[i]=c_geol[j]
        fondos[capas[:j].sum()]
        capas[i:].sum()
        print("i= ",i ," j= ",j," capas= ", capas[:j].sum() )
        
    else:
        print("i= ",i ," j= ",j," capas= ", capas[:j].sum() )
        j+=1
 """       

# Processing monthly rain to GHB condition
df_monthly_rain=df_rain.resample("M", on="Fecha").sum()#resampling rain by month
df_monthly_rain["time"]=df_monthly_rain.index-df_rain["Fecha"][0]#dates since first day
df_monthly_rain["time_s"]=df_monthly_rain["time"].astype("timedelta64[s]")#convert to time in seconds
df_monthly_rain=df_monthly_rain.loc[(df_monthly_rain.index<"2021-12-31")]#filter by date
df_monthly_rain["weight"]=df_monthly_rain["Valor"]/df_monthly_rain["Valor"].mean()

inflow=sf.Reader(os.path.join(path_sh,"Chd_In.shp"))

inflow1=inflow.shapeRecords()[0]
firsti=inflow1.shape.__geo_interface__
shp_geomi=shape(firsti)

for i in range(1,inflow.numRecords):#by default
    inflow1=inflow.shapeRecords()[i]
    firsti=inflow1.shape.__geo_interface__
    shp_geomi=shp_geomi.union(shape(firsti))
    
# assign ghb 
# Antolinez assumes a head variation from 1005 to 1038 in Gravoso constant head, terrain height in ghb varies between 990 and 940, 
# steady maximun head in antolinez was 1032-1038
Lay_hei_off=20#layer height offset between virtual head and ghb terrain height (20+990= 1010? as minimun)

ghb_tr_delta=20#value to add to the average head variation

df_monthly_rain["delta_head"]=df_monthly_rain["weight"]*ghb_tr_delta

resulti=ix.intersect(shp_geomi)
for i in range(resulti.shape[0]):
    for j in range(capas.sum()):#condition for dry cells
        if idom[tuple((j,*resulti["cellids"][i]))]==1 and fondos[tuple((j,*resulti["cellids"][i]))]<dem_Matrix[resulti["cellids"][i]]+Lay_hei_off:#set a variable
            ghb_spd.append([j+1,*resulti["cellids"][i], dem_Matrix[resulti["cellids"][i]]+Lay_hei_off,c[j]/lay_thick[tuple((j,*resulti["cellids"][i]))], 1])#problema de la altura -60 que queda debajo de las celdas, hay que calcular conductancia.
            # I use j+1 because i will write modflow file directly, so flopy doesn't make the transition
            # ghb_spd_tr.append([j+1,*resulti["cellids"][i], dem_Matrix[resulti["cellids"][i]]+Lay_hei_off,c[j]/lay_thick[tuple((j,*resulti["cellids"][i]))],"rain_mult"])
for i in range(len(ghb_spd)):
    # to fix index because i will write modflow file directly, so flopy doesn't make the transition
    ghb_spd[i][1]+=1
    ghb_spd[i][2]+=1

df_ghb = pd.DataFrame(ghb_spd)
df_ghb_tr = pd.DataFrame(ghb_spd)
# df_ghb_tr = pd.DataFrame(ghb_spd_tr)
df_monthly_rain["time"]
ghb_spd_txt={}

# ts_data=[]
# ts_data=[(df_monthly_rain["time_s"].tolist()[i],df_monthly_rain["Valor"].tolist()[i]) for i in range(df_monthly_rain.shape[0])]

for i in range(0,capas.sum()):
    df_ghb_tr[df_ghb[0]==i+1].to_csv(workspace + f"/ghb_{i}_{0}.txt", index=False, header=False, sep=' ')
    for j in range(df_monthly_rain.shape[0]):
        df_ghb.loc[df_ghb[0]==i+1,3]=df_ghb_tr.loc[df_ghb[0]==i+1,3]+df_monthly_rain.reset_index()["delta_head"][j]
        df_ghb[df_ghb[0]==i+1].to_csv(workspace + f"/ghb_{i}_{j+1}.txt", index=False, header=False, sep=' ')
    
    # list of sp and list of names
    ghb_spd_txt[0]={"filename":f"ghb_{i}_{0}.txt"}
    ghb_spd_txt.update( dict((df_monthly_rain.reset_index()["time"][j].days,
                       {"filename":f"ghb_{i}_{j+1}.txt"}) for j in range(df_monthly_rain.shape[0])))#
    
    # ghb_spd_txt={0:{"filename":f"ghb_{i}.txt"},1:{"filename":f"ghb_tr_{i}.txt"}}
    print(i)
    
    # ts_dict = {
    #         "timeseries": ts_data,
    #         "time_series_namerecord": "rain_mult",
    #         "interpolation_methodrecord": "LINEAREND",
    #     }
    
    fp.mf6.ModflowGwfghb(gwf,stress_period_data=ghb_spd_txt,
                         filename=f"{model_name}_{i}.ghb", pname=f"ghb_{i}",
                         # auxiliary="rain_mult",
                         # auxmultname="rain_mult",
                         # timeseries =ts_dict,
                         print_input=True,print_flows=True,save_flows=True)
    print(i)

# gallery construction
gal=sf.Reader(os.path.join(path_sh,"Alineamiento_Galeria.shp"))
gal1=gal.shapeRecords()[0]
firstg=gal1.shape.__geo_interface__
shp_geomg=shape(firstg)

for i in range(1,gal.numRecords):
    gal1=gal.shapeRecords()[i]
    firstg=gal1.shape.__geo_interface__
    shp_geomg=shp_geomg.union(shape(firstg))
    
    
resultg=ix.intersect(shp_geomg)
gal_spd=[]

#conductance per meter of gallery
gal_perim=9.2106+1.3906*2+0.46#ignoring depth channel
gal_thick=0.3
gal_k_conc=1e-7
gal_k_grav=1e-7
gal_drn_space=10#space between sections with drains
len_drn=1.5*7#outside gallery in 2 inches, there are seven drains
gal_drn_cnd=gal_k_grav*len_drn*np.pi*2*(2*0.0254)/0.01#K*L*With/thickness
gal_cond=(gal_perim*gal_k_conc/gal_thick)+gal_drn_cnd/gal_drn_space#
height_0=748
height_f=753.25
gal_len=resultg["lengths"].sum()#525
gal_slp=0.01# gallery slope
gal_speed_exc=2.2 #gallery excavation speed m/d



time_gal_0=365*3+30*3#time when gallery construction begins
gal_spd_tr={}
# write comments
for i in range(resultg.shape[0]):#iterate gallery cells
    for j in range(fondos.shape[0]):#iterate over layer numbers
        if height_0+resultg["lengths"][0:i].sum()*gal_slp > fondos[tuple((j,*resultg["cellids"][i]))]:
            if idom[tuple((0,*resultg["cellids"][i]))]==1:#Cambiar J?
                gal_spd.append([j,*resultg["cellids"][i],
                                height_0+resultg["lengths"][0:i].sum()*gal_slp, gal_cond*resultg["lengths"][i],"gal_flow"])
            print(int(resultg["lengths"][0:i+1].sum()//gal_speed_exc)+1)
            gal_spd_tr[int(resultg["lengths"][0:i+1].sum()//gal_speed_exc)+1+time_gal_0]=list(gal_spd)#review indent
            break
        
gal_obs={"mod_drn_gal_obs.csv":[("gal-flow", "drn", "gal_flow")]}        

drn_gal=fp.mf6.ModflowGwfdrn(gwf,stress_period_data=gal_spd_tr,
                             filename=f"{model_name}_gal.drn",
                             pname="drn_gal", print_input=True,
                             print_flows=True,save_flows=True,
                             boundnames=True, observations=gal_obs)

#Well gallery construction
# i have to check gallery construction to activate wells!! :o
gal_wells_time=np.array([30+9,11,14,15,17,18,20,21])
gal_wells_depth=np.array([25.17,31.04,40.94,42.43,49.04,52.60,57.69,61.30])
gal_wells_depth-=2.5#not 1.5 because gallery size

gal_w_spd=[]
gal_w_spd_tr={}
gal_w=sf.Reader(os.path.join(path_sh,"Pozos_Galeria_Final.shp"))
gal_w1=gal_w.shapeRecords()[0]
firstg_w=gal_w1.shape.__geo_interface__
shp_geomg_w=shape(firstg_w)

for i in range(1,gal_w.numRecords):
    gal_w1=gal_w.shapeRecords()[i]
    firstg_w=gal_w1.shape.__geo_interface__
    shp_geomg_w=shp_geomg_w.union(shape(firstg_w))
    
    
resultg_w=ix.intersect(shp_geomg_w)
# write comments
for i in range(resultg_w.shape[0]):
    for j in range(fondos.shape[0]):
        if fondos[tuple((j,*resultg_w["cellids"][i]))]> dem_Matrix[resultg_w["cellids"][i]] - gal_wells_depth[i]:
            if idom[tuple((0,*resultg_w["cellids"][i]))]==1:#Cambiar J?
                gal_w_spd.append([j,*resultg_w["cellids"][i],fondos[tuple((j,*resultg_w["cellids"][i]))], gal_cond,"gal_w_flow" ])
                print("time= ",gal_wells_time[0:i+1].sum()+time_gal_0)
                print("cell= ",[j,*resultg_w["cellids"][i]], "\n")
                print(gal_w_spd, "\n")
                
    gal_w_spd_tr[int(gal_wells_time[0:i+1].sum()+time_gal_0)]=list(gal_w_spd)
        # break
    
gal_w_obs={"mod_drn_gal_w_obs.csv":[("gal_w-flow", "drn", "gal_w_flow")]}   

drn_gal_w=fp.mf6.ModflowGwfdrn(
    gwf,stress_period_data=gal_w_spd_tr, filename=f"{model_name}_gal_w.drn",
    pname="drn_gal_w", print_input=True,print_flows=True,save_flows=True,
    boundnames=True, observations=gal_w_obs)

# create the initial condition package
start=np.empty((nlay,nrows,ncols))
start[:] = dem_Matrix[0]
ic=fp.mf6.ModflowGwfic(gwf,pname="ic", strt=start)
# ic=fp.mf6.ModflowGwfic(gwf,pname="ic", strt=dem_Matrix)



domin_qd=idom[0]-domin_grav    #this will be necessary for pyemu





# create the output control package
headfile="{}.hds".format(model_name)
head_filerecord =[headfile]
budgetfile = "{}.cbb".format(model_name)
budget_filerecord = [budgetfile]
saverecord = [("HEAD","ALL"),("BUDGET","ALL")]
printrecord = [("HEAD","LAST")]
oc = fp.mf6.ModflowGwfoc(gwf, saverecord=saverecord,
                         head_filerecord=head_filerecord,
                         budget_filerecord=budget_filerecord,
                         printrecord=printrecord)

obslist=[]
import geopandas as gpd
inventory=gpd.read_file("../../05_Vectorial/INV_PAS_V5_DEM.shp")
inv=inventory[inventory["DEPTH_MEA"]>0]
inv.reset_index(drop=True, inplace=True)#because we erased some points

for i in range(inv.shape[0]):
    result_inv=ix.intersect(inv.geometry[i])
    for j in range(fondos.shape[0]):
        print("i= ",str(i)," j= "+str(j))
        print("Z_MEA= ",str(inv["Z_MEA"][i])," dem_lay= "+str(fondos[(0,*result_inv["cellids"][0])]))
        if inv["Z_MEA"][i]>fondos[(j,*result_inv["cellids"][0])]:
            obslist.append([inv["obs_model"][i],"head",(j,*result_inv["cellids"][0])])
            print("yes")
            break
    
    
obsdict={}
obsdict[f"{model_name}.obs.head.csv"]=obslist

obs = fp.mf6.ModflowUtlobs(gwf, print_input= False, continuous=obsdict)

sim.write_simulation()


success,buff=sim.run_simulation()
if not success:
    raise Exception("MODFLOW 6 did not terminate normally.")



# graficos

# no usar, se tira el modelo cuando se rota
# gwf.modelgrid.set_coord_info(angrot=11)

fig=plt.figure()
# gwf.modelgrid.
mapview=fp.plot.PlotMapView(model=gwf)
linecolection = mapview.plot_grid()
quadmesh=mapview.plot_ibound()
quadmesh=mapview.plot_bc("rch", color="purple")
quadmesh=mapview.plot_bc("drn", color="cyan")
quadmesh=mapview.plot_bc("chd", color="blue")
quadmesh=mapview.plot_bc("ghb", color="aquamarine")
quadmesh=mapview.plot_bc("drn_gal", color="brown", plotAll=True,kper=1334+30*3)#
quadmesh=mapview.plot_bc("drn_gal_w", color="olive", plotAll=True,kper=1229+30*3)

linecolection = mapview.plot_grid()

# fig=plt.figure()




