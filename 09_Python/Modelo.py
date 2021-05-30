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
celGlo = 100
celRef = 50

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
import math
def arrayGeneratorRowSmooth(gloRef, locRef, gloSize, locSize):
    print("coordGlob= "+str( gloRef[3]))
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
    print("coordGlob= "+str( gloRef[3]))
    # cloning the list to avoid side effects
    gloRef2=list(gloRef)
    gloRef2[3] =gloRef[3] + delta_coords
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
H=500#borde inferior del modelo
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

fig2=plt.figure(figsize=(12,12))
ax=fig2.add_subplot(1,1,1, aspect="equal")
ax=surface.plot(ax=ax)
plt.colorbar(ax.images[0], shrink=0.7) 
gwf.modelgrid.plot()
# intersecting and resampling raster
dem_Matrix=surface.resample_to_grid(gwf.modelgrid.xcellcenters, gwf.modelgrid.ycellcenters, surface.bands[0], method="nearest")
dis.top=dem_Matrix

botm=np.empty((nlay,nrows,ncols))
print(botm)

botm[0,:,:]=dem_Matrix-100
for i in range(1,nlay):
    botm[i,:,:]=botm[i-1,:,:]-(H/nlay)#Corregir/Revisar/provisional
    
dis.botm=botm

#hydraulic conductivity definition

k_qd_qbg=1e-6*np.ones([nrows,ncols])
k_1="k_qd_qbg.txt"
k_qbo2=1e-7*np.ones([nrows,ncols])
k_2="k_qbo2.txt"
k_qbo1=1e-8*np.ones([nrows,ncols])
k_3="k_qbo1.txt"
k_roc=1e-9*np.ones([nrows,ncols])
k_4="k_roc.txt"

kv_qd_qbg=1e-1*np.ones([nrows,ncols])
kv_1="kv_qd_qbg.txt"
kv_qbo2=1e-1*np.ones([nrows,ncols])
kv_2="kv_qbo2.txt"
kv_qbo1=1e-1*np.ones([nrows,ncols])
kv_3="kv_qbo1.txt"
kv_roc=1e-1*np.ones([nrows,ncols])
kv_4="kv_roc.txt"


# we are including the qbg  by using shapes
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
    domin_grav[tuple(value)]=1 #this will be necessary for pyemu 

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


def dis_layers(path_folder, name_raster, div_layers, kgeo, kvgeo, bottom_model = -1000, min_thick=0):
    """
    dis_layers is made for layer discretization
    this function also helps to define idomain -1 when layer thickness is 0
    raster layers should have a possitive no data value to be corrected by minimun thickness
    
    path_folder is a string with the relative or absolute path of the folder containing raster files
    name_raster is a list of  variables containing the names of raster files, the firxt one has to be model Top
    div_layers is a list of int variables containing the number of division of each layer between raster layers
    it should have the same size as name_raster minus 1(or the same if bottom constant layer is defined[bottom <= -1000])
    bottom is the bottom height  of the model, if negative, the last raster will be used as bottom
    min_thick is a list of int Minimun raster thickness variable allowed for raster layers.
    
    """
    if bottom_model == -1000: 
        botm=np.empty((div_layers.sum(),nrows, ncols))
        demMatrix=np.empty((len(name_raster),nrows, ncols))
        
    else:
        botm=np.empty((div_layers.sum(),nrows, ncols))
        demMatrix=np.empty((len(name_raster) + 1 ,nrows, ncols))
        
    k=div_layers.sum()*[None]
    kv=div_layers.sum()*[None]
    # k=list()
    
    print(demMatrix)
    print(demMatrix.ndim)
    print(demMatrix.shape)
    
    Topo = Raster.load(os.path.join(path_folder, name_raster[0]))
    fig = plt.figure(figsize=(12,12))
    ax=fig.add_subplot(1, 1, 1, aspect = "equal")
    ax = Topo.plot(ax = ax)
    plt.colorbar(ax.images[0], shrink =0.7)
    gwf.modelgrid.plot()
    # intersecting and resampling raster
    demMatrix[0]=Topo.resample_to_grid(gwf.modelgrid.xcellcenters,
                                    gwf.modelgrid.ycellcenters,
                                    Topo.bands[0], method="nearest")
    print(demMatrix)
    print(type(demMatrix))
    print(demMatrix.ndim)
    print(demMatrix.shape)
    # dis.top=demMatrix[0]
    count = 0
    #tengo que colocar K en las capas, y asegurarme de organizar K, idomain abajo

    print("bottom_Shape=",botm.shape)
    thickcells=np.ones(botm.shape)
    for i in range(1, len(capas)+1):
        print("i=",i)
        
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
            demMatrix[i] = bottom.resample_to_grid(gwf.modelgrid.xcellcenters,
                                                  gwf.modelgrid.ycellcenters,
                                                  Topo.bands[0], method="nearest")
        
        
        
        demMatrix[i][demMatrix[i] + min_thick[i-1] > demMatrix[i-1]] = demMatrix[i-1][demMatrix[i] + min_thick[i-1] > demMatrix[i-1]]-min_thick[i-1]
        # demMatrix[i][demMatrix[i]  > demMatrix[i-1]] = demMatrix[i-1][demMatrix[i] > demMatrix[i-1]]
        
        print("count=",count)
        for j in range(count ,div_layers[i-1] + count):
            print("j=",j)
            if div_layers[i-1] == 1:
                botm[j]= demMatrix[i]
                thickcells[j][demMatrix[i] + 0 >= demMatrix[i-1]]=-1
                k[j]=kgeo[i-1]
                kv[j]=kvgeo[i-1]
                break
            
            botm[j] = demMatrix[i-1]+(demMatrix[i]-demMatrix[i-1])*(j-count + 1)/div_layers[i-1]
            thickcells[j][demMatrix[i] + 0>= demMatrix[i-1]]=-1
            k[j]=kgeo[i-1]
            kv[j]=kvgeo[i-1]
        count += div_layers[i-1]
        print("count=",count)
    # if bottom_model != -1000:
    #     botm[-1]= bottom_model
    
    
    
    return botm, demMatrix, thickcells, k, kv
         
raster_names=["R_Topo_Union_Clip.tif", "R_Qbg_Qd.tif", "R_Qbo2.tif", "R_Qbo1.tif"]
capas=np.array([1,2,2,1])
min_thick=([1,0,0,0])
fondos, geol, thickcells, k, kv= dis_layers(path_raster,raster_names, capas, kgeol, kvgeol, bottom_model=500, min_thick=min_thick)
nlay = fondos.shape[0]
dis.nlay = fondos.shape[0]
dis.botm=fondos
        

npf = fp.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, k33overk=True, k33=kv, save_flows=True)
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


#  procedure to include Recharge in shapely

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

rch_spd=[]
for i in range(result.shape[0]):
    rch_spd.append([0,*result["cellids"][i],
                    (0.11/86400)*(result['areas'][i]/
                                  delCArray[result["cellids"][i][0]]/
                                  delRArray[result["cellids"][i][1]])])#hay que revisar si la tupla quedó mal/corregir nombres result, añadir timeseries
# debo revisar la activación de recarga en celdas secas.    
rch=fp.mf6.ModflowGwfrch(gwf, stress_period_data=rch_spd,
                            filename=f"{model_name}.rch",pname="RCH",
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
for i, value in enumerate(domain_cells):
    # print(i)
    # print(value)
    idom[tuple(value)]=1
    
    
# cloning in all layers
idom[:]=idom[0]
idom=idom*thickcells
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
inflow=sf.Reader(os.path.join(path_sh,"Chd_In.shp"))


inflow1=inflow.shapeRecords()[0]
firsti=inflow1.shape.__geo_interface__
shp_geomi=shape(firsti)

for i in range(1,inflow.numRecords):
    inflow1=inflow.shapeRecords()[i]
    firsti=inflow1.shape.__geo_interface__
    shp_geomi=shp_geomi.union(shape(firsti))
    
    
resulti=ix.intersect(shp_geomi)
for i in range(resulti.shape[0]):
    if idom[tuple((0,*resulti["cellids"][i]))]==1:
        ghb_spd.append([0,*resulti["cellids"][i], dem_Matrix[resulti["cellids"][i]]-0,(1e-6*20*5/20)])#problema de la altura -60 que queda debajo de las celdas, hay que calcular conductancia.


ghb=fp.mf6.ModflowGwfghb(gwf,stress_period_data=ghb_spd, filename=f"{model_name}.ghb", pname="ghb", print_input=True,print_flows=True,save_flows=True)



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
quadmesh=mapview.plot_bc("ghb", color="red")
linecolection = mapview.plot_grid()

# fig=plt.figure()





