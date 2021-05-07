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

from flopy.utils import Raster



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
celGlo = 25
celRef = 5

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
    accumCoordinate =  gloRef2[3] - cellArray.sum()
    print("cellarray:", cellArray)
    i=1

    while gloRef2[3] - cellArray.sum() > locRef[3] + celGlo + lenght_afected:
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + celGlo + lenght_afected, gloRef2[3] - cellArray.sum() > locRef[3] + celGlo + lenght_afected)
        cellArray = np.append(cellArray,[gloSize])
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + celGlo + lenght_afected, gloRef2[3] - cellArray.sum() > locRef[3] + celGlo + lenght_afected)
        
    
    i=1    
    while gloRef2[3] - cellArray.sum() > locRef[3] + celGlo:
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + celGlo, gloRef2[3] - cellArray.sum() > locRef[3] + celGlo)
        # print(i," i/cellcom" ,cellcomplete)
        cellArray = np.append(cellArray,smooth_cell[cellcomplete-i])
        # print("smoothCell: "+str(smooth_cell[cellcomplete-i]))
        # print(gloRef2[3] - cellArray.sum()," comparar" ,locRef[3] + celGlo, gloRef2[3] - cellArray.sum() > locRef[3] + celGlo)
        i += 1
        
        
    while gloRef2[3] - cellArray.sum() < locRef[3] + celGlo and gloRef2[3] - cellArray.sum() > locRef[1] - celGlo:
        cellArray = np.append(cellArray,[locSize])
        
    i=0      
    while gloRef2[3] - cellArray.sum() < locRef[1] + celGlo and gloRef2[3] - cellArray.sum() > locRef[1] - lenght_afected:
        cellArray = np.append(cellArray,smooth_cell[i])
        i += 1
        
    while gloRef2[3] - cellArray.sum() < locRef[1] - celGlo and gloRef2[3] - cellArray.sum() > gloRef2[1]:
        cellArray = np.append(cellArray,[gloSize])

    return cellArray



#And DELC is the space between columns, so its the row dimension array
delCArray = arrayGeneratorRowSmooth(GloRefBox, LocRefBox, celGlo, celRef)
print("delCols= \n"+str(delCArray))


#Remember that DELR is space between rows, so it is the column dimension array
delRArray = arrayGeneratorCol(GloRefBox, LocRefBox, celGlo, celRef)
print("delRows= \n"+str(delRArray))



# delCArray = arrayGeneratorRow(GloRefBox, LocRefBox, celGlo, celRef)
# print("delCols= \n"+str(delCArray))

#Calculating number or rows and cols since they are dependant from the discretization
nrows = delCArray.shape[0]
ncols = delRArray.shape[0]
print('Number of rows: %d and number of cols: %d' % (nrows,ncols))

#Define some parameters and values for the spatial and temporal discretization (DIS package)

#Number of layers and layer elevations
nlay = 20
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
path_raster=r"D:\OneDrive - UNIVERSIDAD INDUSTRIAL DE SANTANDER\Maestria\06_Tesis\01_Tesis_Dev\03_Raster\Input_ModelR\Superficies_R_Tiff"

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

botm[0,:,:]=dem_Matrix-100
for i in range(1,nlay):
    botm[i,:,:]=botm[i-1,:,:]-(H/nlay)#Corregir/Revisar/provisional
    
dis.botm=botm

#  procedure to include Recharge in shapely

recarga=sf.Reader(path_sh+"/Zonas_Rec3")
recarga1=recarga.shapeRecords()[0]
first=recarga1.shape.__geo_interface__
# print(first)
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


# creating the main chd inflow boundary condition

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
        chd_spd.append([0,*resulti["cellids"][i], dem_Matrix[resulti["cellids"][i]]-60])


chd=fp.mf6.ModflowGwfchd(gwf,stress_period_data=chd_spd, filename=f"{model_name}.chd", pname="chd", print_input=True,print_flows=True,save_flows=True)



# create the initial condition package
start=np.empty((nlay,nrows,ncols))
start[:] = dem_Matrix[0]
ic=fp.mf6.ModflowGwfic(gwf,pname="ic", strt=start)
# ic=fp.mf6.ModflowGwfic(gwf,pname="ic", strt=dem_Matrix)

k=1e-6
npf = fp.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, save_flows=True)

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
linecolection = mapview.plot_grid()

# fig=plt.figure()





