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
import  shapefile as sf #tomado de github flopy --pip install pyshp
import pandas as pd
import shapely
from shapely.geometry import shape, MultiLineString#lastone from StackO.F.
from flopy.utils.gridintersect import GridIntersect
from flopy.utils import Raster
import geopandas as gpd

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

def import_shapefile(path):
    pyshp_obj=sf.Reader(path)
    shp_obj=pyshp_obj.shapeRecords()[0].shape
    geom_obj= shape(shp_obj.__geo_interface__)
    #TODO look after a better way to do this
    #Intersect shape and grid to get cell info(selection)
    # sel_domain  = set(grid_obj.intersect(shap_domain).cellids)
    return shp_obj, geom_obj

 
    
def def_strd_grid(wd:str, cell_size_big, cell_size_ref, n_domain_shp:str, n_refined_shp:str, change_ratio=1.5, plot=False):# open shapefiles of limits and refinement
    #define structured grid, for dis packages usually
    path_dom=os.path.join(wd, n_domain_shp)    
    pyshp_dom=sf.Reader(path_dom)
    shp_dom=pyshp_dom.shapeRecords()[0].shape
    # geom_dom= shape(shp_dom.__geo_interface__)
    #TODO check if I can imporove this
    path_ref=os.path.join(wd, n_refined_shp)    
    pyshp_ref=sf.Reader(path_ref)
    shp_ref=pyshp_ref.shapeRecords()[0].shape
    # geom_ref= shape(shp_ref.__geo_interface__)
    
    
    dom_array = np.array(shp_dom.points)
    ref_area_array = np.array(shp_ref.points)
    # ref_area_array=np.array([point.shape.points[0] for point in GalleryShp.shapeRecords()])
    
    GloRefBox=np.array(pyshp_dom.bbox)
    LocRefBox=np.array(pyshp_ref.bbox)
    
    if plot:
        # To plot Area model
        fig = plt.figure()
        plt.plot(dom_array[:,0],dom_array[:,1])
        plt.arrow(ref_area_array[0,0],ref_area_array[0,1],-ref_area_array[0,0]+ref_area_array[1,0],-ref_area_array[0,1]+ref_area_array[1,1], shape="full")
        plt.plot([LocRefBox[0],LocRefBox[0],LocRefBox[2],LocRefBox[2],LocRefBox[0]],[LocRefBox[1],LocRefBox[3],LocRefBox[3],LocRefBox[1],LocRefBox[1]])
        plt.plot([GloRefBox[0],GloRefBox[0],GloRefBox[2],GloRefBox[2],GloRefBox[0]],[GloRefBox[1],GloRefBox[3],GloRefBox[3],GloRefBox[1],GloRefBox[1]])
        plt.show()
            
        # if you want to see it on a map
        # crs={'init' : 'epsg:3116'}#deprecated
        # mplleaflet.show(fig, epsg=3116)#muy lento por eso lo comento
    
    print(f'Global coords are: {GloRefBox}')
    print(f'Local coords are: {LocRefBox}')
    
    #Calculating Global Model (Glo) and Local Refinement (Loc) dimensions
    GloLx = GloRefBox[2] - GloRefBox[0] #x_max - x_min
    GloLy = GloRefBox[3] - GloRefBox[1] #y_max - y_min
    print('Global Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (GloLx,GloLy))
    
    LocLx = LocRefBox[2] - LocRefBox[0] #x_max - x_min
    LocLy = LocRefBox[3] - LocRefBox[1] #y_max - y_min
    print('Local Refinement Dimension. Easting Dimension: %8.1f, Northing Dimension: %8.1f' % (LocLx,LocLy))
    
    # y si quiero dejar diagonal la malla?t
    #Defining Global and Local Refinements, for purpose of simplicity cell x and y dimension will be the same
    celGlo = cell_size_big
    celRef = cell_size_ref
    assert celGlo >=celRef, 'you cant refine with a bigger cell size genius'
    if change_ratio==None or 0:
        #And DELC is the space between columns, so its the row dimension array
        delCArray = arrayGeneratorRow(GloRefBox, LocRefBox, celGlo, celRef)
        print("delCols= \n"+str(delCArray))
        
        
        #Remember that DELR is space between rows, so it is the column dimension array
        delRArray = arrayGeneratorCol(GloRefBox, LocRefBox, celGlo, celRef)
        print("delRows= \n"+str(delRArray))
    else:
        #And DELC is the space between columns, so its the row dimension array
        delCArray = arrayGeneratorRowSmooth(GloRefBox, LocRefBox, celGlo, celRef)
        print("delCols= \n"+str(delCArray))
        
        
        #Remember that DELR is space between rows, so it is the column dimension array
        delRArray = arrayGeneratorColSmooth(GloRefBox, LocRefBox, celGlo, celRef)
        print("delRows= \n"+str(delRArray))
        
        
    #Calculating number or rows and cols since they are dependant from the discretization
    nrows = delCArray.shape[0]
    ncols = delRArray.shape[0]
    print('Number of rows: %d and number of cols: %d' % (nrows,ncols))
    return delCArray, delRArray, nrows, ncols, GloRefBox[0], GloRefBox[1]

def setup_dis_layers(path_folder, name_raster, div_layers, kgeo, kvgeo, ssgeo, sygeo, gwf, bottom_model = -1000, min_thick=0):
    """
    dis_layers is made for layer discretization
    this function also helps to define idomain -1 when layer thickness is 0
    also helps to define hydraulic conductivity
    raster layers should have a possitive no data value to be corrected by minimun thickness
    
    path_folder is a string with the relative or absolute path of the folder containing raster files
    name_raster is a list of  variables containing the names of raster files, the firxt one has to be model Top
    div_layers is a list of int variables containing the number of division of each layer between raster layers
    it should have the same size as name_raster minus 1(or the same if bottom constant layer is defined[bottom <= -1000]) ->I need to check this.
    bottom is the bottom height  of the model, if negative, the last raster will be used as bottom
    min_thick is a list of int Minimun raster thickness variable allowed for raster layers.
    
    Return
    botm is a matrix of each layer bottom
    demMatrix is a matrix of geological layers bottoms, with the first value as top
    thickcells is a matrix used for identify no thick cells, later can be used to build an idomain matrix
    k, kv hydraulic condivity matrix for every layer
    
    """
    #TODO solve issue of namelayer not ending in txt?or define if there is an issue
    
    
    nrows = gwf.modelgrid.nrow
    ncols = gwf.modelgrid.ncol
    botm=np.empty((div_layers.sum(),nrows, ncols))
    # deciding if there is an extra layer for model bottom
    if bottom_model == None: 
        
        demMatrix=np.empty((len(name_raster),nrows, ncols))
        
    else:
        
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
    for i in range(1, len(div_layers)+1):#This for will define "geological layers"--I changed capas for div_layers
        print("i=",i)
        # assign geologicalbottom layers to matrix
        if i  == len(div_layers) and bottom_model != None:
            demMatrix[i] = bottom_model
            print("opt1")
        elif i  == len(div_layers) and bottom_model == None:
            print("you probably didn't define a suitable bottom model!!")
            demMatrix[i] = bottom_model#take a look of this, maybe I could change it for something better
            print("opt2")#TODO delete this
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
    # if bottom_model != None#-1000:
    #     botm[-1]= bottom_model
    
    
    
    return botm, demMatrix, thickcells, k, kv, ss, sy    

def setup_dis(delCArray, delRArray, nrows, ncols, xorg, yorg,
              layer_data, model_bottom, mod_hyd_data, shp_dom, gwf,
              path_raster, path_vector,
              ws='.', one_file_per_geolayer=True):
    
    nlay=layer_data['num_lay_per_geo'].sum()#it was one more, mistake? or had a purpose?
    mtop = 1000#dummy parameter just to start seting up
    botm = np.linspace(mtop-((mtop-model_bottom)/nlay), model_bottom, nlay)#dummy botm just to start
    
    # Apply the spatial and temporal discretization parameters to the DIS package
    
    dis = fp.mf6.ModflowGwfdis(gwf, pname= "dis", nlay=nlay,
                                nrow=nrows, ncol=ncols,  delr=delRArray,
                                delc=delCArray, top=mtop, botm=botm,
                                filename= f"{gwf.name}.dis",xorigin=xorg,
                                yorigin=yorg)
    
    #To modify the array of hyd props
    # now we use shapely to instersect
    # plt.figure()
    gwf.modelgrid.plot()
    ix = GridIntersect(gwf.modelgrid, method="structured", rtree=False)
    
    dic_int_cells={}
    for i, j in zip(mod_hyd_data['shp_nam'], mod_hyd_data['mod_geo_lay']):
        _, geom_obj = import_shapefile(os.path.join(path_vector, i))
        dic_int_cells[str(j)] = ix.intersect(geom_obj)

    def modif_mat_zon(matrix, geo_layer, dic_cells_int, key):
        if str(geo_layer) in dic_int_cells.keys():
            for i in dic_cells_int[str(geo_layer)]['cellids']:
                matrix[i]=mod_hyd_data[key][geo_layer]
        return matrix
    

    if one_file_per_geolayer:
        hkey=list(layer_data.keys())[:4]#hydraulic key['kh','kv'... ]
        geol_suff=[x+'.txt' for x in layer_data['nam']]#suffix
        hfiles={}
        for i in hkey:                        
            hfiles[i]=[i+'_'+f for f in geol_suff]#hydraulic file names
            for j, fname in enumerate(hfiles[i]):
                hyd_mat=layer_data[i][j]*np.ones([nrows,ncols])#matrix of values in every geol layer
                #to modify array if there are other geological elements in this layer
                hyd_mat=modif_mat_zon(hyd_mat, j, dic_int_cells, i)
                np.savetxt(os.path.join(ws,fname), hyd_mat)
                       
                
    else:
        hkey=list(layer_data.keys())[:4]#hydraulic key
        # geol_suff=[x+'.txt' for x in layer_data['nam']]#suffix
        hfiles={}
        for i in hkey:                        
            hfiles[i]=layer_data[i]#hydraulic data       
            #TODO I need to finish this to use flopy external files
    
    botm, geol_top_botm, thickcells, k, kv, ss, sy = setup_dis_layers(
        path_raster, layer_data['raster_files'],
        layer_data['num_lay_per_geo'], hfiles['kh'], hfiles['kv'],
        hfiles['ss'], hfiles['sy'], gwf, 
        bottom_model=model_bottom, min_thick=layer_data['min_thick'])
    
    #now correct values according to the new layers      
    dis.top = geol_top_botm[0]
    nlay = botm.shape[0]# is updating something, or I already solved it before
    dis.botm = botm
    #TODO will I use the ? thickcells, k, kv, ss, sy
    #TODO IDOMAIN is missing, I should set up here
    #use IDOMAIN to track geo layers
    _, dom_geom_obj = import_shapefile(os.path.join(path_vector,shp_dom))
    dom_cid=ix.intersect(dom_geom_obj)
    idom=np.zeros((nlay,nrows,ncols))
    for i in dom_cid['cellids']:
        idom[0][i]=1
    
    idom[:]=idom[0]
    idom=idom*thickcells#in case there are pinchout layers
    dis.idomain=idom
    
    return dis, k, kv, ss, sy, ix

def setup_rech(workspace, path_sh,shp_name, rain_file, ix, gwf):
    #  procedure to include Recharge 
    df_rain=pd.read_csv(rain_file, sep=",")#load file    
    df_rain["Fecha"] = pd.to_datetime(df_rain.Fecha, dayfirst=False)#format of of dates setted
    df_rain["time"] = df_rain["Fecha"]-df_rain["Fecha"][0]#days since the first day
    # df_rain["time_s"]=df_rain["time"].astype("timedelta64[s]")#count in seconds for modflow
    df_rain["time_s"] = df_rain["time"].dt.total_seconds()
    #shapely process 
    _, shp_geom=import_shapefile(os.path.join(path_sh, shp_name))#TODO allow more flexibility here
    
    # print(shp_geom)
    print(type(shp_geom))

    delCArray=gwf.modelgrid.delc
    delRArray=gwf.modelgrid.delr

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
                                filename=f"{gwf.name}.rch",pname="RCH",
                                auxiliary="rain_mult",
                                auxmultname="rain_mult",
                                timeseries= ts_dict,
                                print_input=False,print_flows=False,save_flows=True)
    
    return rch

def setup_drn_creeks(workspace, path_sh, creek_file, ix, gwf):
    # creating drains from shapefiles

    queb=sf.Reader(os.path.join(path_sh,"Queb_corr.shp"))


    queb1=queb.shapeRecords()[0]
    firstq=queb1.shape.__geo_interface__
    shp_geomq=shape(firstq)
    # list_shp =[]  
    for i in range(1,queb.numRecords):
        queb1=queb.shapeRecords()[i]
        firstq=queb1.shape.__geo_interface__
        shp_geomq=shp_geomq.union(shape(firstq))
        # shp_geomq=shapely.ops.unary_union(shape(firstq))
        # list_shp.append(shape(firstq))
    # shp_geomq=MultiLineString(lines=list_shp)
    # shp_geomq=shapely.ops.unary_union(shp_geomq)
    # shp_geomq=shapely.ops.unary_union(shp_geomq)
        #change conductance accordingly to area/lenght?
    resultq=ix.intersect(shp_geomq)
    drn_spd=[]
    idom=gwf.modelgrid.idomain
    delCArray=gwf.modelgrid.delc
    delRArray=gwf.modelgrid.delr
    dem_Matrix=gwf.modelgrid.top
    for i in range(resultq.shape[0]):
        if idom[tuple((0,*resultq["cellids"][i]))]==1:
            drn_spd.append([0,*resultq["cellids"][i], dem_Matrix[resultq["cellids"][i]]+1, 1e-6*delCArray[resultq["cellids"][i][0]]*delRArray[resultq["cellids"][i][1]]])#falta agregar valores de quebradas, I need terrain

    drn=fp.mf6.ModflowGwfdrn(gwf,stress_period_data=drn_spd, filename=f"{gwf.name}.drn",
                              pname="drn", print_input=False,print_flows=False,save_flows=True)

    return drn

def setup_chd_rivers(workspace, path_sh, river_file, ix, gwf):#TODO check if works
    # # cretaing Rivers from shapefiles as Constant heads

    rio=sf.Reader(os.path.join(path_sh,river_file))


    rio1=rio.shapeRecords()[0]
    firstr=rio1.shape.__geo_interface__
    shp_geomr=shape(firstr)

    for i in range(1,rio.numRecords):
        rio1=rio.shapeRecords()[i]
        firstr=rio1.shape.__geo_interface__
        shp_geomr=shp_geomr.union(shape(firstr))
    
    idom=gwf.modelgrid.idomain    
    dem_Matrix=gwf.modelgrid.top    
    resultr=ix.intersect(shp_geomr)
    chd_spd=[]
    for i in range(resultr.shape[0]):
        if idom[tuple((0,*resultr["cellids"][i]))]==1:
            chd_spd.append([0,*resultr["cellids"][i], dem_Matrix[resultr["cellids"][i]]+1 ])#falta agregar valores de quebradas, I need terrain

    chd_riv=fp.mf6.ModflowGwfchd(gwf,stress_period_data=chd_spd, filename=f"{gwf.name}.chd", pname="chd", print_input=False,print_flows=True,save_flows=True)
    return chd_riv

def setup_ghb_lat_inflow(workspace, path_sh,shp_name, rain_file, ix, gwf, k_geo, capas):
    # creating the main ghc inflow boundary condition

    ghb_spd=[]
    ghb_spd_tr=[]
    # defining conductances per thickness unit in general head boundary condition
    celGlo=gwf.modelgrid.delr.max()
    fondos=gwf.modelgrid.botm
    dem_Matrix=gwf.modelgrid.top 
    idom=gwf.modelgrid.idomain
    #this is repeated, could be better a different approach
    df_rain=pd.read_csv(rain_file, sep=",")#load file    
    df_rain["Fecha"]=pd.to_datetime(df_rain.Fecha, dayfirst=False)#format of of dates setted
    df_rain["time"]=df_rain["Fecha"]-df_rain["Fecha"][0]#days since the first day
    # df_rain["time_s"]=df_rain["time"].astype("timedelta64[s]")#count in seconds for modflow
    df_rain["time_s"] = df_rain["time"].dt.total_seconds()
    
    
    #TODO changing this approach for a k cell based one?
    c_1=k_geo[0]*(celGlo/2000)#*thickness_layer to get conductance

    c_2=k_geo[1]*(celGlo/2000)

    c_3=k_geo[2]*(celGlo/2000)

    c_4=k_geo[3]*(celGlo/2000)

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
    # df_rain["time_s"]=df_rain["time"].astype("timedelta64[s]")#count in seconds for modflow
    df_rain["time_s"] = df_rain["time"].dt.total_seconds()
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
    ghb_i=[]
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
        
        ghb=fp.mf6.ModflowGwfghb(gwf,stress_period_data=ghb_spd_txt,
                              filename=f"{gwf.name}_{i}.ghb", pname=f"ghb_{i}",
                              # auxiliary="rain_mult",
                              # auxmultname="rain_mult",
                              # timeseries =ts_dict,
                              print_input=False,print_flows=False,save_flows=True)
        print(i)
        ghb_i.append(ghb)
    return ghb_i
        
        
def setup_gallery_drn(workspace, path_sh, gal_file, ix, gwf):
    # # gallery construction
    fondos=gwf.modelgrid.botm
    idom=gwf.modelgrid.idomain
    
    gal = sf.Reader(os.path.join(path_sh, gal_file))
    gal1 = gal.shapeRecords()[0]
    firstg = gal1.shape.__geo_interface__
    shp_geomg = shape(firstg)

    for i in range(1,gal.numRecords):
        gal1=gal.shapeRecords()[i]
        firstg=gal1.shape.__geo_interface__
        shp_geomg=shp_geomg.union(shape(firstg))
        
        
    resultg=ix.intersect(shp_geomg)
    

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
    gal_spd=[]
    gal_spd_tr={}
    # write comments
    for i in range(resultg.shape[0]):#iterate gallery cells
        for j in range(fondos.shape[0]):#iterate over layer numbers
            if height_0+resultg["lengths"][0:i].sum()*gal_slp > fondos[tuple((j,*resultg["cellids"][i]))]:
                if idom[tuple((0,*resultg["cellids"][i]))]==1:#Cambiar J?
                    gal_spd.append([j,*resultg["cellids"][i],
                                    height_0+resultg["lengths"][0:i].sum()*gal_slp, gal_cond*resultg["lengths"][i],"gal_flow"])
                    print('gal time = ',int(resultg["lengths"][0:i+1].sum()//gal_speed_exc)+1+time_gal_0)
                    gal_spd_tr[int(resultg["lengths"][0:i+1].sum()//gal_speed_exc)+1+time_gal_0]=list(gal_spd)#review indent
                    break
    gal_spd_tr_2 = {}
    for i in range(1, len(gal_spd_tr)):#to set up a gradual increasing of conductance for a smooth activation of the BC
        strt = list(gal_spd_tr.keys())[i-1]
        end = list(gal_spd_tr.keys())[i]
        for j in range(strt + 1, end): #TODO check I wrote this sleepy
            time = end - strt
            count = j - strt
            # h_delta = (gal_spd_tr[end][-1][3] - gal_spd_tr[strt][-1][3])/time
            h = gal_spd_tr[end][-1][3]#gal_spd_tr[strt][-1][3] + h_delta*count
            cond_delta = gal_spd_tr[end][-1][4]/time
            cond =cond_delta*count
            gal_spd_tr_2[j] = gal_spd_tr[strt].copy()
            gal_spd_tr_2[j].append([*gal_spd_tr[end][-1][:3], h, cond, 'gal_flow'])#TODO add previous cells
        
    gal_spd_tr_2.update(gal_spd_tr)
    gal_spd_tr = {x:gal_spd_tr_2[x] for x in range(min(gal_spd_tr_2.keys()),max(gal_spd_tr_2.keys()))}     
    gal_obs={"mod_drn_gal_obs.csv":[("gal-flow", "drn", "gal_flow")]}        

    drn_gal=fp.mf6.ModflowGwfdrn(gwf,stress_period_data=gal_spd_tr,
                                  filename=f"{gwf.name}_gal.drn",
                                  pname="drn_gal", print_input=False,
                                  print_flows=True,save_flows=True,
                                  boundnames=True, observations=gal_obs)
    return drn_gal

def setup_gal_wel(workspace, path_sh, gal_w_file, ix, gwf):
    #Well gallery construction
    fondos=gwf.modelgrid.botm
    idom=gwf.modelgrid.idomain
    dem_Matrix=gwf.modelgrid.top 
    nlay = gwf.modelgrid.nlay
    nrows = gwf.modelgrid.nrow
    ncols = gwf.modelgrid.ncol
    
    # i have to check gallery construction to activate wells!! :o
    gal_wells_time=np.array([30+9,11,14,15,17,18,20,21])
    gal_wells_depth=np.array([25.17,31.04,40.94,42.43,49.04,52.60,57.69,61.30])
    gal_wells_depth-=2.5#not 1.5 because gallery size
    #TODO time_gal_0 is repeated, and the following variables as well
    time_gal_0 = 365*3+30*3#time when gallery construction begins
    gal_perim = 9.2106+1.3906*2+0.46#ignoring depth channel
    gal_thick = 0.3
    gal_k_conc = 1e-7
    gal_k_grav = 1e-7
    gal_drn_space = 10#space between sections with drains
    len_drn = 1.5*7#outside gallery in 2 inches, there are seven drains
    gal_drn_cnd = gal_k_grav*len_drn*np.pi*2*(2*0.0254)/0.01#K*L*With/thickness
    gal_cond = (gal_perim*gal_k_conc/gal_thick)+gal_drn_cnd/gal_drn_space#
    

    gal_w_spd=[]
    gal_w_spd_tr={}
    gal_w=sf.Reader(os.path.join(path_sh, gal_w_file))
    gal_w1=gal_w.shapeRecords()[0]
    firstg_w=gal_w1.shape.__geo_interface__
    shp_geomg_w=shape(firstg_w)

    for i in range(1,gal_w.numRecords):
        gal_w1=gal_w.shapeRecords()[i]
        firstg_w=gal_w1.shape.__geo_interface__
        shp_geomg_w=shp_geomg_w.union(shape(firstg_w))
        
        
    resultg_w=ix.intersect(shp_geomg_w)
    # write comments
    for i in range(resultg_w.shape[0]):#iterate every well
        print('i=',i, "\n")
        for j in range(fondos.shape[0]):#iterate every layer bottom
            print('j=', j, "\n")
            if fondos[tuple((j,*resultg_w["cellids"][i]))]> dem_Matrix[resultg_w["cellids"][i]] - gal_wells_depth[i]:
                if idom[tuple((0,*resultg_w["cellids"][i]))]==1:#Cambiar J?
                    gal_w_spd.append([j,*resultg_w["cellids"][i],fondos[tuple((j,*resultg_w["cellids"][i]))], gal_cond,"gal_w_flow" ])
                    print('i=',i,', j=', j, "\n")
                    print("time= ",gal_wells_time[0:i+1].sum()+time_gal_0)
                    print("cell= ",[j,*resultg_w["cellids"][i]], "\n")
                    print(gal_w_spd, "\n")
                    
        gal_w_spd_tr[int(gal_wells_time[0:i+1].sum()+time_gal_0)]=list(gal_w_spd)
            # break
        
    gal_w_obs={"mod_drn_gal_w_obs.csv":[("gal_w-flow", "drn", "gal_w_flow")]}   

    drn_gal_w=fp.mf6.ModflowGwfdrn(
        gwf,stress_period_data=gal_w_spd_tr, filename=f"{gwf.name}_gal_w.drn",
        pname="drn_gal_w", print_input=False,print_flows=True,save_flows=True,
        boundnames=True, observations=gal_w_obs)

    return drn_gal_w

def setup_obs(workspace, path_sh, obs_h_file, ix, gwf):
    fondos=gwf.modelgrid.botm
    obslist=[]
    
    inventory=gpd.read_file(os.path.join(path_sh,'..' ,obs_h_file))
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
    obsdict[f"{gwf.name}.obs.head.csv"]=obslist

    obs = fp.mf6.ModflowUtlobs(gwf, print_input= False, continuous=obsdict)
    return obs
         
def setup_mf6(workspace, model_name,
              nsper,cell_size_big,cell_size_ref,
              hbot=500, exe_name='mf6',
              path_sh= os.path.join("..","02_Shp_Vect","Input_Model"),
              path_raster=os.path.join("..", "03_Raster", "Input_ModelR", "Superficies_R_Tiff")):
#creating general simulation
    
    sim =fp.mf6.MFSimulation(sim_name=model_name, version="mf6",
                              exe_name=exe_name,
                              sim_ws=workspace)
    
    #setting modflow time
    # [perlen, nstp, tsmult]
    # 2/21/2021 real rain ideam until there
    
    
    #there are 86400 seconds
    time_disc=[(86400, 1, 1.0) for _ in range(nsper-1)]#[(1,1,1.0)]
    time_disc.insert(0,(1,1,1.0))#inserting the steady stress period at the beginning of list
    tdis= fp.mf6.ModflowTdis(sim, pname="tdis",
                              time_units="SECONDS", 
                              nper=nsper, perioddata=time_disc)
    period_ats =[(i, 86400, 1.0e-5, 86400, 2.0, 5.0) for i in range(1, nsper)]
    ats = fp.mf6.ModflowUtlats(tdis,maxats=len(period_ats), perioddata= period_ats, pname="ats"   )
    # Create gwf model
    model_nam_file=f"{model_name}.nam"
    gwf = fp.mf6.ModflowGwf(sim, modelname=model_name,
                            model_nam_file=model_nam_file,
                            newtonoptions="UNDER_RELAXATION", save_flows=True,)
    
    
    # Setting the solver
    ims=  fp.mf6.modflow.mfims.ModflowIms(sim, pname="ims",
                                          complexity= "MODERATE",print_option="SUMMARY",
                                          outer_maximum=350, outer_dvclose=0.1,
                                          under_relaxation="DBD", under_relaxation_gamma=0.1,
                                          under_relaxation_theta=0.7, under_relaxation_kappa=0.1,
                                          backtracking_number=20, backtracking_tolerance=20,
                                          backtracking_reduction_factor=0.1, backtracking_residual_limit=0.002, 
                                          inner_maximum=500, inner_dvclose=0.01,
                                          linear_acceleration="bicgstab", preconditioner_levels=20,
                                          preconditioner_drop_tolerance=0.0001,ats_outer_maximum_fraction=0.05,
                                          scaling_method="DIAGONAL", reordering_method="MD",
                                          number_orthogonalizations=2)
    
    delCArray, delRArray, nrows, ncols, xorg, yorg = def_strd_grid(path_sh,
                                                      cell_size_big, cell_size_ref, 'Domin_Mod.shp',
                                                      'Alineamiento_Galeria.shp',
                                                      1.5, True)
    mdl_dom_namfile='Domin_Mod'#TODO define model domain, in every layer?
    #define the file names of every geological layer
    raster_names = ["R_Topo_resamp.tif", "R_Qbg_Qd.tif", "R_Qbo2.tif", "R_Qbo1.tif"]
    geo_layer_names = ['qd_qbg','qbo2','qbo1','roc']
    # define other data about layering:
    capas = np.array([1,2,2,1])
    min_thick = ([5,5,5,5])    
    #define some prior mean values 
    k_geo_values = [1e-7, 1e-7, 1e-8, 1e-7]
    kv_geo_values = [1e-1, 1e-1, 1e-1, 1e-1]
    sy_geo_values = [1e-1, 1e-1, 1e-1, 1e-1]
    ss_geo_values = [1e-4, 1e-4, 1e-4, 1e-5]
    
    dic_geo_layer={'kh':k_geo_values,'kv':kv_geo_values, 'sy':sy_geo_values, 'ss':k_geo_values,
                   'nam':geo_layer_names, 'raster_files':raster_names,
                   'min_thick':min_thick, 'num_lay_per_geo': capas}
    #TODO check the definition of external .txt hydraulic property files, and the missing of .txt
    #setting up a modification of x layer  hydraulic properties not matching the layers in dic_geo_layer

    k_geo_mod = [1e-5]
    kv_geo_mod = [1e-1]
    sy_geo_mod = [1e-1]
    ss_geo_mod = [1e-4]
    geo_lay_mod =   [0] 
    shp_names = [os.path.join('Superficies','Cont_Qbg2.shp')]
    
    dic_mod_geo_layer={'kh':k_geo_mod,'kv':kv_geo_mod, 'sy':sy_geo_mod, 'ss':k_geo_mod,
                       'mod_geo_lay':geo_lay_mod, 'shp_nam':shp_names}
    
    
    #TODO define a plot variable
    dis, kh, kv, ss, sy, grd_int = setup_dis(delCArray, delRArray, nrows, ncols,
                                             xorg, yorg, dic_geo_layer, hbot,
                                             dic_mod_geo_layer, mdl_dom_namfile, gwf,
                                             path_raster, path_sh, ws=workspace,
                                             one_file_per_geolayer=True)#TODO setup one file per geo layer false

    #TODO need to setup Gravoso properties


    
    # define node property flow package
    npf = fp.mf6.ModflowGwfnpf(gwf, icelltype=1, k=kh, k33overk=True, k33=kv,
                                save_flows=True, save_specific_discharge = True)
    # specifies storage
    sto = fp.mf6.ModflowGwfsto(gwf, pname="sto", save_flows=True, iconvert=1, ss=ss, sy=sy, steady_state= {0:True}, transient={1:True})
    #later i will need it to pyemu
    np.savetxt(os.path.join(workspace,"layers"), capas)
    
    # #start Setting up Boundary conditions
    path_rain_hist=os.path.join("..","04_Xls","Lluvia_Ideam.csv")
    rain_area_nam_shp="Zonas_Rec3"
    rch=setup_rech(workspace, path_sh,rain_area_nam_shp, path_rain_hist, grd_int, gwf)
    
    path_creeks_nam_shp="Queb_corr.shp"
    drn=setup_drn_creeks(workspace, path_sh, path_creeks_nam_shp, grd_int, gwf)


    path_rivers_nam_shp="Rios.shp"
    chd_riv=setup_chd_rivers(workspace, path_sh, path_rivers_nam_shp, grd_int, gwf)
        
    k_ghb=k_geo_values.copy()
    k_ghb[0]= k_geo_mod[0]
    path_bound_ghb_nam_shp="Chd_In.shp"
    ghb_i = setup_ghb_lat_inflow(workspace, path_sh, path_bound_ghb_nam_shp,
                             path_rain_hist, grd_int, gwf, k_ghb, capas)
    
    path_gallery_nam_shp="Alineamiento_Galeria.shp"
    drn_gal = setup_gallery_drn(workspace, path_sh, path_gallery_nam_shp, grd_int, gwf)
    
    path_gal_wel_nam_shp = "Pozos_Galeria_Final.shp"
    drn_gal_w = setup_gal_wel(workspace, path_sh, path_gal_wel_nam_shp, grd_int, gwf)

    # create the initial condition package
    start=np.empty((gwf.modelgrid.nlay, gwf.modelgrid.nrow, gwf.modelgrid.ncol))
    start[:] = gwf.modelgrid.top 
    ic=fp.mf6.ModflowGwfic(gwf, pname="ic", strt=start)
    # ic=fp.mf6.ModflowGwfic(gwf,pname="ic", strt=dem_Matrix)
    
    head_obs_nam_shp = "INV_PAS_V5_DEM.shp"
    obs = setup_obs(workspace, path_sh, head_obs_nam_shp, grd_int, gwf)

    #TODO CHECK THIS
    # domin_qd=idom[0]-domin_grav    #this will be necessary for pyemu


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
    # success,buff=sim.run_simulation()
    # if not success:
    #     raise Exception("MODFLOW 6 did not terminate normally.")
    return sim

    
def plot_model(sim, model_name):
    gwf=sim.get_model(model_name)
    # graficos
    
    # no usar, se tira el modelo cuando se rota
    # gwf.modelgrid.set_coord_info(angrot=11)
    
    fig=plt.figure()
    # gwf.modelgrid.
    mapview=fp.plot.PlotMapView(model=gwf)
    # linecolection = mapview.plot_grid()
    quadmesh=mapview.plot_ibound()
    quadmesh=mapview.plot_bc("rch", color="purple")
    quadmesh=mapview.plot_bc("drn", color="cyan")
    quadmesh=mapview.plot_bc("chd", color="blue")
    quadmesh=mapview.plot_bc("ghb", color="aquamarine")
    quadmesh=mapview.plot_bc("drn_gal", color="brown", plotAll=True,kper=1334+30*3)#
    quadmesh=mapview.plot_bc("drn_gal_w", color="olive", plotAll=True,kper=1229+30*3)
    
    # linecolection = mapview.plot_grid()
    
    # fig=plt.figure()
def main():
    model_name= "modelo_Norte"
    
    # Creation of workspace folder
    workspace = os.path.join("data", model_name)
    if not os.path.exists(workspace):
        os.makedirs(workspace)
     
    exe_name=os.path.join('..',"10_Exe","mf6.exe")
    nsper=365*4#int(365*4.5)#number of  stress periods, 365 per year
    cell_size_big=80
    cell_size_ref=10
    sim = setup_mf6(workspace, model_name, nsper,
                    cell_size_big=cell_size_big,cell_size_ref=cell_size_ref,
                    exe_name=exe_name)
    success,buff=sim.run_simulation()
    if not success:
        raise Exception("MODFLOW 6 did not terminate normally.")
    plot_model(sim, model_name)
    
if __name__ == "__main__":
    main()
    # path_dom= os.path.join("..","02_Shp_Vect","Input_Model", "Domin_Mod.shp").
    path_sh= os.path.join("..","02_Shp_Vect","Input_Model")
    # geom_dom, shape_dom = import_shapefile(path_dom)
    # delCArray, delRArray, nrows, ncols = def_dis_grid(path_sh,
    #                                                   40, 3, 'Domin_Mod.shp',
    #                                                   'Alineamiento_Galeria.shp',
    #                                                   1.5, True)


