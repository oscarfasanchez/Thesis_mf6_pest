# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 07:24:36 2021

@author: Oscar
"""

import flopy as fp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from flopy.export import vtk
import pandas as pd

import matplotlib.colors as colors
# workspace="data/modelo_Norte"
workspace="template"
model_name= "modelo_Norte"
sim_name="mfsim.nam"
budget_file = model_name + '.cbb'
head_file = model_name + '.hds'
budget_file = os.path.join(workspace, budget_file)
head_file = os.path.join(workspace, head_file)
#loading model
sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
gwf=sim.get_model(model_name)



def plot_model(
        layer, row, column,
        BC=True, Elv_mdl=False, cr_sect= False, cr_sect_hd=False, path=None, time_sp=0):#, w_s=None):
    """
    

    Parameters
    ----------
    layer : int
        DESCRIPTION.
    row : int
        DESCRIPTION.
    column : int
        DESCRIPTION.
    BC : Boolean, optional
        DESCRIPTION. The default is True.
    Elv_mdl : Boolean, optional
        DESCRIPTION. The default is False.
    cr_sect : Boolean, optional
        DESCRIPTION. The default is False.
    cr_sect_hd : Boolean, optional
        DESCRIPTION. The default is False.
    w_s : str, optional
        DESCRIPTION. The default is None.
    path : str, optional
        folder to be saved image. The default is None.
    time_sp : int, optional
        time stress period to show. The default is 0.

    Returns
    -------
    Cool graphs.

    """
    
    sect_col=column#int(gwf.modelgrid.ncol/2)
    sect_row=row#gwf.modelgrid.nrow-1#int(gwf.modelgrid.nrow/2)
    sec_lay=layer
    # layer=0 #layer to show
    alpha=0.4
    lw=0.3
    arrow_stp=3
    
    # if type(w_s)==str:
    #     workspace=w_s
    #     model_name= "modelo_Norte"
    #     sim_name="mfsim.nam"
    #     budget_file = model_name + '.cbb'
    #     head_file = model_name + '.hds'
    #     budget_file = os.path.join(workspace, budget_file)
    #     head_file = os.path.join(workspace, head_file)
    #     #loading model
    #     sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
    #     gwf=sim.get_model(model_name)
    
    if BC:
        
        #setting graph
        fig=plt.figure(figsize=(15,10))
        
        
        # plotting boundary condition
        # ax=fig.add_subplot(2,1,1, aspect="equal")
        mapview=fp.plot.PlotMapView(gwf)
        
        linecolection = mapview.plot_grid(alpha=alpha, lw=lw)
        quadmesh=mapview.plot_ibound()
        quadmesh=mapview.plot_bc("rch", color="purple")
        quadmesh=mapview.plot_bc("drn", color="cyan")
        quadmesh=mapview.plot_bc("chd", color="blue")
        quadmesh=mapview.plot_bc("ghb", color="aquamarine")
        quadmesh=mapview.plot_bc("drn_gal", color="brown", plotAll=True,kper=1334+90)#
        quadmesh=mapview.plot_bc("drn_gal_w", color="olive", plotAll=True,kper=1229+90)
        # ax.set_title("Plot boundary conditions")
        if path!=None:
            plt.savefig(path+"/BC", dpi=300)
        
    if Elv_mdl:
    
        # plot model bottom elevations
        fig=plt.figure(figsize=(15,10))
        a = gwf.dis.botm.array
        
        # ax = fig.add_subplot(1, 2, 2, aspect='equal')
        # ax.set_title('Model Bottom Elevations')
        mapview = fp.plot.PlotMapView(model=gwf, layer=layer)#change for layer desired
        quadmesh = mapview.plot_array(a)
        inactive = mapview.plot_inactive()
        linecollection = mapview.plot_grid(alpha=alpha)
        cb = plt.colorbar(quadmesh, shrink=0.5)
    
        
        # Contouring Arrays for model elevations
        levels=np.arange(500,1000,10)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        ax.set_title('Model Bottom Elevations')
        mapview = fp.plot.PlotMapView(model=gwf, layer=layer)
        contour_set = mapview.contour_array(a, levels=levels)
        linecollection = mapview.plot_grid(alpha=alpha, lw=lw)
        
        # set up and plot a continuous colorbar in matplotlib for a contour plot
        norm= mpl.colors.Normalize(vmin=contour_set.cvalues.min(), 
                                   vmax=contour_set.cvalues.max())
        sm = plt.cm.ScalarMappable(norm=norm, cmap=contour_set.cmap)
        sm.set_array([])
        fig.colorbar(sm, shrink=0.75);
        
        if path!=None:
            plt.savefig(path+"/Elv", dpi=300)

    if cr_sect or cr_sect_hd:
        
        # Cross Section part
        #Show plan view of cross sections
        # print("numrows= ",gwf.modelgrid.nrow)
        # print("numcols= ",gwf.modelgrid.ncol)
        # print("numcols= ",gwf.modelgrid.nlay)

        
        line = np.array([(gwf.modelgrid.xcellcenters[0,sect_col],
                        gwf.modelgrid.ycellcenters[0,sect_col]),
                        (gwf.modelgrid.xcellcenters[gwf.modelgrid.nrow-1,sect_col],
                        gwf.modelgrid.ycellcenters[gwf.modelgrid.nrow-1,sect_col])])
        
        line_row = np.array([(gwf.modelgrid.xcellcenters[sect_row,0],
                        gwf.modelgrid.ycellcenters[sect_row,0]),
                        (gwf.modelgrid.xcellcenters[sect_row,gwf.modelgrid.ncol-1],
                        gwf.modelgrid.ycellcenters[sect_row,gwf.modelgrid.ncol-1])])
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        ax.set_title(f"Grid layer {sec_lay} BCs (DIS) with cross sectional line")
         # use PlotMapView to plot a DIS model
        mapview = fp.plot.PlotMapView(gwf, layer=sec_lay)
        linecolection = mapview.plot_grid(alpha=alpha, lw=lw)
        quadmesh=mapview.plot_ibound()
        quadmesh=mapview.plot_bc("rch", color="purple")
        quadmesh=mapview.plot_bc("drn", color="cyan")
        quadmesh=mapview.plot_bc("chd", color="blue")
        quadmesh=mapview.plot_bc("ghb", color="aquamarine")
        quadmesh=mapview.plot_bc("drn_gal", color="brown",kper=1334)#
        quadmesh=mapview.plot_bc("drn_gal_w", color="olive",kper=1229)
        lc = plt.plot(line.T[0], line.T[1], 'r--', lw=2)
        linecollection = mapview.plot_grid(alpha=alpha, lw=lw)
        lc2 = plt.plot(line_row.T[0], line_row.T[1], 'b--', lw=2)
        if path!=None:
            plt.savefig(path+f"/plan_{sec_lay}_{sect_row}_{sect_col}_cs", dpi=300)
        
    if cr_sect:
        
        
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        ax.set_title(f"Grid layer {sec_lay} K (DIS) with cross sectional line")
         # use PlotMapView to plot a DIS model
        mapview = fp.plot.PlotMapView(gwf, layer=sec_lay)
        linecolection = mapview.plot_grid(alpha=alpha, lw=lw)
        
        a = gwf.npf.k.array
        csa = mapview.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
        quadmesh=mapview.plot_ibound()
        cb = plt.colorbar(csa, shrink=0.75)
        lc = plt.plot(line.T[0], line.T[1], 'r--', lw=2)
        linecollection = mapview.plot_grid(alpha=alpha, lw=lw)
        lc2 = plt.plot(line_row.T[0], line_row.T[1], 'b--', lw=2)
        if path!=None:
            plt.savefig(path+f"/geo_plan_{sec_lay}_{sect_row}_{sect_col}_cs", dpi=400)

        
        #Column Cross section
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(2, 1, 1)
        # plot boundary condition
        xsect=fp.plot.PlotCrossSection(model=gwf, line={"column":sect_col})
        patches = xsect.plot_bc("rch", color="purple")
        patches = xsect.plot_bc("drn", color="cyan")
        patches = xsect.plot_bc("chd", color="blue")
        patches = xsect.plot_bc("ghb", color="aquamarine")
        patches = xsect.plot_bc("drn_gal", color="brown",kper=1334)#
        patches = xsect.plot_bc("drn_gal_w", color="olive",kper=1229)
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        cb = plt.colorbar(mappable=(patches),shrink=0.75,)
        t = ax.set_title(f'Column {sect_col} Cross-Section with Boundary Conditions')
        # plot xxxx
        ax = fig.add_subplot(2, 1, 2)
        # plot the horizontal hydraulic conductivities
        a = gwf.npf.k.array
        xsect = fp.plot.PlotCrossSection(model=gwf, line={'Column': sect_col})
        csa = xsect.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        t = ax.set_title(f'Column {sect_col} Cross-Section with Horizontal hydraulic conductivity')
        # ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
        cb = plt.colorbar(csa, shrink=0.75)
        if path!=None:
            plt.savefig(path+f"/geo_col_{sect_col}_cs", dpi=400)
        
    

        #Row Cross section
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(2, 1, 1)
        # plot boundary condition
        xsect=fp.plot.PlotCrossSection(model=gwf, line={"row":sect_row})
        patches = xsect.plot_bc("rch", color="purple")
        patches = xsect.plot_bc("drn", color="cyan")
        patches = xsect.plot_bc("chd", color="blue")
        patches = xsect.plot_bc("ghb", color="aquamarine")
        patches = xsect.plot_bc("drn_gal", color="brown",kper=1334)#
        patches = xsect.plot_bc("drn_gal_w", color="olive",kper=1229)
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        cb = plt.colorbar(mappable=(patches),shrink=0.75,)
        t = ax.set_title(f'row {sect_row} Cross-Section with Boundary Conditions')
        # plot xxxx
        ax = fig.add_subplot(2, 1, 2)
        # plot the horizontal hydraulic conductivities
        a = gwf.npf.k.array
        xsect = fp.plot.PlotCrossSection(model=gwf, line={'row': sect_row})
        csa = xsect.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        t = ax.set_title(f'Column {sect_row} Cross-Section with Horizontal hydraulic conductivity')
        # ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
        cb = plt.colorbar(csa, shrink=0.75)
        if path!=None:
            plt.savefig(path+f"/geo_row_{sect_row}_cs", dpi=400)
    if cr_sect_hd:

        #Plotting Specific discharge and head
        # get the specific discharge from the cell budget file
        time=time_sp
        # cbc=fp.utils.CellBudgetFile(budget_file)#cell budget
        # spdis = cbc.get_data(text="SPDIS")[0]
        bud = gwf.output.budget()
        spdis = bud.get_data(text='DATA-SPDIS')[time]
        qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(model=gwf, vectors=spdis)
        
        # get the head from the head file
        
        # head = fp.utils.HeadFile(head_file)
        # hdata = head.get_alldata()[0]
        head = gwf.output.head().get_alldata()
        # plot specific discharge using PlotMapView
        fig = plt.figure(figsize=(8, 8))
        mapview = fp.plot.PlotMapView(model=gwf, layer=layer)
        linecollection = mapview.plot_grid(color="black", alpha=alpha, lw=lw)
        # quadmesh = mapview.plot_array(a=head, alpha=0.5, masked_values=[1e30,500,-1e30])
        quadmesh =mapview.plot_array(head[time],masked_values=[1e30, -1e30])
        quiver = mapview.plot_vector(qx, qy,scale=60, normalize=True,istep=arrow_stp, jstep=arrow_stp, alpha=0.5)
        inactive = mapview.plot_inactive()
        cont=mapview.contour_array(head[time][layer] ,levels=np.linspace(600,1000,10), masked_values=[1e30, -1e30], colors="white")#cmap="seismic")
        plt.title(f"Head and Spcfc Discharge layer {layer} in time {time}")
        plt.colorbar(quadmesh, shrink=0.75 )
        # plt.colorbar(cont, shrink=0.75 )
        plt.clabel(cont,fmt="%1.0f")
        # plt.savefig(fname='pic_head')
        if path!=None:
            plt.savefig(path+f"/head_lay_{layer}_tim_{time}",dpi=400 )
        
        # plotting head and discharge in cross sections
        #Column Cross section
        
        levels=np.arange(500,2400,20)
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(2, 1, 1)
        # plot boundary condition
        xsect=fp.plot.PlotCrossSection(model=gwf, line={"column":sect_col})
        patches = xsect.plot_bc("rch", color="purple")
        patches = xsect.plot_bc("drn", color="cyan")
        patches = xsect.plot_bc("chd", color="blue")
        patches = xsect.plot_bc("ghb", color="aquamarine")
        patches = xsect.plot_bc("drn_gal", color="brown",kper=1334)#Write bug with time.
        patches = xsect.plot_bc("drn_gal_w", color="olive",kper=1229)
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        cb = plt.colorbar(mappable=(patches),shrink=0.75,)
        contour_set = xsect.contour_array(head[time],head=head[time], levels=levels, colors='k', masked_values=[1e30, -1e30])
        plt.clabel(contour_set, fmt='%.1f', colors='k', fontsize=11)
        # cb = plt.colorbar(mappable=contour_set,shrink=0.75,)
        t = ax.set_title(f'Column {sect_col} Cross-Section with BC and contours in time {time}')

        
        # plot xxxx
        ax = fig.add_subplot(2, 1, 2)
        # plot the horizontal hydraulic conductivities
        a = gwf.npf.k.array
        xsect = fp.plot.PlotCrossSection(model=gwf, line={'Column': sect_col})
        csa = xsect.plot_array(head[time],head=head[time] , masked_values=[1e30, -1e30], alpha=0.5 )
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        quiver=xsect.plot_vector(qx, qy, qz,head=head[time],
                                 hstep=arrow_stp, scale= 45,headwidth=2,
                                 headlength=2,headaxislength=2, normalize=True, alpha=0.5)
        # xsect.plot_surface(a=head[0][0], head=head[0], masked_values=[1e30, -1e30])#bug
        t = ax.set_title(f'Column {sect_col} Cross-Section with BC , Heads and flow in time {time}')
        # ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
        cb = plt.colorbar(csa, shrink=0.75)
        if path!=None:
            plt.savefig(path+f"/head_col_{sect_col}_tim_{time}_cs", dpi=400)

        
        #Row Cross section
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(2, 1, 1)
        # plot boundary condition
        xsect=fp.plot.PlotCrossSection(model=gwf, line={"row":sect_row})
        patches = xsect.plot_bc("rch", color="purple")
        patches = xsect.plot_bc("drn", color="cyan")
        patches = xsect.plot_bc("chd", color="blue")
        patches = xsect.plot_bc("ghb", color="aquamarine")
        patches = xsect.plot_bc("drn_gal", color="brown",kper=1334)#
        patches = xsect.plot_bc("drn_gal_w", color="olive",kper=1229)
        contour_set = xsect.contour_array(head[time],head=head[time], levels=levels, colors='k', masked_values=[1e30, -1e30])
        plt.clabel(contour_set, fmt='%.1f', colors='k', fontsize=11)
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        cb = plt.colorbar(mappable=(patches),shrink=0.75,)
        # cb = plt.colorbar(mappable=head[time],shrink=0.75,)
        t = ax.set_title(f'row {sect_row} Cross-Section with BC and contours in time {time}')
        # plot xxxx
        ax = fig.add_subplot(2, 1, 2)
        # plot the horizontal hydraulic conductivities
        xsect = fp.plot.PlotCrossSection(model=gwf, line={'row': sect_row})
        csa = xsect.plot_array(head[time],head=head[0] , masked_values=[1e30, -1e30], alpha=0.5 )
        patches = xsect.plot_ibound()
        linecollection = xsect.plot_grid(alpha=alpha, lw=lw)
        quiver=xsect.plot_vector(qx, qy, qz,head=head[time],
                                 hstep=arrow_stp, scale= 45,headwidth=2,
                                 headlength=2,headaxislength=2, normalize=True, alpha=0.5)
        t = ax.set_title(f'Column {sect_row} Cross-Section with  with Heads and flow in time {time}')
        # ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
        cb = plt.colorbar(csa, shrink=0.75)
        if path!=None:
            plt.savefig(path+f"/head_col_{sect_row}_tim_{time}_cs", dpi=400)
        
        

def paraview_export():
    pass

    """
    Export to paraview
    """
    pv_folder = os.path.join("../" , '11_pv_test')
    if not os.path.exists(pv_folder):
        os.mkdir(pv_folder)
    
    
    gwf.dis.top.export(pv_folder, fmt='vtk')
    
    # 3D Array export
    # export model bottoms
    gwf.dis.botm.export(pv_folder, fmt='vtk')
    
    # transient 2d array
    # export recharge
    gwf.rch.export(pv_folder, fmt='vtk')#doesn't work
    # gwf.rch.rech.export(pv_folder, fmt='vtk')#doesn't work
    # 3D Array export
    # hk export, with points
    gwf.npf.k.export(pv_folder, smooth=True, fmt='vtk', name='HK', point_scalars=True)
    
    # npf export, with points
    gwf.npf.export(pv_folder, smooth=True, fmt='vtk', name='NPF', point_scalars=True)
    
    # DRN export, with points
    gwf.drn_gal.export(pv_folder, fmt='vtk', name='drn_gal_w', point_scalars=False)#doesn't work
    # fp.export.utils.package_export(pv_folder, gwf.drn_gal,fmt='vtk')
    
    # 3D Array export
    # hk export, with points
    gwf.sto.export(pv_folder, smooth=True, fmt='vtk', name='STO', point_scalars=True)
    
    # ghb export, with points
    gwf.ghb_0.export(pv_folder, smooth=True, fmt='vtk', name='ghb_0', point_scalars=True)#doesn't work
    
    # chd export, with points
    gwf.chd.export(pv_folder, smooth=True, fmt='vtk', name='CHD', point_scalars=True)#doesn't work
    
    # model export
    gwf.export(pv_folder, fmt='vtk', binary=True) #works for dis, ic, npf, sto
    
    # head export
    
    heads_output_folder = os.path.join(pv_folder, 'heads_output_test')
    vtk.export_heads(gwf, head_file, heads_output_folder, binary=True, nanval=-1e30)#doesn't work beru well,pvd error, many values nan
    
    
    #with points
    vtk.export_heads(gwf, head_file,
                      heads_output_folder,
                      kstpkper=[(0,0), (0, 49), (0, 99), (0, 999)],#review time steps
                      point_scalars=True, nanval=1e30)#doesn't work many values nan
    
    # Export output cell by cell file to .vtu
    # vtk.export_cbc(gwf, budget_file, pv_folder, kstpkper=[(0, 0), (0, 9), (0, 10), (0, 11)],#review time steps
    #                text=['CONSTANT HEAD', 'STORAGE'], point_scalars=True, binary=True)#doen't work, ERROR
    
    # revisar # https://flopy.readthedocs.io/en/latest/source/flopy.export.vtk.html?highlight=vtk.export_heads#flopy.export.vtk.export_heads
def gal_time_series(path=None):
    df_gal=pd.read_csv(os.path.join(workspace,"mod_drn_gal_obs.csv"))
    df_gal_w=pd.read_csv(os.path.join(workspace,"mod_drn_gal_w_obs.csv"))
    df_gal.replace([1e30,3e30,-1e30],0, inplace=True)
    df_gal_w.replace([1e30,3e30,-1e30],0, inplace=True)
    df_gal_tot=df_gal.copy()
    df_gal_tot["GAL-FLOW"]=(df_gal_tot["GAL-FLOW"]+df_gal_w["GAL_W-FLOW"])*-1000
    df_gal_tot["time"]=df_gal_tot["time"]/86400
    df_gal_tot.plot(x="time", y="GAL-FLOW", kind="line",
                    title="Base realization gallery flow",
                    xlabel="time[d]",ylabel="l/s",
                    figsize=(15,4), grid=True,xlim=[1000,365*4+30*3])

    if path!=None:
            plt.savefig(path+f"/gal_base_time series")
    
    
if __name__ == '__main__':
    nlay=gwf.modelgrid.nlay
    nrow=gwf.modelgrid.nrow
    ncol=gwf.modelgrid.ncol
    print("numlays= ",nlay)
    print("numrows= ",nrow)
    print("numcols= ",ncol)
    layers=np.loadtxt(os.path.join(workspace,"layers"))
    
    
    folder= "v2_post/base"#os.path.join("../06_jpg/", ) 
    full_path=os.path.join("../06_Jpg/",folder)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    
    gal_time_series(path=full_path)
    
    plot_model(int(layers[0:0].sum()), int(nrow*0.7), int(ncol*0.7), BC=True, Elv_mdl=False, cr_sect= True, cr_sect_hd=True, path=full_path, time_sp=0)
    plot_model(int(layers[0:1].sum()), int(nrow*0.7), int(ncol*0.7), BC=False, Elv_mdl=False, cr_sect= False, cr_sect_hd=True, path=full_path, time_sp=1185)
    plot_model(int(layers[0:1].sum()), int(nrow*0.7), int(ncol*0.7), BC=False, Elv_mdl=False, cr_sect= False, cr_sect_hd=True, path=full_path, time_sp=365*4)
    plot_model(int(layers[0:2].sum()), int(nrow*0.7), int(ncol*0.7), BC=False, Elv_mdl=False, cr_sect= False, cr_sect_hd=True, path=full_path, time_sp=1185)
    plot_model(int(layers[0:2].sum()), int(nrow*0.7), int(ncol*0.7), BC=False, Elv_mdl=False, cr_sect= False, cr_sect_hd=True, path=full_path, time_sp=365*4)

    