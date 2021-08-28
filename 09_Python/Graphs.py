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

import matplotlib.colors as colors

workspace="data/modelo_Norte"
model_name= "modelo_Norte"
sim_name="mfsim.nam"
budget_file = model_name + '.cbb'
head_file = model_name + '.hds'
budget_file = os.path.join(workspace, budget_file)
head_file = os.path.join(workspace, head_file)
#loading model
sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
gwf=sim.get_model(model_name)

layer=0 #layer to show

#setting graph
fig=plt.figure(figsize=(15,10))


# plotting boundary condition
# ax=fig.add_subplot(2,1,1, aspect="equal")
mapview=fp.plot.PlotMapView(gwf)

linecolection = mapview.plot_grid()
quadmesh=mapview.plot_ibound()
quadmesh=mapview.plot_bc("rch", color="purple")
quadmesh=mapview.plot_bc("drn", color="cyan")
quadmesh=mapview.plot_bc("chd", color="blue")
quadmesh=mapview.plot_bc("ghb", color="aquamarine")
quadmesh=mapview.plot_bc("drn_gal", color="brown", plotAll=True,kper=1334)#
quadmesh=mapview.plot_bc("drn_gal_w", color="olive", plotAll=True,kper=1229)
# ax.set_title("Plot boundary conditions")

# plot model bottom elevations
fig=plt.figure(figsize=(15,10))
a = gwf.dis.botm.array

# ax = fig.add_subplot(1, 2, 2, aspect='equal')
# ax.set_title('Model Bottom Elevations')
mapview = fp.plot.PlotMapView(model=gwf, layer=layer)#change for layer desired
quadmesh = mapview.plot_array(a)
inactive = mapview.plot_inactive()
linecollection = mapview.plot_grid()
cb = plt.colorbar(quadmesh, shrink=0.5)

# Contouring Arrays for model elevations
levels=np.arange(500,1000,10)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_title('Model Bottom Elevations')
mapview = fp.plot.PlotMapView(model=gwf, layer=layer)
contour_set = mapview.contour_array(a, levels=levels)
linecollection = mapview.plot_grid()

# set up and plot a continuous colorbar in matplotlib for a contour plot
norm= mpl.colors.Normalize(vmin=contour_set.cvalues.min(), 
                           vmax=contour_set.cvalues.max())
sm = plt.cm.ScalarMappable(norm=norm, cmap=contour_set.cmap)
sm.set_array([])
fig.colorbar(sm, shrink=0.75);


#Plotting Specific discharge and head
# get the specific discharge from the cell budget file

# cbc=fp.utils.CellBudgetFile(budget_file)#cell budget
# spdis = cbc.get_data(text="SPDIS")[0]
bud = gwf.output.budget()
spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(model=gwf, vectors=spdis)

# get the head from the head file

# head = fp.utils.HeadFile(head_file)
# hdata = head.get_alldata()[0]
head = gwf.output.head().get_alldata()
# plot specific discharge using PlotMapView
fig = plt.figure(figsize=(8, 8))
mapview = fp.plot.PlotMapView(model=gwf, layer=layer)
linecollection = mapview.plot_grid(color="black")
# quadmesh = mapview.plot_array(a=head, alpha=0.5, masked_values=[1e30,500,-1e30])
quadmesh =mapview.plot_array(head[0],masked_values=[1e30, -1e30])
quiver = mapview.plot_vector(qx, qy, normalize=True,)
inactive = mapview.plot_inactive()
cont=mapview.contour_array(head[0][layer] ,levels=np.linspace(600,1000,10), masked_values=[1e30, -1e30], colors="white")#cmap="seismic")
plt.title("Specific Discharge (" + r'$L/T$' + ') layer {} in time 0'.format(layer))
plt.colorbar(quadmesh, shrink=0.75 )
# plt.colorbar(cont, shrink=0.75 )
plt.clabel(cont,fmt="%1.0f")
# plt.savefig(fname='pic_head')


# Cross Section part
#Show plan view of cross sections
print("numrows= ",gwf.modelgrid.nrow)
print("numcols= ",gwf.modelgrid.ncol)
print("numcols= ",gwf.modelgrid.nlay)
sect_col=12
sect_row=10
sec_lay=0

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
linecolection = mapview.plot_grid()
quadmesh=mapview.plot_ibound()
quadmesh=mapview.plot_bc("rch", color="purple")
quadmesh=mapview.plot_bc("drn", color="cyan")
quadmesh=mapview.plot_bc("chd", color="blue")
quadmesh=mapview.plot_bc("ghb", color="aquamarine")
quadmesh=mapview.plot_bc("drn_gal", color="brown",kper=1334)#
quadmesh=mapview.plot_bc("drn_gal_w", color="olive",kper=1229)
lc = plt.plot(line.T[0], line.T[1], 'r--', lw=2)
linecollection = mapview.plot_grid()
lc2 = plt.plot(line_row.T[0], line_row.T[1], 'b--', lw=2)



fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_title(f"Grid layer {sec_lay} K (DIS) with cross sectional line")
 # use PlotMapView to plot a DIS model
mapview = fp.plot.PlotMapView(gwf, layer=sec_lay)
linecolection = mapview.plot_grid()

a = gwf.npf.k.array
csa = mapview.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
quadmesh=mapview.plot_ibound()
cb = plt.colorbar(csa, shrink=0.75)
lc = plt.plot(line.T[0], line.T[1], 'r--', lw=2)
linecollection = mapview.plot_grid()
lc2 = plt.plot(line_row.T[0], line_row.T[1], 'b--', lw=2)


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
linecollection = xsect.plot_grid()
cb = plt.colorbar(mappable=(patches),shrink=0.75,)
t = ax.set_title(f'Column {sect_col} Cross-Section with Boundary Conditions')
# plot xxxx
ax = fig.add_subplot(2, 1, 2)
# plot the horizontal hydraulic conductivities
a = gwf.npf.k.array
xsect = fp.plot.PlotCrossSection(model=gwf, line={'Column': sect_col})
csa = xsect.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
patches = xsect.plot_ibound()
linecollection = xsect.plot_grid()
t = ax.set_title(f'Column {sect_col} Cross-Section with Horizontal hydraulic conductivity')
# ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
cb = plt.colorbar(csa, shrink=0.75)


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
linecollection = xsect.plot_grid()
cb = plt.colorbar(mappable=(patches),shrink=0.75,)
t = ax.set_title(f'row {sect_row} Cross-Section with Boundary Conditions')
# plot xxxx
ax = fig.add_subplot(2, 1, 2)
# plot the horizontal hydraulic conductivities
a = gwf.npf.k.array
xsect = fp.plot.PlotCrossSection(model=gwf, line={'row': sect_row})
csa = xsect.plot_array(a, norm=colors.LogNorm(vmin=a.min(), vmax=a.max()) )
patches = xsect.plot_ibound()
linecollection = xsect.plot_grid()
t = ax.set_title(f'Column {sect_row} Cross-Section with Horizontal hydraulic conductivity')
# ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
cb = plt.colorbar(csa, shrink=0.75)


# plotting head and discharge in cross sections
#Column Cross section
time=1400
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
linecollection = xsect.plot_grid()
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
csa = xsect.plot_array(head[time],head=head[0] , masked_values=[1e30, -1e30], alpha=0.5 )
patches = xsect.plot_ibound()
linecollection = xsect.plot_grid()
quiver=xsect.plot_vector(qx, qy, qz,head=head[time],
                         hstep=2, scale= 30,headwidth=2,
                         headlength=2,headaxislength=2, normalize=True, alpha=0.5)
# xsect.plot_surface(a=head[0][0], head=head[0], masked_values=[1e30, -1e30])#bug
t = ax.set_title(f'Column {sect_col} Cross-Section with  with Heads and flow in time {time}')
# ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
cb = plt.colorbar(csa, shrink=0.75)

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
linecollection = xsect.plot_grid()
cb = plt.colorbar(mappable=(patches),shrink=0.75,)
# cb = plt.colorbar(mappable=head[time],shrink=0.75,)
t = ax.set_title(f'row {sect_row} Cross-Section with BC and contours in time {time}')
# plot xxxx
ax = fig.add_subplot(2, 1, 2)
# plot the horizontal hydraulic conductivities
xsect = fp.plot.PlotCrossSection(model=gwf, line={'row': sect_row})
csa = xsect.plot_array(head[time],head=head[0] , masked_values=[1e30, -1e30], alpha=0.5 )
patches = xsect.plot_ibound()
linecollection = xsect.plot_grid()
quiver=xsect.plot_vector(qx, qy, qz,head=head[time],
                         hstep=2, scale= 30,headwidth=2,
                         headlength=2,headaxislength=2, normalize=True, alpha=0.5)
t = ax.set_title(f'Column {sect_row} Cross-Section with  with Heads and flow in time {time}')
# ax.pcolor(norm=colors.LogNorm(vmin=a.min(), vmax=a.max()))
cb = plt.colorbar(csa, shrink=0.75)


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
gwf.rch.export(pv_folder, fmt='vtk')#doen't work

# 3D Array export
# hk export, with points
gwf.npf.k.export(pv_folder, smooth=True, fmt='vtk', name='HK', point_scalars=True)

# npf export, with points
gwf.npf.export(pv_folder, smooth=True, fmt='vtk', name='NPF', point_scalars=True)

# DRN export, with points
gwf.drn.export(pv_folder, smooth=True, fmt='vtk', name='DRN', point_scalars=True)#doen't work

# 3D Array export
# hk export, with points
gwf.sto.export(pv_folder, smooth=True, fmt='vtk', name='STO', point_scalars=True)

# ghb export, with points
gwf.ghb_0.export(pv_folder, smooth=True, fmt='vtk', name='ghb_0', point_scalars=True)#doen't work

# chd export, with points
gwf.chd.export(pv_folder, smooth=True, fmt='vtk', name='CHD', point_scalars=True)#doen't work

# head export, with points
# vtk.export_heads(gwf, head_file,
#                  pv_folder,
#                  kstpkper=[(0,0), (0, 49), (0, 99), (0, 999)],
#                  point_scalars=True, nanval=(1e30,-1e30))

# vtk.export_heads(gwf, cellbycellfile, pv_folder)