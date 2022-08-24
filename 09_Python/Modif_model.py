# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:05:11 2022

@author: Oscar
"""

import flopy as fp
import numpy as np
import os
import shutil

def modif_ats_model(run_path="data/modelo_Norte",update=True):
    """
    Function to update MF6 model, when the new flopy version failed, due to some
    incompability between packages. 
    

    Parameters
    ----------
    run_path : TYPE, optional
        Path where the mf6 model will be updated using the ATS package.
        The default is "data/modelo_Norte".
    update : TYPE, optional
        Variable to define if there is a ATS package to update, otherwise a new ATS package will be defined.
        The default is True.

    Returns
    -------
    None.

    """
    

    # run_path="E:/backup"
    run_path="data/modelo_Norte"
    m_d=os.path.join(run_path)#,"master")
    
    workspace=m_d
    # workspace="template"
    
    model_name= "modelo_Norte"
    sim_name="mfsim.nam"
    # # budget_file = model_name + '.cbb'
    # # head_file = model_name + '.hds'
    # budget_file = os.path.join(workspace, budget_file)
    # head_file = os.path.join(workspace, head_file)
    #loading model
    sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
    gwf=sim.get_model(model_name)
    # nsper
    
    tdis=sim.get_package("tdis")
    ims=sim.get_package("ims")
    nsper=tdis.nper.array
    period_ats =[(i, 86400, 1.0e-5, 86400, 2.0, 5.0) for i in range(1, nsper)]   
    
    ims.ats_outer_maximum_fraction=0.05
    if update:
        ats=sim.get_package("ats")
        ats.perioddata=period_ats
        # ats.pname="ats"
    
    else:
        ats = fp.mf6.ModflowUtlats(tdis,maxats=len(period_ats), perioddata=period_ats,pname="ats" )
      
    ats.write()
    tdis.write()
    ims.write()
    
        
def export_head(run_path="data/modelo_Norte"):
    
    import flopy as fp
    """
    Function to get the array of Heads results in a Mf6 model

    Parameters
    ----------
    run_path : TYPE, optional
        path where the model is. The default is ""data/modelo_Norte"".

    Returns 
    -------
    head
    TYPE: 3D numpy array
        Array of heads in every layer for the first stress period.
    """
    
    
    # run_path="E:/backup"
    # run_path="data/modelo_Norte"
    # m_d=os.path.join(run_path,"master")
    m_d=os.path.join(run_path)
    workspace=m_d
    # workspace="template"
    
    model_name= "modelo_Norte"
    sim_name="mfsim.nam"
    budget_file = model_name + '.cbb'
    head_file = model_name + '.hds'
    budget_file = os.path.join(workspace, budget_file)
    head_file = os.path.join(workspace, head_file)
    #loading model
    sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
    gwf=sim.get_model(model_name)
    
    head = gwf.output.head().get_alldata()
    return head[0,:,:,:]+0.1, head.shape[0]
    
    
def write_ext_str_heads(head_0, run_path="E:/Thesis_Runs"):
    """
    Function to write a numpy array in external files for modflow6
    it is expected that run path be the parent folder where the model(s) is(are)

    Parameters
    ----------
    head_0 : TYPE:numpy array
        3D numpy array of heads.
    run_path : TYPE: str, optional
        path where the numpy array will be written . The default is "E:/Thesis_Runs".

    Returns
    -------
    None.

    """
     
    for i in range(head_0.shape[0]):
        np.savetxt(os.path.join(run_path,f"strt_head_{i}.txt"), head_0[i])
        
        
def update_start_head(run_path="E:/Thesis_Runs",folder_update=".", first_update=False):
    import shutil
    import flopy as fp
    """
    Function to use strt external files in run_path to copy in folder_update

    Parameters
    ----------
    run_path : TYPE:str, optional
        path where the parent folder with external strt files are located.
        The default is "E:/Thesis_Runs".
    folder_update : TYPE:str, optional
        name of the folder inside run_path where the model to be updated is located. The default is "\.".
    first_update : TYPE_Boolean, optional
        Variable to define if the model is not set with external files for
        the STRT package, only the first update needs to change strt package.
        The default is False.

    Returns
    -------
    None.

    """
    if folder_update==".":
        workspace="\."
    else:
        workspace=os.path.join(run_path,folder_update)#problem here!!!!!!!!!!
    # workspace="template"
    
    model_name= "modelo_Norte"
    sim_name="mfsim.nam"
    
    strt_files = [f for f in os.listdir(run_path) if "strt" in f and f.endswith(".txt")]
    for i in strt_files:
        shutil.copyfile(os.path.join(run_path,i),
                        os.path.join(workspace, i))
    if first_update:
        
        #loading model
        sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
        gwf=sim.get_model(model_name)
        strt=[f"strt_head_{i}.txt" for i in range(gwf.modelgrid.nlay)] 
        ic=gwf.get_package("ic")
        ic.strt=strt
        ic.write()

def copy_from_out():
    update_start_head(run_path="E:/Thesis_Runs",folder_update=".",first_update=False)

def copy_to_out():
    head_0, numsp=export_head(run_path=".")
    if numsp>1:
        write_ext_str_heads(head_0)
        
def update_ims(orig_wd=os.path.join("E:/Thesis_Runs","worker_0"),run_path="E:/Thesis_Runs", nworkers=10):
    """
    Function to copy quickly an ims file into other model folders,(pest workers) this to test the same solver settings in many models

    Parameters
    ----------
    orig_wd : TYPE:str, optional
        path where the original working directory 
        is located. The default is os.path.join("E:/Thesis_Runs","worker_0").
    run_path : TYPE:str, optional
        paths where the workers are. The default is "E:/Thesis_Runs".
    nworkers : TYPE:int, optional
        Number of workers available to be updated. The default is 10.

    Returns
    -------
    None.

    """
    ims_file = [f for f in os.listdir(orig_wd) if f.endswith(".ims")]
    for i in range(nworkers):
        if os.path.join(orig_wd, ims_file[0]) != os.path.join(run_path, f"worker_{i}", ims_file[0]):
            shutil.copyfile(os.path.join(orig_wd, ims_file[0]),
                                os.path.join(run_path, f"worker_{i}", ims_file[0]))

    m_d=os.path.join(run_path,"worker_0")
    

    print("ready ims")
        
if __name__ == '__main__':
    head_0, numsp =export_head(run_path="data/modelo_Norte")
    # write_ext_str_heads(head_0, run_path="E:/Thesis_Runs")
    # update_start_head(run_path="E:/Thesis_Runs",folder_update="template",first_update=True)
    # update_ims(orig_wd=os.path.join("E:/Thesis_Runs","worker_0"),run_path="E:/Thesis_Runs", nworkers=10)
    
    
    
    
    