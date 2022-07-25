# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:05:11 2022

@author: Oscar
"""

import flopy as fp
import numpy as np
import os
import shutil

def modif_ats_model(update=True):
    

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
    
        
def export_head():
        run_path="E:/backup"
        # run_path="data/modelo_Norte"
        m_d=os.path.join(run_path,"master")
        
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
        return head[0,:,:,:]
    
    
def write_ext_str_heads(head_0, run_path="E:/Thesis_Runs"):
    
               

    for i in range(head_0.shape[0]):
        np.savetxt(os.path.join(run_path,f"strt_head_{i}.txt"), head_0[i])
        
def update_start_head(run_path="E:/Thesis_Runs"):
        # run_path="E:/backup"
        # run_path="data/modelo_Norte"
        m_d=os.path.join(run_path,"worker_0")
        
        workspace=m_d
        # workspace="template"
        
        model_name= "modelo_Norte"
        sim_name="mfsim.nam"
        
        strt_files = [f for f in os.listdir(run_path) if "strt" in f and f.endswith(".txt")]
        for i in strt_files:
            shutil.copyfile(os.path.join(run_path,i),
                            os.path.join(workspace, i))
 
        #loading model
        sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
        gwf=sim.get_model(model_name)
        strt=[f"strt_head_{i}.txt" for i in range(6)] 
        ic=gwf.get_package("ic")
        ic.strt=strt
        ic.write()
        
def update_drn(run_path="E:/Thesis_Runs", mult=10):
        m_d=os.path.join(run_path,"worker_0")
        
        workspace=m_d
        # workspace="template"
        
        model_name= "modelo_Norte"
        sim_name="mfsim.nam"
        
        sim = fp.mf6.MFSimulation.load(sim_name=sim_name, sim_ws=workspace, exe_name= r"C:\WRDAPP\mf6.2.0\bin\mf6")
        gwf=sim.get_model(model_name)  
        
        drn=gwf.get_package("drn")
        drn.set_all_data_external
        drn.write()
        print("ready")
        
if __name__ == '__main__':
    # head_0=export_head()
    # write_ext_str_heads(head_0)
    # update_start_head()
    update_drn
    
    
    
    
    