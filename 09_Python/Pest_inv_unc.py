# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:03:39 2021

@author: oscar.sanchez
"""
import os
import sys
import flopy as fp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pyemu
import shutil
from shutil import copyfile

def setup_obs():
    
    #procedure to import dipper log head measuments
    obs_path="../04_Xls/Observ"#"../04_Xls/Observ"
    dip_log=pd.read_excel( os.path.join(obs_path,"dip_log.xlsx"),sheet_name=None)
    
    #filter for bad dates due to piezometer leakage
    dip_log["P-Gal-S3"].loc[dip_log["P-Gal-S3"]["P-Gal-S3"]<9.4,"P-Gal-S3"]=None
    
    for i in dip_log.keys():
        dip_log[i]=dip_log[i].groupby('date').mean()
        
    dip_log_d=pd.concat(
        dip_log,
        ignore_index=False).reset_index(level=0).drop(columns="level_0").groupby('date').mean()
    
    #procedure to import vibrating wire head measuments
    vib_wir=pd.read_excel( os.path.join(obs_path,"vib_wire.xlsx"),sheet_name=None)
    for i in vib_wir.keys():
        vib_wir[i]=vib_wir[i].groupby('date').mean()
    
    vib_wir_d=pd.concat(vib_wir,ignore_index=False).groupby('date').mean()
    vib_wir_d_2=vib_wir_d.resample("d").mean()
    
    #filter for bad dates due to piezometer installation and one anomaly
    
    
    vib_wir_d_2.loc[vib_wir_d_2.index<"7/11/2017","P-RN-S7_1"]=None
    vib_wir_d_2.loc[vib_wir_d_2.index<"8/3/2017","P-RN-S7_4"]=None
    vib_wir_d_2.loc[vib_wir_d_2.index<"7/16/2017","P-E3-S6_1"]=None
    vib_wir_d_2.loc[vib_wir_d_2.index<"10/13/2017","P-E3-S6_2"]=None
    vib_wir_d_2.loc[vib_wir_d_2.index<"10/12/2017","P-E3-S6_3"]=None
    vib_wir_d_2.loc[vib_wir_d_2.index<"8/5/2017","P-VR-S_1"]=None
    vib_wir_d_2.loc[(vib_wir_d_2.index >"5/10/2018") & (vib_wir_d_2.index<"7/7/2018"),
                    "P-VR-S_1"]=None#anomaly
    vib_wir_d_2.loc[vib_wir_d_2.index<"6/25/2017","P-VR-S_2"]=None
    
    #procedure to import water level meter head measuments
    wlm=pd.read_excel( os.path.join(obs_path,"WLM.xlsx"),sheet_name=0, index_col="date")
    
    #procedure to import flow measuments
    flow=pd.read_excel( os.path.join(obs_path,"Flow_m_s.xlsx"),sheet_name=0, index_col="date")
    
    
    #merge all dataframes, warning 
    
    obs_heads2=pd.concat([dip_log_d, wlm], ignore_index=False, axis=0, sort=False)#,join="inner")#, verify_integrity=True)
    
    obs_heads3=pd.concat([vib_wir_d_2, obs_heads2], ignore_index=False, axis=1, sort=False)#, verify_integrity=True)
    return obs_heads3

def modif_obs_csv(csv_file):
    tmp_model_ws="temp_pst_model"#erase later
    df_obs=pd.read_csv(os.path.join(tmp_model_ws,csv_file),index_col=0)
    df_obs.loc[:,:]=None
    df_field_mea["time"]=df_field_mea.index-pd.Timestamp('2017-01-01 00:00:00')
    df_field_mea["time"]=df_field_mea["time"].astype("timedelta64[s]")+1 #becausesteady time shift every stress period
    df_field_mea.columns=df_field_mea.columns.str.upper()
    df_obs_final=pd.concat([df_obs,df_field_mea.set_index("TIME")], join="outer", axis=0 )
    
    import geopandas as gpd
    inventory=gpd.read_file("../../05_Vectorial/INV_PAS_V5_DEM.shp")
    inv=inventory[inventory["DEPTH_MEA"]>0]
    inv.reset_index(drop=True, inplace=True)#because we erased some points
    
    df_depth=inv.loc[:,["obs_model","SAMPLE_DEM"]]
    df_depth["obs_model"]=df_depth["obs_model"].str.upper()
    df_depth.set_index("obs_model", inplace=True)
    # df_depth.transpose()
    df_obs_final=-df_obs_final+df_depth.squeeze()
    df_obs_final.index.name="time"
    df_obs_final.fillna(" ", inplace=True)
    df_obs_final.to_csv(os.path.join(tmp_model_ws,csv_file))
    return df_obs_final
    
    
    
    
    
def setup_inv_model(org_ws):
    # print(os.listdir(org_ws))
    exe_name=r"C:\WRDAPP\mf6.2.0\bin\mf6"
    # pyemu.os_utils.run(exe_name, cwd=org_ws)
    tmp_model_ws="temp_pst_model"
    #copy files to avoid corruption
    if os.path.exists(tmp_model_ws):
        shutil.rmtree(tmp_model_ws)
    shutil.copytree(org_ws, tmp_model_ws)
    #identify spatial parameters
    sim = fp.mf6.MFSimulation.load(sim_ws=tmp_model_ws)
    m = sim.get_model("modelo_Norte")
    sr = pyemu.helpers.SpatialReference(delr=m.dis.delr.array,
                                        delc=m.dis.delc.array,
                                        xll=m.dis.xorigin.array,
                                        yll=m.dis.yorigin.array,
                                        epsg="3116"
                                        )
    # print(sr)
    #create instance of pstfrom for pest++
    template_ws = "template"
    df=  modif_obs_csv("modelo_Norte.obs.head.csv")#create file with real obs before cloning folder
    pf = pyemu.utils.PstFrom(original_d=tmp_model_ws,
                             new_d=template_ws,
                             remove_existing=True,
                             longnames=True,
                             spatial_reference=sr,
                             zero_based=False,
                             start_datetime="01-01-2019" #REVIEW
                             )
   
    #add , check real head name, because there may be an error
    
    # df = pd.read_csv(
    #     os.path.join(tmp_model_ws,"modelo_Norte.obs.head.csv"),
    #     index_col=0)
    # df=  modif_obs_csv("modelo_Norte.obs.head.csv")
    
    
    hds_df=pf.add_observations(
        "modelo_Norte.obs.head.csv",
        insfile="heads.csv.ins",
        index_cols="time",
        use_cols=list(df.columns.values),
        prefix="hds")
    
    print(hds_df)
    
    print([f for f in os.listdir(template_ws) if f.endswith(".ins")])

    
    # add parameters
    #set variogram parameters
    pp_cell_space= 2 #each x cells a point is placed
    pp_v = pyemu.geostats.ExpVario(contribution= 1.0,
                                a=max(m.dis.delr.array)*pp_cell_space*3)
    pp_rain_v = pyemu.geostats.GauVario(contribution= 1.0,
                                a=max(m.dis.delr.array)*pp_cell_space*3)
    #set geostructu
    grid_gs = pyemu.geostats.GeoStruct(variograms=pp_v,
                                       transform = "log")
    rech_gs = pyemu.geostats.GeoStruct(variograms=pp_rain_v)
    #what? my recharge doesn't have spatial variation
    fig, ax = plt.subplots(1,1,figsize=(12,8))
    ax.set_title("spatial variogram")
    grid_gs.plot(ax=ax)
    
    #select 
    # maybe separate qbg from qg?
    #this charge the # of layers of every geological layer
    layers=np.loadtxt(os.path.join(tmp_model_ws,"layers"))
    ib=np.zeros(m.dis.idomain.array.shape)
    # assign domain accordingly to the number of layers
    for i in range(layers.size):
        ib[i] = m.dis.idomain.array[int(layers[0:i].sum())]
        ib[i][ib[i]==-1]=0
        print("i= ",i," and idom(j)= ",int(layers[0:i].sum()) )
        
    # list files to modify in calibration/uncertainty
    hk_arr_files = [f for f in os.listdir(tmp_model_ws) if "k_" in f and f.endswith(".txt")]
    vk_arr_files = [f for f in os.listdir(tmp_model_ws) if "kv_" in f and f.endswith(".txt")]
    ss_arr_files = [f for f in os.listdir(tmp_model_ws) if "ss_" in f and f.endswith(".txt")]
    sy_arr_files = [f for f in os.listdir(tmp_model_ws) if "sy_" in f and f.endswith(".txt")]
    # ghb_arr_files = [f for f in os.listdir(tmp_model_ws) if "ghb_" in f and f.endswith(".txt")]#remember why
    rch_arr_files = [f for f in os.listdir(tmp_model_ws) if "rch_" in f and f.endswith(".txt")]
    print(vk_arr_files)
    # assign parameters for instruction file
    for i in range(len(hk_arr_files)):
        pf.add_parameters(filenames=hk_arr_files[i],
                          par_type="pilotpoint",
                          pp_space=pp_cell_space,
                          par_name_base=f"hk_glayer_{i}",
                          pargp=f"hk_glayer_{i}",
                          zone_array=ib[i],
                          upper_bound=10.,
                          lower_bound=0.01,
                          ult_ubound=100,
                          ult_lbound=0.001,
                          spatial_reference=sr,
                          geostruct=grid_gs)
    # for i in range(len(vk_arr_files)):
    #     pf.add_parameters(filenames=hk_arr_files[i],
    #                       par_type="pilotpoint",
    #                       pp_space=pp_cell_space,
    #                       par_name_base=f"hk_glayer_{i}",
    #                       pargp=f"vk_glayer_{i}",
    #                       zone_array=ib[i],
    #                       upper_bound=10.,
    #                       lower_bound=0.1,
    #                       ult_ubound=100,
    #                       ult_lbound=0.01,
    #                       spatial_reference=sr,
    #                       geostruct=grid_gs)
    # for i in range(len(ss_arr_files)):#warning in bounds
    #     pf.add_parameters(filenames=hk_arr_files[i],
    #                       par_type="pilotpoint",
    #                       pp_space=pp_cell_space,
    #                       par_name_base=f"ss_glayer_{i}",
    #                       pargp=f"ss_glayer_{i}",
    #                       zone_array=ib[i],
    #                       upper_bound=10,
    #                       lower_bound=0.1,
    #                       ult_ubound=50,
    #                       ult_lbound=0.01,
    #                       spatial_reference=sr,
    #                       geostruct=grid_gs)
    # for i in range(len(sy_arr_files)):#warning in bounds
    #     pf.add_parameters(filenames=sy_arr_files[i],
    #                       par_type="pilotpoint",
    #                       pp_space=pp_cell_space,
    #                       par_name_base=f"sy_glayer_{i}",
    #                       pargp=f"sy_glayer_{i}",
    #                       zone_array=ib[i],
    #                       upper_bound=3.,
    #                       lower_bound=0.1,
    #                       ult_ubound=5,
    #                       ult_lbound=0.01,
    #                       spatial_reference=sr,
    #                       geostruct=grid_gs)
    # ghb_list=[]
    # for i in range(int(layers.sum())):
    #     ghb_list.append([f for f in os.listdir(tmp_model_ws) if f"ghb_{i}" in f and f.endswith(".txt")])
        
    
    # for i in range(len(ghb_list)):#warning in bounds
    #     pf.add_parameters(filenames=ghb_list[i],
    #                       par_type="constant",
    #                       par_name_base=f"ghb_glayer_{i}",
    #                       pargp=f"ghb _glayer_{i}",
    #                       index_cols=[0,1,2],
    #                       use_cols=[4],
    #                       upper_bound=10.,
    #                       lower_bound=0.1,
    #                       ult_ubound=100,
    #                       ult_lbound=0.01,
    #                       spatial_reference=sr
    #                       )
        
    # #be careful about this setup, review PLEASE
    # pf.add_parameters(filenames=rch_arr_files,
    #                   par_type="constant",
    #                   par_name_base=f"rch_rain",
    #                   pargp=f"rch_rain",
    #                   index_cols=[0,1,2],
    #                   use_cols=[3],
    #                   upper_bound=3,
    #                   lower_bound=0.3,
    #                   ult_ubound=5,
    #                   ult_lbound=0.1,
    #                   spatial_reference=sr,
    #                   )
    
    #add run model command(run once only?), review later
    pf.mod_sys_cmds.append(exe_name)
    pst = pf.build_pst("model_pest.pst")
    
    
    # pst =pf.build_pst()#i suspect, this is necessary only when you modify something, check later
    cov=pf.build_prior()
    x = cov.x.copy()
    x[x==0.0] = np.NaN
    fig = plt.figure(figsize=(12,12))
    im = plt.imshow(x, interpolation="none")
    plt.gca().set_facecolor("k")
    fig.savefig('kcov.png')
    
    # pestpp-ies
    # draw from the prior and save the ensemble in binary format
    # it doesn't work because i have an irregular structured grid
    # pe = pf.draw(300, use_specsim=True)     
    # pe.to_binary(os.path.join(template_ws, "prior.jcb"))
    
    pst = pf.build_pst("model_pest.pst")
    
    # set up control file
    pst.control_data.noptmax=0
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.write(os.path.join(pf.new_d, "model_pest.pst"))
    
    # run with noptmax = 0
    exe_p_name=r"C:\WRDAPP\bin\pestpp-ies"
    pyemu.os_utils.run(exe_p_name + " model_pest.pst", cwd=pf.new_d)
    
    # make sure it ran... what?
    res_file = os.path.join(pf.new_d, "model_pest.base.rei")
    assert os.path.exists(res_file), res_file
    pst.set_res(res_file)
    print(pst.phi)
    
    # now I use noptmax -1 to run prior monte carlo
    pst.control_data.noptmax=-1
    #update files
    pst.write(os.path.join(pf.new_d, "model_pest.pst"))
    
    
    
def run_pest(t_d):
    num_workers=6
    exe_p_name=r"C:\WRDAPP\bin\pestpp-ies"
    pyemu.os_utils.start_workers(t_d,
                                  exe_p_name,#"../10_exe/pestpp-ies.exe",
                                  "model_pest.pst",
                                  num_workers=num_workers,
                                  worker_root=".",
                                  silent_master=False,
                                  verbose=True,
                                  master_dir="master")#silent_master?
    

    
    
    
       
        
    

    
    
    
    
    
if __name__ == "__main__":
    df_field_mea=setup_obs()
    # df=  modif_obs_csv("modelo_Norte.obs.head.csv")
    setup_inv_model("data/modelo_Norte")
    run_pest("template")
    
    