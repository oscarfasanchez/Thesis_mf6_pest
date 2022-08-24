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
sys.path.append(os.path.abspath(r"D:\OneDrive - UNIVERSIDAD INDUSTRIAL DE SANTANDER\Maestria\06_Tesis\01_Tesis_Dev\09_Python"))
import Modif_model
print(sys.path)
case="model_pest"

def setup_obs(obs_path="../04_Xls/Observ"):
    """
         it creates a data frame of observations with dates(rows) and
         measument names(columns) using raw files coming from excel. 
         those raw filescomes from dipper logs, vibrating wires and manual
         measurement with water level meters
         
    Parameters
    ----------
    obs_path : string
        path of raw xls field observations.         

    Returns
    -------
    obs_heads3 : Dataframe with time as index and columns as model observations
    (all of them, even model observations without field observations)


    """
    #procedure to import dipper log head measuments
    # obs_path="../04_Xls/Observ"#"../04_Xls/Observ"
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

def modif_obs_csv(csv_file, df_field_mea , lower=False):
    """
    

    Parameters
    ----------
    csv_file : string
        csv name of model(modflow6) output observations.
    df_field_mea : Dataframe
        dataframe of field measurament with obs name in colums and dates in index.
    lower : Boolean, optional
        it decides if column are in uppercase or lowercase. The default is False.

    Returns
    -------
    df_obs_final : TYPE
        dataframe of field observations in mf6 style.

    """
    tmp_model_ws="temp_pst_model"#erase later
    df_obs=pd.read_csv(os.path.join(tmp_model_ws,csv_file),index_col=0)
    df_obs.loc[:,:]=None
    df_field_mea["time"]=df_field_mea.index-pd.Timestamp('2017-01-01 00:00:00')
    df_field_mea["time"]=df_field_mea["time"].astype("timedelta64[s]")+1 #becausesteady time shift every stress period
    df_field_mea.columns=df_field_mea.columns.str.upper()
    df_obs_final=pd.concat([df_obs,df_field_mea.set_index("TIME")], join="outer", axis=0 )
    df_obs_final=df_obs_final.reset_index().groupby("index").max()
    
    import geopandas as gpd
    inventory=gpd.read_file("../../05_Vectorial/INV_PAS_V5_DEM.shp")
    inv=inventory[inventory["DEPTH_MEA"]>0]
    inv.reset_index(drop=True, inplace=True)#because we erased some points
    
    df_depth=inv.loc[:,["obs_model","SAMPLE_DEM"]]
    df_depth["obs_model"]=df_depth["obs_model"].str.upper()#to operate later
    df_depth.set_index("obs_model", inplace=True)
    
    df_obs_final=-df_obs_final+df_depth.squeeze()
    df_obs_final.index.name="time"
    df_obs_final.fillna(0, inplace=True)
    if lower==True:
        df_obs_final.columns = df_obs_final.columns.str.lower()
    # df_obs_final.to_csv(os.path.join(tmp_model_ws,csv_file))#used for modify original csv output
    
    return df_obs_final
    
    
    
    
    
def setup_inv_model(org_ws, template_ws, df_field_meas_2, updt_obs_field=True, run_type="ies", restart_strt=False):
    """
    

    Parameters
    ----------
    org_ws : TYPE
        DESCRIPTION.
    template_ws : TYPE
        DESCRIPTION.
    df_field_meas_2 : TYPE
        DESCRIPTION.
    updt_obs_field : TYPE, optional
        DESCRIPTION. The default is True.
    run_type : TYPE, optional
        DESCRIPTION. The default is "ies".
    restart_strt : TYPE:boolean, optional
        flag variable to update strt in every new model. The default is "False".

    Returns
    -------
    None.

    """
    print("path_setup_func ",sys.path)
    # print(os.listdir(org_ws))
    exe_name=r"C:\WRDAPP\mf6.3.0\bin\mf6"
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
    template_ws =  os.path.join(template_ws,"template")
    # df=  modif_obs_csv("modelo_Norte.obs.head.csv")#create file with real obs before cloning folder
    pf = pyemu.utils.PstFrom(original_d=tmp_model_ws,
                             new_d=template_ws,
                             remove_existing=True,
                             longnames=True,
                             spatial_reference=sr,
                             zero_based=False,
                             start_datetime="01-01-2019" #REVIEW
                             )
   
    #add , check real head name, because there may be an error
    
    df = pd.read_csv(
        os.path.join(tmp_model_ws,"modelo_Norte.obs.head.csv"),
        index_col=0)
    # df=  modif_obs_csv("modelo_Norte.obs.head.csv")
    
    # add observations
    
    hds_df=pf.add_observations(
        "modelo_Norte.obs.head.csv",
        insfile="heads.csv.ins",
        index_cols="time",
        use_cols=list(df.columns.values),
        prefix="hds")
    # hds_df.loc[hds_df["obsval"]==0,["weight"]]=0 # to disable un existent observations
    # pf.obs_dfs[0]["weight"]=hds_df
    print(hds_df)
    

    df_gal = pd.read_csv(
        os.path.join(tmp_model_ws,"mod_drn_gal_obs.csv"),
        index_col=0)
    gal_df=pf.add_observations(
        "mod_drn_gal_obs.csv",
        insfile="gal_flow.csv.ins",
        index_cols="time",
        use_cols=list(df_gal.columns.values),
        prefix="flow_gal")   
    
    df_w_gal = pd.read_csv(
        os.path.join(tmp_model_ws,"mod_drn_gal_w_obs.csv"),
        index_col=0)
    gal_w_df=pf.add_observations(
        "mod_drn_gal_w_obs.csv",
        insfile="gal_w_flow.csv.ins",
        index_cols="time",
        use_cols=list(df_w_gal.columns.values),
        prefix="flow_gal_w")

    print([f for f in os.listdir(template_ws) if f.endswith(".ins")])

    # add parameters
    #set variogram parameters
    pp_cell_space= 8 #each x cells a point is placed
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
    ib=np.zeros(tuple([layers.shape[0],int(m.dis.nrow.data), m.dis.ncol.data]))#m.dis.idomain.array.shape)
    # assign domain accordingly to the number of layers
    for i in range(layers.size):
        ib[i] = m.dis.idomain.array[int(layers[0:i].sum())]
        ib[i][ib[i]==-1]=0
        print("i= ",i," and idom(j)= ",int(layers[0:i].sum()) )
    
   

        
    
    # list files to modify in calibration/uncertainty
    hk_arr_files = [f for f in os.listdir(tmp_model_ws) if "k_" in f and f.endswith(".txt")]
    hk_arr_files[0], hk_arr_files[2] = hk_arr_files[2], hk_arr_files[0]
    
    vk_arr_files = [f for f in os.listdir(tmp_model_ws) if "kv_" in f and f.endswith(".txt")]
    vk_arr_files[0], vk_arr_files[2] = vk_arr_files[2], vk_arr_files[0]
    
    ss_arr_files = [f for f in os.listdir(tmp_model_ws) if "ss_" in f and f.endswith(".txt")]
    ss_arr_files[0], ss_arr_files[2] = ss_arr_files[2], ss_arr_files[0]
    
    sy_arr_files = [f for f in os.listdir(tmp_model_ws) if "sy_" in f and f.endswith(".txt")]
    sy_arr_files[0], sy_arr_files[2] = sy_arr_files[2], sy_arr_files[0]
    
    # ghb_arr_files = [f for f in os.listdir(tmp_model_ws) if "ghb_" in f and f.endswith(".txt")]#remember why
    rch_arr_files = [f for f in os.listdir(tmp_model_ws) if "rch_" in f and f.endswith(".txt")]
    print(vk_arr_files)
    
    # list files to define pilot points 
    #if we GLM we need less parameters, so we use custom PP and parameterization
    if run_type== "glm_fosm":
        pp_model_ws="../04_Xls/PP"
        pp_arr_files = [f for f in os.listdir(pp_model_ws) if "pp_" in f and f.endswith(".csv")]
        pp_cell_space = pp_arr_files
        for i in range(len(pp_arr_files)):
            shutil.copyfile(os.path.join(pp_model_ws,pp_cell_space[i]),
                            os.path.join(template_ws, pp_cell_space[i]))
        # define type of parameterization
        ss_par_type=["constant"]*layers.size
        num_sy_pp=2
        sy_par_type=["pilotpoint"]*num_sy_pp+["constant"]*(layers.size-num_sy_pp)
        
    else:
        pp_cell_space=[pp_cell_space]*layers.size    
        # define type of parameterization
        ss_par_type=["pilotpoint"]*layers.size
        sy_par_type=["pilotpoint"]*layers.size
        
    # assign parameters for instruction file
    for i in range(len(hk_arr_files)):
        pf.add_parameters(filenames=hk_arr_files[i],
                          par_type="pilotpoint",
                          pp_space=pp_cell_space[i],
                          par_name_base=f"hk_glayer_{i}",
                          pargp=f"hk_glayer_{i}",
                          zone_array=ib[i],
                          upper_bound=500.,
                          lower_bound=0.005,
                          ult_ubound=1e-4,
                          ult_lbound=1e-10,
                          spatial_reference=sr,
                          geostruct=grid_gs)
    for i in range(len(vk_arr_files)):
        pf.add_parameters(filenames=vk_arr_files[i],
                          par_type="pilotpoint",
                          pp_space=pp_cell_space[i],
                          par_name_base=f"vk_glayer_{i}",
                          pargp=f"vk_glayer_{i}",
                          zone_array=ib[i],
                          upper_bound=10.,
                          lower_bound=0.1,
                          ult_ubound=100,
                          ult_lbound=0.01,
                          spatial_reference=sr,
                          geostruct=grid_gs)
    for i in range(len(ss_arr_files)):#warning in bounds
        pf.add_parameters(filenames=ss_arr_files[i],
                          par_type=ss_par_type[i],
                          pp_space=pp_cell_space[i],
                          par_name_base=f"ss_glayer_{i}",
                          pargp=f"ss_glayer_{i}",
                          zone_array=ib[i],
                          upper_bound=10,
                          lower_bound=0.1,
                          ult_ubound=1e-2,
                          ult_lbound=1e-6,
                          spatial_reference=sr,
                          geostruct=grid_gs)
    for i in range(len(sy_arr_files)):#warning in bounds
        pf.add_parameters(filenames=sy_arr_files[i],
                          par_type=sy_par_type[i],
                          pp_space=pp_cell_space[i],
                          par_name_base=f"sy_glayer_{i}",
                          pargp=f"sy_glayer_{i}",
                          zone_array=ib[i],
                          upper_bound=3,
                          lower_bound=0.3,
                          ult_ubound=0.5,
                          ult_lbound=0.01,
                          spatial_reference=sr,
                          geostruct=grid_gs)
    ghb_list=[]
    for i in range(int(layers.sum())):
        ghb_list.append([f for f in os.listdir(tmp_model_ws) if f"ghb_{i}" in f and f.endswith(".txt")])
        
    
    for i in range(len(ghb_list)):#warning in bounds
        pf.add_parameters(filenames=ghb_list[i],
                          par_type="constant",
                          par_name_base=f"ghb_glayer_{i}",
                          pargp=f"ghb_glayer_{i}",
                          index_cols=[0,1,2],
                          use_cols=[4],
                          upper_bound=100.,
                          lower_bound=0.01,
                          ult_ubound=100,
                          ult_lbound=0,
                          spatial_reference=sr
                          )
    ubnd=[1.157e-6,3.48e-8]    # transient 100% daily percol, steady 1100 mm/yr
    #be careful about this setup, review PLEASE
    for i in range(len(rch_arr_files)):        
        pf.add_parameters(filenames=rch_arr_files[i],
                          par_type="constant",
                          par_name_base=f"rch_rain",
                          pargp=f"rch_rain",
                          index_cols=[0,1,2],
                          use_cols=[3],
                          upper_bound=4,
                          lower_bound=0.25,
                          ult_ubound=ubnd[i],
                          ult_lbound=0,
                          spatial_reference=sr,
                          )
    
    #add run model command(run once only?), review later
    pf.mod_sys_cmds.append(exe_name)
    
    if restart_strt:
        pf.add_py_function("Modif_model.py",
                           'update_start_head(run_path="E:/Thesis_Runs",folder_update=".",first_update=False)', is_pre_cmd=None)
        pf.add_py_function("Modif_model.py",
                           'export_head(run_path=".")', is_pre_cmd=None)
        pf.add_py_function("Modif_model.py",
                           'write_ext_str_heads(head_0)', is_pre_cmd=None)   
        pf.add_py_function("Modif_model.py",
                           'copy_from_out()', is_pre_cmd=True)
        pf.add_py_function("Modif_model.py",
                           'copy_to_out()', is_pre_cmd=False)
    
    
    pst = pf.build_pst(f"{case}.pst")
    
    
    # pst =pf.build_pst()#i suspect, this is necessary only when you modify something, check later
    pf.build_prior()
    # cov=pf.build_prior()
    # x = cov.x.copy()
    # x[x==0.0] = np.NaN
    # fig = plt.figure(figsize=(12,12))
    # im = plt.imshow(x, interpolation="none")
    # plt.gca().set_facecolor("k")
    # fig.savefig('kcov.png')
    
    # pestpp-ies
    # draw from the prior and save the ensemble in binary format
    # it doesn't work because i have an irregular structured grid
    # pe = pf.draw(300, use_specsim=True)     
    # pe.to_binary(os.path.join(template_ws, "prior.jcb"))
    
    if updt_obs_field:
        #process to update observations with field values
        
        # we use a functionto set field obs
        df_field_mf6=  modif_obs_csv("modelo_Norte.obs.head.csv",
                                     df_field_meas_2 ,lower=True)
        df_field_mf6=df_field_mf6.reset_index().melt(id_vars="time", var_name="usecol2")
        

        df_weight = pd.read_excel("../04_Xls/obs_error.xlsx")
        df_weight.usecol2 = df_weight.usecol2.str.lower()
        df_weight["weight2"]=1/(df_weight.desvstd**2)       
        df_field_mf6=df_field_mf6.merge(df_weight, on="usecol2", how="left")
        
        
        df_field_mf6=df_field_mf6.sort_values(["usecol2", "time"])
        df_field_mf6=df_field_mf6.astype("str")
        df_field_mf6=df_field_mf6.astype("object")
        
        df_pst_obs=pst.observation_data.copy()# i need to modify this
        df_pst_obs["usecol2"]=None
        for i in range(df_pst_obs.shape[0]):
            df_pst_obs["usecol2"][i]= df_pst_obs.obgnme[i].split(":")[3]
            
        df_field_mf6["weight2"]=df_field_mf6.weight2.astype("float64")
        df_pst_obs2=df_pst_obs.merge(df_field_mf6, how= "outer", on=["usecol2", "time"])
    

        
        
        df_pst_obs2.loc[df_pst_obs2.value=="0.0","weight"]=0
        df_pst_obs2.value = df_pst_obs2.value.fillna("0.0")
        df_pst_obs2.weight2 = df_pst_obs2.weight2.fillna(0.0)
        
        df_pst_obs2.loc[df_pst_obs2.usecol2=="gal-flow","weight"]=0.0
        df_pst_obs2.loc[df_pst_obs2.usecol2=="gal_w-flow","weight"]=0.0
        df_pst_obs2.loc[df_pst_obs2.obsval==3e+30,"obsval"]=0.0
        
        df_pst_obs2.loc[df_pst_obs2.value!="0.0","weight"]=df_pst_obs2["weight2"]
        df_pst_obs2.loc[df_pst_obs2.value!="0.0","obsval"]=df_pst_obs2["value"]

        df_pst_obs2=df_pst_obs2.set_index("obsnme")
        df_pst_obs2.obsval=df_pst_obs2.obsval.astype("float64")#to avoid problems with pst.set_res
        
        
        
        #modifying actual pest values finallY!!!
        #Arreglar
        pst.observation_data.loc[pst.observation_data.index==df_pst_obs2.index,"obsval"]=df_pst_obs2["obsval"]#update values
        pst.observation_data["weight"]=df_pst_obs2["weight"]#update weights
    
    # pst = pf.build_pst(f"{case}.pst")
    
    # set up control file
    
    pst.control_data.noptmax=0
    if run_type== "ies":
        # pst.svd_data.maxsing = 150
        pst.pestpp_options["ies_autoadaloc"]=True
        pst.pestpp_options["ies_num_threads"]=6
        pst.pestpp_options["additional_ins_delimiters"] = ","
        # pst.pestpp_options["ies_bad_phi"]=1e25
        pst.pestpp_options["ies_num_reals"]=250
        pst.pestpp_options["parcov"] = "{}.prior.cov".format(case) 
        pst.write(os.path.join(pf.new_d, f"{case}.pst"))
        
        # run with noptmax = 0 '''??
        exe_p_name=r"C:\WRDAPP\bin\pestpp-ies"
        # pyemu.os_utils.run(exe_p_name + f" {case}.pst", cwd=pf.new_d)
    
        # make sure it ran... what?
        # res_file = os.path.join(pf.new_d, f"{case}.base.rei")
        # assert os.path.exists(res_file), res_file
        # pst.set_res(res_file)
        # print(pst.phi)
        
    elif run_type=="glm_fosm":
        
        
        pst.control_data.pestmode = "regularization"
        # pst.pestpp_options["n_iter_base"] = 1
        # pst.pestpp_options["n_iter_super"] = 2 
        pst.pestpp_options["glm_num_reals"] = 120 
        pst.pestpp_options["glm_iter_mc"] = "true"
        pst.pestpp_options["glm_accept_mc_phi"] = "true"
        pst.pestpp_options["max_run_fail"] = 1
        
        # pst.pestpp_options["base_jacobian"] = os.path.join("jcb_glm","{}.jcb".format(case))#I'll put it manually when i need it
        pst.pestpp_options["parcov"] = "{}.prior.cov".format(case) 
        pyemu.helpers.zero_order_tikhonov(pst)
        pst.write(os.path.join(pf.new_d, f"{case}.pst"))
        exe_p_name=r"C:\WRDAPP\bin\pestpp-glm"
        
    
    #debugging solver issues
   
    pst.pestpp_options["panther_agent_freeze_on_fail"] = "False"
    # now I use noptmax -1 to run prior monte carlo
    #noptmax 0 JUST run once
    pst.control_data.noptmax=4
    #update files
    pst.write(os.path.join(pf.new_d, f"{case}.pst"))
    
    
    
def run_pest(t_d, run_type="ies"):
    num_workers=11
    if run_type== "ies":
        exe_p_name=r"C:\WRDAPP\bin\pestpp-ies"
    elif run_type=="glm_fosm":
        exe_p_name=r"C:\WRDAPP\bin\pestpp-glm"
        
    pyemu.os_utils.start_workers(os.path.join(t_d,"template"),
                                  exe_p_name,#"../10_exe/pestpp-ies.exe",
                                  os.path.join("{}.pst".format(case)),
                                  num_workers=num_workers,
                                  worker_root=t_d,
                                  silent_master=False,
                                  verbose=True,
                                  master_dir=os.path.join(t_d,"master"))#silent_master?
    
def pest_graphs(m_d):
    pst_a = pyemu.Pst(os.path.join(m_d,"{}.pst".format(case)))
    pst_a.plot(kind='1to1',)
    plt.savefig("1to1")
    plt.savefig("obs_v_sim")
    pst_a.plot(kind="prior")
    pst_a.plot(kind="phi_pie")

    # pst_a.plot()
    pst_a.get_res_stats()
    pst_a.phi
    pst_a.phi_components
    pst_a.phi_components_normalized
    df_residuals=pst_a.res
    obs_summary=pst_a.write_obs_summary_table()
    par_summary=pst_a.write_par_summary_table()
    
    
def update_strt_first_time():
    head_0, numsp =Modif_model.export_head(run_path="data/modelo_Norte")
    Modif_model.write_ext_str_heads(head_0, run_path="E:/Thesis_Runs")
    Modif_model.update_start_head(run_path="E:/Thesis_Runs",folder_update="template",first_update=True)
       
        
    

    
    
    
    
    
if __name__ == "__main__":
    
    
    run_path="E:/Thesis_Runs"
    run_type="glm_fosm"#"ies" , "glm_fosm"
    # df_field_meas=setup_obs("../04_Xls/Observ")
    # setup_inv_model("data/modelo_Norte",run_path, df_field_meas, updt_obs_field=True, run_type=run_type, restart_strt=True)
    update_strt_first_time()
    run_pest(run_path, run_type=run_type)
    # pest_graphs(os.path.join(run_path,"master"))
    
    