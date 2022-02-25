# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 11:50:13 2021

@author: oscar.sanchez
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
plt.rcParams['font.size'] = 12


case="model_pest"

m_d="master"
t_d = "template"
assert os.path.exists(m_d)," need to run the Pest_inv_unc file first!"
pst = pyemu.Pst(os.path.join(m_d,f"{case}.pst"))



     
obs_filename =case+".0.obs.csv"
obs_df = pd.read_csv(os.path.join(m_d,obs_filename),index_col=0)#read modelled obs
obs = pst.observation_data

def time_series_graphs(prior=True, posterior=False, iter=1, path=None):
    if posterior:
        obs_post_filename =case+f".{iter}.obs.csv"
        obs_post_df = pd.read_csv(os.path.join(m_d,obs_post_filename),index_col=0)#read modelled obs   
    
        
    for obs_group in pst.obs_groups:
        obs_g = obs.loc[obs.obgnme==obs_group,:].copy()#select field/init obs from control file of one obs_group
        obs_g.time=obs_g.time.astype("float64")
        obs_g.sort_values(by="time",inplace=True)#organize values
        fig,ax = plt.subplots(1,1,figsize=(15,2))

        if prior:
            obs_g_df = obs_df.loc[:,obs_g.obsnme]#select from modelled values those that match with obs_group and field_obs(Inlude 0 weight)
            obs_g_df.replace([1e30,3e30,-1e30],np.nan, inplace=True)#mask empty values in GAL or dry cells
            
            [ax.plot((obs_g.time.astype("float64"))/86400,obs_g_df.loc[i,obs_g.obsnme],
                     color='0.6',alpha=0.7,lw=0.1) for i in obs_g_df.index]#prior graphs

        if posterior:
        
            obs_post_g_df = obs_post_df.loc[:,obs_g.obsnme]#select from modelled values those that match with obs_group and field_obs(Inlude 0 weight)
            obs_post_g_df.replace([1e30,3e30,-1e30],np.nan, inplace=True)#mask empty values in GAL or dry cells
            [ax.plot(obs_g.time.astype("float64")/86400,obs_post_g_df.loc[i,obs_g.obsnme],
                 color='blue',alpha=0.5,lw=0.1) for i in obs_post_g_df.index]#posterior graphs

            
        obs_g_f=obs_g.loc[obs_g.weight!=0]#taking out non existen obs
        ax.plot(obs_g_f.time.astype("float64")/86400,obs_g_f.obsval,"x--",
                color='red', ms="1.0",alpha=1,lw=1)#plotting field obs

        ax.set_title("realizations "+obs_group)
        ax.set_xlabel("time[d]")
        ax.set_ylabel(obs_group.split("_")[0])
        # ax.set_xlim(0,21513601)
        plt.setp( ax.xaxis.get_majorticklabels()[::50], rotation=45, horizontalalignment='right' )
        # plt.show()

        if not path==None:
            plt.savefig(path+f"/t_s_{obs_group.split(':')[1]}", dpi=300)
            plt.show()
    
# ['gal_flow_usecol:gal-flow',
#  'gal_flow_usecol:gal_w-flow']
#setting graphs of frecuency flow
def frequency_graphs(prior=True, posterior=False, iter=1, norm=False,path=None):
    obs_gal_nam=obs.loc[obs.obgnme=='flow_gal_usecol:gal-flow'].copy()
    obs_gal_w_nam=obs.loc[obs.obgnme=='flow_gal_w_usecol:gal_w-flow'].copy()
    obs_gal_nam.sort_values(by="time",inplace=True)#organize values
    obs_gal_w_nam.sort_values(by="time",inplace=True)#organize values
    pr=""
    post=""
    legend=[]
    fig=plt.figure()
    num_reals=[]
    if prior:
        pr="_pr"
        obs_gal_df = obs_df.loc[:,obs_gal_nam.obsnme]#select from modelled values those that match with obs_group and field_obs(Inlude 0 weight)
        obs_gal_w_df = obs_df.loc[:,obs_gal_w_nam.obsnme]
        
        obs_gal_df.replace([1e30,3e30,-1e30], 0, inplace=True)#mask empty values
        obs_gal_w_df.replace([1e30,3e30,-1e30], 0, inplace=True)
    
        obs_gal_tot=obs_gal_df.copy()
        for i in obs_gal_df.columns:
            i_j=i.split("gal")
            k=i_j[0]+"gal_w"+i_j[1]+"gal_w"+i_j[2]
            
            obs_gal_tot.loc[:,i]=obs_gal_df.loc[:,i]*1000+obs_gal_w_df.loc[:,k]*1000
    
    
        x_1=obs_gal_tot.min(axis=1)*-1
        plt.hist(x_1, np.arange(int(min(x_1)-1),int(max(x_1)+1),1),density=norm, facecolor='0.5', alpha=1, lw=1, ec="black")
        legend.append("prior")
        num_reals.append(len(x_1))

    if posterior:
        post="_post"
        obs_post_filename =case+f".{iter}.obs.csv"
        obs_post_df = pd.read_csv(os.path.join(m_d,obs_post_filename),index_col=0)#read modelled obs    
        
        obs_post_gal_df = obs_post_df.loc[:,obs_gal_nam.obsnme]#select from modelled values those that match with obs_group and field_obs(Inlude 0 weight)
        obs_post_gal_w_df = obs_post_df.loc[:,obs_gal_w_nam.obsnme]
        
        obs_post_gal_df.replace([1e30,3e30,-1e30], 0, inplace=True)#mask empty values
        obs_post_gal_w_df.replace([1e30,3e30,-1e30], 0, inplace=True)
    
        obs_post_gal_tot=obs_post_gal_df.copy()
        for i in obs_post_gal_df.columns:
            i_j=i.split("gal")
            k=i_j[0]+"gal_w"+i_j[1]+"gal_w"+i_j[2]
            
            obs_post_gal_tot.loc[:,i]=obs_post_gal_df.loc[:,i]*1000+obs_post_gal_w_df.loc[:,k]*1000
    
    
        x_2=obs_post_gal_tot.min(axis=1)*-1
        plt.hist(x_2, np.arange(int(min(x_2)-1),int(max(x_2)+1),1),density=norm, facecolor='b', alpha=0.3, ls='dotted', lw=3, ec="red") 
        legend.append("posterior")
        num_reals.append(len(x_2))

        
        
    plt.xlabel("flow [l/s]")
    plt.ylabel("Frequency ")
    if norm:
        plt.ylabel("Probability density ")
    
    plt.title(f" frequency maximum total gal_flow \n num_reals={num_reals}")
    
    plt.legend(legend)
    # plt.text(30, 6, f"num_reals={num_reals}")
    if not path==None:
            plt.savefig(path+"/hist"+pr+post)    
    plt.plot()




#setting graph for 1to1 of all ensembles
# prior,1to1,phi_pie,phi_progress, and someday  obs_v_sim?
def oneto1_graph(iter=None, ofilename=None, post=False):
    
    
    print_1to1={'0.5': obs_df.replace(-1e30, np.nan)}
    if post:
        obs_filename_post =case+f".{iter}.obs.csv"
        obs_df_post = pd.read_csv(os.path.join(m_d,obs_filename_post),index_col=0)#
        # obs_df_post = obs_df_post.replace(-1e30, np.nan)
        print_1to1["b"]=obs_df_post.replace(-1e30, np.nan)
    
    obsplusnoise_nam =case+".obs+noise.csv"
    
    obsplusnoise_df = pd.read_csv(os.path.join(m_d,obsplusnoise_nam),index_col=0)
    
    if not ofilename==None:
        ofilename=os.path.join(ofilename,"1to1.pdf" )
    pyemu.plot_utils.ensemble_res_1to1(
         print_1to1,
         pst=pst,  #  pyemu.Pst() object -- this defines the groupings of the obs (so you can regroup obs here if you want more/less grnularity on the plots) 
         filename=ofilename,
         base_ensemble=obsplusnoise_df, alpha=0.5
     )
    
    
if __name__ == '__main__':
    it=0
    folder= f"v3_post/i{it}"#os.path.join("../06_jpg/", ) 
    full_path=os.path.join("../06_Jpg/",folder)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    
    oneto1_graph(iter=it, ofilename=full_path, post=True)

    time_series_graphs(prior=True,posterior=True, iter=it, path=full_path)
    
  
    frequency_graphs(prior=True, posterior=True, iter=it, path=full_path)    


