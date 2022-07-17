# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 13:51:56 2022

@author: Oscar
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
plt.rcParams['font.size'] = 12


case="model_pest"
run_path="E:/backup"
m_d=os.path.join(run_path,"master")
t_d = os.path.join(run_path,"template")
assert os.path.exists(m_d)," need to run the Pest_inv_unc file first!"
pst = pyemu.Pst(os.path.join(m_d,f"{case}.pst"))

def glm_ensemble_plots(m_d, it):
    pst_a = pyemu.Pst(os.path.join(m_d,"{}.pst".format(case)))
    
    # #NSMC
    # df = df=pd.read_csv(os.path.join(m_d,"{}.post.obsen.csv".format(case)),index_col=0)
    # oe = pyemu.ObservationEnsemble.from_dataframe(pst=pst,df=df)
    # ax = oe.phi_vector.hist()
    
    #FOSM
    dfp = df=pd.read_csv(os.path.join(m_d,"{}.{}.par.usum.csv".format(case, it)),index_col=0)
    print(10**dfp["post_mean"])
    plt.hist(10**dfp["post_mean"],bins=30)
    plt.show()
    # pst_a.plot()
    
    
if __name__ == '__main__':
    it=0
    folder= f"Glm_v3/i{it}"#os.path.join("../06_jpg/", ) 
    full_path=os.path.join("../06_Jpg/",folder)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
        
    glm_ensemble_plots(m_d,3)