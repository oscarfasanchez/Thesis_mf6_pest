# -*- coding: utf-8 -*-
"""
Created on Wed May 11 18:56:10 2022

@author: Oscar
"""

import pyemu
import os
import subprocess as sp
import time
import shutil

case="model_pest"
exe_p_name=r"C:\WRDAPP\bin\pestpp-glm"
cwd="master"#os.path.join("Thesis_Runs","template")
os.chdir(cwd)
with open("output.txt","w") as f:
    sp.Popen([exe_p_name, case, "/r", "/h", ":4004"], stdout=f, text=True)

time.sleep(1.5)
num_workers=11
worker_root=".."
for i in range(num_workers):
    new_worker_dir = os.path.join(worker_root, "worker_{0}".format(i))
    if os.path.exists(new_worker_dir):
        shutil.rmtree(new_worker_dir)
    shutil.copytree(os.path.join(worker_root,"template"), new_worker_dir)    
    os.chdir(new_worker_dir)
    sp.Popen([exe_p_name, case, "/h", "localhost:4004"])
    