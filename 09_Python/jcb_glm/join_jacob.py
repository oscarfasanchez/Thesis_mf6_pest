import pyemu
import os
path1=os.path.join("01_orig" ,"model_pest.jcb" )
jco = pyemu.Jco.from_binary(path1)


path2=os.path.join("02_complement" ,"model_pest.jcb" )
jco2 = pyemu.Jco.from_binary(path2)

mats=[jco, jco2]
jco_new=pyemu.mat.mat_handler.concat(mats)
jco_new.to_coo("model_pest2.jcb")
