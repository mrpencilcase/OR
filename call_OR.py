# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:56:44 2016

@author: Lukas
"""


import numpy as np
import or_main
import or_fkt
import time
import os          
pathAl2O3a = "C:/Users/Lukas/Documents/Diplomarbeit/Results/BaseStructures/alpha_Al2O3"
pathCr5Al8 = "C:\Users\Lukas\Documents\Diplomarbeit\Results\BaseStructures\Cr5Al8_non_mag"       
cellAl2O3a =[1,1,1]
cellCr5Al8 =[1,1,1]
Al2O3a = or_fkt.read_data(pathAl2O3a,cellAl2O3a)
Cr5Al8 = or_fkt.read_data(pathCr5Al8,cellCr5Al8)
Al2O3a.name="alphaAl2O3"
Cr5Al8.name = "Cr5Al8"

testMatA = or_fkt.def_lattice("CoO",4.252,4.252,4.252,90,90,90)
testMatB = or_fkt.def_lattice("YSZ",5.127,5.127,5.127,90,90,90)





# settings for OR calculations: Angles
orset  = or_fkt.or_setting() 
orset.alpha_max = np.deg2rad(90)
orset.alpha_inc = np.deg2rad(5)
orset.beta_max = np.deg2rad(90)
orset.beta_inc = np.deg2rad(5)
orset.gamma_max = np.deg2rad(90)
orset.gamma_inc = np.deg2rad(5)



R =[2,3,5,10]
r =[0.1,0.2,0.3]
path_base="C:\Users\Lukas\Documents\Diplomarbeit\Results\OrientationalRelations\\"
orset.path_save="{}_{}_{}".format(path_base,testMatA.name,testMatB.name)
if not os.path.exists(orset.path_save):
    os.makedirs(orset.path_save)

for x in R:
    for y in r:
        print("R       : {}".format(x))
        print("r       : {}".format(y))
        orset.r_scale = y
        orset.R_scale = x
        start = time.time()        
        or_main.or_main_fkt(orset,testMatA,testMatB)
        print ("Runtime : {:.2f} sec".format(time.time()-start))
        print ("Time    : {}".format(time.strftime("%a, %d %b %Y %H:%M:%S ", time.gmtime())))        
        print ("")