# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:56:44 2016

@author: Lukas
"""


import numpy as np
import or_main
import or_fkt
import or_class
import time
import os     
import scattering_factors_tool as sftool
import or_gautam

settings = or_class.or_setting()
settings.alpha_inc = 5
settings.alpha_max = 90
settings.beta_inc = 5
settings.beta_max = 90
settings.gamma_inc = 5
settings.gamma_max = 90

path_folder = "C:\Users\Lukas\Documents\Diplomarbeit\\findOR\CONTCAR\\"
fuckpath = "C:\Users\Lukas\Documents\Diplomarbeit\fuckyoundOR\CONTCAR/"
print(fuckpath)
V_lattice, V_unit_cell = or_fkt.read_data(path_folder+"V",[1,1,1])
MnO_lattice, MnO_unit_cell = or_fkt.read_data(path_folder+"MgO",[1,1,1])


start = time.time()

or_gautam.or_gautam_meth(settings,MnO_lattice, MnO_unit_cell, V_lattice,V_unit_cell)

print("Total Time: {}".format(time.time()-start))
