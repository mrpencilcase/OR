# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:56:44 2016

@author: Lukas
"""

import os
import time

import numpy as np

import or_gautam
import plot_angel_Volume as pltI
import or_main
import or_fkt
import or_class
select = 1
# Settings for the calculations
path_skript = os.path.dirname(os.path.realpath(__file__))

settings = or_class.or_setting()
settings.alpha_inc = np.deg2rad(5)
settings.alpha_max = np.deg2rad(90)
settings.beta_inc = np.deg2rad(5)
settings.beta_max = np.deg2rad(90)
settings.gamma_inc = np.deg2rad(5)
settings.gamma_max = np.deg2rad(0)
settings.R_scale = 5
settings.r_scale = 0.25
settings.path_save = os.path.abspath(os.path.join(path_skript,os.pardir)) + "\\Results\\"

# Read Data for both materials form CONTCAR like files
path_folder = path_skript + "\\Cont\\"

MnO_lattice, MnO_unit_cell = or_fkt.read_data(path_folder+"MgO",[1,1,1])
V_lattice, V_unit_cell = or_fkt.read_data(path_folder+"V",[1,1,1])

MnO_lattice.name = "MgO"
V_lattice.name = "V"


if select == 1:
    hkl = [ 
        [4,4,4]
        ]
    start = time.time()

    for ent in hkl:
        print("Start calculation with hkl = {}{}{}".format(ent[0],ent[1],ent[2]))
        path_plt = or_gautam.or_gautam_meth(settings,MnO_lattice, MnO_unit_cell, V_lattice,V_unit_cell,ent)
        pltI.plot_intens(path_plt,0,select)
    
        print("Total Time: {}".format(time.time()-start))

elif select ==2:

    start = time.time()
    R_r = [5]
    r_r = [0.2]
    for R in R_r:
        for r in r_r:
            print("Start calculation with R={:.2f} and r={:.2f}".format(R,r))
            settings.R_scale = R
            settings.r_scale = r
            MnOred = MnO_lattice
            MnOred.a = MnOred.a/2
            MnOred.b = MnOred.b/2
            MnOred.c = MnOred.c/2

            Vred = V_lattice
            Va = Vred.anorm/2
            Vred.a = [-Va,Va,Va]
            Vred.b = [Va,-Va,Va]
            Vred.c = [Va,Va,-Va]

            path_plt = or_main.or_main_fkt(settings,MnOred,Vred)
            pltI.plot_intens(path_plt,0,select)
    
    print("Total Time: {}".format(time.time()-start))
