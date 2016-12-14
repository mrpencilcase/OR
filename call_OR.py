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
import os
import plot_angel_Volume as pltI

select = 1

if select == 1:
    path_skript = os.path.dirname(os.path.realpath(__file__))
    path_folder = path_skript + "\\Cont\\"
    settings = or_class.or_setting()
    settings.alpha_inc = np.deg2rad(5)
    settings.alpha_max = np.deg2rad(90)
    settings.beta_inc = np.deg2rad(5)
    settings.beta_max = np.deg2rad(90)
    settings.gamma_inc = np.deg2rad(5)
    settings.gamma_max = np.deg2rad(0)
    settings.path_save = os.path.abspath(os.path.join(path_skript,os.pardir)) + "\\Results\\"

    V_lattice, V_unit_cell = or_fkt.read_data(path_folder+"V",[1,1,1])
    MnO_lattice, MnO_unit_cell = or_fkt.read_data(path_folder+"MgO",[1,1,1])
    MnO_lattice.name = "MgO"
    V_lattice.name = "V"

    hkl = [ [2,2,2],
            [3,3,3],
                    ]


    start = time.time()
    for ent in hkl:
        or_gautam.or_gautam_meth(settings,MnO_lattice, MnO_unit_cell, V_lattice,V_unit_cell,ent)

    print("Total Time: {}".format(time.time()-start))

elif select == 2:
    path_skript = os.path.dirname(os.path.realpath(__file__))
    path_save = os.path.abspath(os.path.join(path_skript,os.pardir)) + "\\Results\\"
    dat_name = "MgO_V_HKL_333_d0_20.0"
    pltI.plot_intens(path_save + dat_name + ".dat",0)
