# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:56:44 2016

@author: Lukas
"""

import os
import time

import numpy as np

import or_main
import or_fkt
import or_class
import or_plot_intensity as pltI

# Settings for the calculations
path_skript = os.path.dirname(os.path.realpath(__file__))

settings = or_class.OrSetting()
settings.alpha_inc = np.deg2rad(30)
settings.alpha_start = np.deg2rad(0)
settings.alpha_end = np.deg2rad(90)

settings.beta_inc = np.deg2rad(30)
settings.beta_start = np.deg2rad(0)
settings.beta_end = np.deg2rad(90)

settings.gamma_inc = np.deg2rad(5)
settings.gamma_start = np.deg2rad(0)
settings.gamma_end = np.deg2rad(0)

settings.path_save = os.path.abspath(os.path.join(path_skript,os.pardir)) + "\\Results\\"

# Read Data for both materials form CONTCAR like files
path_folder = path_skript + "\\Cont\\"

MgO_lattice, MgO_cell = or_fkt.read_data(path_folder + "MgO" ,[1,1,1])
MgO_lattice.name = "MgO"

V_lattice, V_cell = or_fkt.read_data(path_folder + "V" ,[1,1,1])
V_lattice.name = "V"

Al2O3_lattice, Al2O3_cell = or_fkt.read_data(path_folder + "Al2O3_cell.txt" ,[1,1,1])
Al2O3_lattice.name = "Al2O3"

Cr2Al_lattice, Cr2Al_cell = or_fkt.read_data(path_folder + "Cr2Al_cell.txt" ,[1,1,1])
Cr2Al_lattice.name = "Cr2Al"

Cr5Al8_lattice, Cr5Al8_cell = or_fkt.read_data(path_folder + "Cr5Al8_cell.txt" ,[1,1,1])
Cr5Al8_lattice.name = "Cr5Al8"


MatB_lattice = V_lattice
MatB_cell = V_cell
MatA_lattice = MgO_lattice
MatA_cell = MgO_cell



hkl = [ 
    [2,2,2]
    ]
start = time.time()
print("Find Orientational Relationship between {} and {}".format(MatA_lattice.name,MatB_lattice.name))
for ent in hkl:
    print("hkl set to {}{}{}".format(ent[0],ent[1],ent[2]))

    path_plt, file_name = or_main.or_gautam_meth(settings,MatA_lattice, MatA_cell, MatB_lattice,MatB_cell,ent)

    max_intens = or_fkt.find_max_inten(path_plt,5)

    or_main.io_gautam_meth(max_intens,MatA_lattice, MatA_cell, MatB_lattice,MatB_cell,ent)

    pltI.plot_intens(path_plt,0,file_name)
    
    print("Total Time: {}".format(time.time()-start))

