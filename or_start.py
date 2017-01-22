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
settings.alpha_inc = np.deg2rad(5)
settings.alpha_start = np.deg2rad(0)
settings.alpha_end = np.deg2rad(90)

settings.beta_inc = np.deg2rad(5)
settings.beta_start = np.deg2rad(0)
settings.beta_end = np.deg2rad(90)

settings.gamma_inc = np.deg2rad(5)
settings.gamma_start = np.deg2rad(0)
settings.gamma_end = np.deg2rad(90)

settings.path_save = os.path.abspath(os.path.join(path_skript,os.pardir)) + "/Results"

settings.numb_intens = 5 

# Read Data for both materials form CONTCAR like files
path_folder = path_skript + "/Cont/"

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

Al_lattice, Al_cell = or_fkt.read_data(path_folder + "Al" ,[1,1,1])
Al_lattice.name = "Al"

Ag_lattice, Ag_cell = or_fkt.read_data(path_folder + "Ag" ,[1,1,1])
Ag_lattice.name = "Ag"


#MatB_lattice = Cr2Al_lattice
#MatB_cell = Cr2Al_cell
#MatA_lattice = Al2O3_lattice
#MatA_cell = Al2O3_cell


#MatB_lattice = Cr5Al8_lattice
#MatB_cell = Cr5Al8_cell
#MatA_lattice = Al2O3_lattice
#MatA_cell = Al2O3_cell

MatB_lattice = V_lattice
MatB_cell = V_cell
MatA_lattice = MgO_lattice
MatA_cell = MgO_cell

#MatB_lattice = Ag_lattice
#MatB_cell = Ag_cell
#MatA_lattice = Al_lattice
#MatA_cell = Al_cell



hkl = [
    [2,2,2]
    ]

start = time.time()

print("Find Orientational Relationship between {} and {}".format(MatA_lattice.name,MatB_lattice.name))

for ent in hkl:
    print("hkl set to {}{}{}".format(ent[0],ent[1],ent[2]))

    # determin the overlapping intensities 
    path_plt, settings.name = or_main.or_gautam_meth(settings,MatA_lattice, MatA_cell, MatB_lattice,MatB_cell,ent)
    pltI.plot_intens(path_plt,0,settings.name)

    # closer look into the orientaions with the largest intensities
    max_intens = or_fkt.find_max_inten(path_plt,5)
    io_data = or_main.io_gautam_meth(settings,max_intens,MatA_lattice, MatA_cell, MatB_lattice,MatB_cell,ent)
        
    
    print("Total Time: {}".format(time.time()-start))

