# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:23:51 2016

@author: Lukas
"""

import or_fkt 
import numpy as np

sfa = 4.252
aA = [sfa/2,     0,     0]
bA = [    0, sfa/2,     0]
cA = [    0,     0, sfa/2]

sfB = 5.127
aB = [sfB/4, sfB/4, sfB/4]
bB = [sfB/4, sfB/4,-sfB/4]
cB = [sfB/4,-sfB/4, sfB/4]

testMatA = or_fkt.def_lattice("CoO",aA,bA,cA)
testMatB = or_fkt.def_lattice("YSZ",aB,bB,cB)

latticeA = np.vstack((testMatA.a,testMatA.b,testMatA.c))
latticeB = np.vstack((testMatB.a,testMatB.b,testMatB.c))

latticeA_rec = or_fkt.reziprocal_lattice(latticeA)
hkl = [1,1,1]
g = hkl[0]*latticeA_rec[0]+hkl[1]*latticeA_rec[1]+hkl[2]*latticeA_rec[2]

AlO = [ ["Al" , 0.342,0.2553,1.3240] ,
        ["O"  , 0.342,0.2553,1.3240] ]
        
        
pathAl2O3a = "C:/Users/Lukas/Documents/Diplomarbeit/Results/BaseStructures/alpha_Al2O3"
pathCr5Al8 = "C:\Users\Lukas\Documents\Diplomarbeit\Results\BaseStructures\Cr5Al8_non_mag"       
cellAl2O3a =[1,1,1]
cellCr5Al8 =[1,1,1]
Al2O3_lattice, Al2O3_unitcell = or_fkt.read_data(pathAl2O3a,cellAl2O3a)
Cr5Al8_lattice, Cr5Al8_unitcell = or_fkt.read_data(pathCr5Al8,cellCr5Al8)
Al2O3_lattice.name="alphaAl2O3"
Cr5Al8_lattice.name = "Cr5Al8"

Al2O3_lattice_rec = or_fkt.reziprocal_lattice(np.vstack((Al2O3_lattice.a,Al2O3_lattice.b,Al2O3_lattice.c)))

g = g = hkl[0]*Al2O3_lattice_rec[0]+hkl[1]*Al2O3_lattice_rec[1]+hkl[2]*Al2O3_lattice_rec[2]

fi = or_fkt.reziprocal_lattice_Gautam(Al2O3_lattice)