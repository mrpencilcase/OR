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

Al = ["Al" , 0.342,0.2553,1.3240]

fi = or_fkt.intensity_lattice_point(Al,g)