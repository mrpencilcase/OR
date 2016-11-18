# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 09:41:42 2016

@author: Lukas
"""

import numpy as np
import postprocess_lukas as pp
import or_fkt
import math

class Rlp:
    def __init__(self,h,k,l,g,r):
        self.h = h
        self.k = k
        self.l = l 
        self.g = g 
        self.r = r
# paths to the two crystal phases directories need OUTCAR and CONTCAR
pathA = "C:/Users/Lukas/Documents/Diplomarbeit/Results/BaseStructures/alpha_Al2O3"
pathB = "C:/Users/Lukas/Documents/Diplomarbeit/Results/BaseStructures/Cr5Al8_non_mag"
nameA = "alphaAl2O3"
nameB = "Cr5Al8"
supercellA =[1,1,1]
supercellB =[1,1,1]

# Read geometrical datat of the two strucutres
crystalA = pp.read_data(pathA,nameA,supercellA)
aA     = pp.select_data("a",crystalA)
bA     = pp.select_data("b",crystalA)
cA     = pp.select_data("c",crystalA)
alphaA = math.radians(pp.select_data("alpha",crystalA)[0])
betaA  = math.radians(pp.select_data("beta",crystalA)[0])
gammaA = math.radians(pp.select_data("gamma",crystalA)[0])

alA = np.sqrt(sum(np.square(aA)))
blA = np.sqrt(sum(np.square(bA)))
clA = np.sqrt(sum(np.square(cA)))
latticeA_orth = or_fkt.orthon_trans(alA,blA,clA,alphaA,betaA,gammaA)
latticeA_rec = or_fkt.reziprocal_lattice(latticeA_orth)

A_rlp=[]
h_max = or_fkt.hkl_max(latticeA_rec[0],R)
k_max = or_fkt.hkl_max(latticeA_rec[1],R)
l_max = or_fkt.hkl_max(latticeA_rec[2],R)

h=1
while h <= h_max:
    k = 0 
    while k <= k_max:
        l = 0        
        while l <= l_max:
            g = latticeA_rec[0]*h + latticeA_rec[1]*k + latticeA_rec[2]*l
            x = Rlp(h,k,l,g,r)
            if np.sqrt(sum(np.square(g))) <= R: 
                 A_rlp.append(x)          
            l =+ 1
        k =+ 1
    h =+ 1 



crystalB = pp.read_data(pathB,nameB,supercellB)
aB     = pp.select_data("a",crystalB)
bB     = pp.select_data("b",crystalB)
cB     = pp.select_data("c",crystalB)
alphaB = math.radians(pp.select_data("alpha",crystalB)[0])
betaB  = math.radians(pp.select_data("beta",crystalB)[0])
gammaB  = math.radians(pp.select_data("gamma",crystalB)[0])

alB = np.sqrt(sum(np.square(aB)))
blB = np.sqrt(sum(np.square(bB)))
clB = np.sqrt(sum(np.square(cB)))
latticeB_orth = or_fkt.orthon_trans(alB,blB,clB,alphaB,betaB,gammaB)
latticeB_rec = or_fkt.reziprocal_lattice(latticeB_orth)


# set the maximal radius of R and r
# R defines how many lattice planes are looked at
# r sets the volumes around the reciprocal lattice points.
R = 10  * np.sqrt(sum(np.square(latticeA_rec[0])))
r = 0.2 * np.sqrt(sum(np.square(latticeA_rec[0])))

# find the maximal h,k and l values.

B_rlp=[]

    
h=1
while h <= h_max:
    k = 0 
    while k <= k_max:
        l = 0        
        while l <= l_max:
            g = latticeB_rec[0]*h + latticeB_rec[1]*k + latticeB_rec[2]*l
            x = Rlp(h,k,l,g,r)
            if np.sqrt(sum(np.square(g))) <= R: 
                 B_rlp.append(x)          
            l =+ 1
        k =+ 1
    h =+ 1     
    
    
    
# determin all reciprocal lattice points within R for both crystals
# calculate the the overlaping volumes of the reciprocal lattice points
# and summ it up
Vab = 0
for i in A_rlp:
    for j in B_rlp:
        d = np.sqrt(sum(np.square(i.g-j.g)))
        ri = i.r
        rj = j.r
        if d <= ri+rj:
            Vij = np.pi/(12*d)*np.square(ri+rj-d)*(np.square(d)+2*d*(ri+rj)-3*np.square(ri-rj))
            Vab = Vab + Vij
