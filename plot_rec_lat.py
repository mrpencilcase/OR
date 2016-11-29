# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 09:56:26 2016

@author: Lukas
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import os
import or_fkt



def plot_rec(lattice):
    
    hkl =[3,3,3]    
    

    lattice_rec = or_fkt.reziprocal_lattice(lattice)
    
    
    
    x = []
    y = []
    labels = []    
    h = -hkl[0]
    while h <= hkl[0]:
        
        k = -hkl[1]
        while k <= hkl[1]:
            l = -hkl[2]
            while l <= hkl[2]:
                g = g = lattice_rec[0]*h + lattice_rec[1]*k + lattice_rec[2]*l
                x.append(g[0])
                y.append(g[1])
                labels.append("{} {} {}".format(h,k,l))                
                l += 1 
            k += 1
        h += 1
    plt.plot(x,y,"o")
    for label,xt,yt in zip(labels,x,y):
        plt.annotate(
        label, 
        xy=(xt,yt), xytext = (-20, 20))
    
    
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

lattice = np.vstack((testMatB.a,testMatB.b,testMatB.c))
lattice_rot = or_fkt.rot_z(lattice,np.deg2rad(30))
plot_rec(lattice_rot)