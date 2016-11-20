# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:57:20 2016

@author: Lukas
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import os
# Read data from file
path = "C:/Users/Lukas/Documents/Diplomarbeit/findOR/OrientationalRelations/_CoO_YSZ"

R_select = 3
gamma_select = 0.0

files = os.listdir(path)

for filename in files:
    if "R{}_".format(R_select) in filename:
        alpha = []
        beta = []
        V = []
        
        with open("{}/{}".format(path,filename),"r") as Vfile:
            data = list(Vfile)
            for ent in data:
                if "#" in ent:
                    print(ent)
                else:
                    line = [x for x in ent.split()]
                    if float(line[3]) == gamma_select:            
                        V.append(float(line[0]))
                        alpha.append(float(line[1]))
                        beta.append(float(line[2]))

        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        ax.plot_trisurf(alpha, beta, V, cmap=cm.jet, linewidth=0.2)
        
        plt.show()
