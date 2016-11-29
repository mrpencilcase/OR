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


path = "C:\Users\Lukas\Documents\Diplomarbeit\Results\OrientationalRelations\CoO_YSZ"

R_select = 3 
r_select = [0.1,0.15,0.2,0.25,0.3,0.35]
gamma_select = 0.0
beta_select = 0.0

files = os.listdir(path)

for filename in files:
    if "R{}_".format(R_select) in filename:
        alpha = []
        beta = []
        V = []
        
        with open("{}/{}".format(path,filename),"r") as Vfile:
            data = list(Vfile)
            matA =data[0].split()[3]
            matB =data[0].split()[5]
            R = data[4].split()[3]
            r = data[4].split()[6]
            for ent in data:
                if "#" in ent:
                    
                    print(ent.rstrip())
                else:
                    line = [x for x in ent.split()]
                    if float(line[3]) == gamma_select:            
                        V.append(float(line[0]))
                        alpha.append(float(line[1]))
                        beta.append(float(line[2]))
        
                
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        ax.plot_trisurf(alpha, beta, V, cmap=cm.jet, linewidth=0.2)
        ax.set_title("OR between {} and {} (R={}, r={})".format(matA,matB,R,r))
        ax.set_xlabel("Alpha")
        ax.set_ylabel("Beta")
        ax.set_zlabel("Overlaping Volume")
        
#        plt.title("".format(R,r))
#        plt.suptitle("OR between {} and {}", y=1.05, fontsize=18)
        
        if not os.path.exists("{}/plots".format(path)):
                os.makedirs("{}/plots".format(path))
                
        fig.savefig("{}/plots/{}_{}_R{}_r{}.pdf".format(path,matA,matB,R,r))        
        plt.close(fig)
#        fig = plt.figure()
#        
#        plt.plot(alpha,V)
#        plt.show()
