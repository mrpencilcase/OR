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

def plot_intens(path,gamma_select,select):

   
    alpha = []
    beta = []
    V = []
    if select == 1:     
        with open(path,"r") as Vfile:
            data = list(Vfile)
            matA = data[0].split()[3]
            matB = data[0].split()[5]
            delta0 = data[4].split()[3] 
            hkl = [data[5].split()[5],data[5].split()[6],data[5].split()[7]]
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
        ax.set_title("OR between {} and {})".format(matA,matB))
        ax.set_xlabel("Alpha")
        ax.set_ylabel("Beta")
        ax.set_zlabel("Overlaping Volume")
        
        #        plt.title("".format(R,r))
        #        plt.suptitle("OR between {} and {}", y=1.05, fontsize=18)
    
        path_dir,filename = os.path.split(path)
        if not os.path.exists("{}/plots".format(path_dir)):
                os.makedirs("{}/plots".format(path_dir))
                
        fig.savefig("{}/plots/{}_{}_{}_{}{}{}.pdf".format(path_dir,matA,matB,delta0,hkl[0],hkl[1],hkl[2]))    
        plt.close(fig)

    elif select == 2:
        with open(path,"r") as Vfile:
            data = list(Vfile)
            matA = data[0].split()[3]
            matB = data[0].split()[5]
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
        ax.set_title("OR between {} and {} with R={:.1f} and r={:.2f})".format(matA,matB,float(R),float(r)))
        ax.set_xlabel("Alpha")
        ax.set_ylabel("Beta")
        ax.set_zlabel("Overlaping Volume")
        
        #        plt.title("".format(R,r))
        #        plt.suptitle("OR between {} and {}", y=1.05, fontsize=18)
    
        path_dir,filename = os.path.split(path)
        if not os.path.exists("{}/plots".format(path_dir)):
                os.makedirs("{}/plots".format(path_dir))
                
        fig.savefig("{}/plots/{}_{}_R{}_r{}.pdf".format(path_dir,matA,matB,float(R),float(r)))    
        plt.close(fig)