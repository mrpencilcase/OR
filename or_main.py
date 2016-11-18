# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 09:41:42 2016

@author: Lukas
"""  
import numpy as np
import or_fkt
import math
import os 
import time

def or_main_fkt(settings,matA,matB):
 
    
    # paths to the two crystal phases directories need OUTCAR and CONTCAR
    alpha_max = settings.alpha_max
    alpha_inc = settings.alpha_inc
    beta_max = settings.beta_max
    beta_inc = settings.beta_inc
    gamma_max = settings.gamma_max
    gamma_inc = settings.gamma_inc
    R_scale = settings.R_scale
    r_scale = settings.r_scale
    # Read geometrical datat of the two strucutres
    aA     = matA.a
    bA     = matA.b
    cA     = matA.c
    alphaA = matA.alpha
    betaA  = matA.beta
    gammaA = matA.gamma
    
 
    # calculate the lattice of material A in an orthonormal coordinate system
    # and calculate the reviprcal lattice of the new lattice
    latticeA_orth = or_fkt.orthon_trans(aA,bA,cA,alphaA,betaA,gammaA)
    latticeA_rec = or_fkt.reziprocal_lattice(latticeA_orth)
    
    aB     = matB.a
    bB     = matB.b
    cB     = matB.c
    alphaB = matB.alpha
    betaB  = matB.beta
    gammaB = matB.gamma
    
    # calculate the lattice of material B in an orthonormal coordinate system
    # the reciprokal lattice of B will be calculated later since it is the 
    # one that will be rotated
    latticeB_orth = or_fkt.orthon_trans(aB,bB,cB,alphaB,betaB,gammaB)
    # set the maximal radius of R and r
    # R defines how many lattice planes are looked at
    # r sets the volumes around the reciprocal lattice points.
    R = R_scale * np.sqrt(sum(np.square(latticeA_rec[0])))
    r = r_scale * np.sqrt(sum(np.square(latticeA_rec[0])))
    V_RLP = np.pi*4*np.power(r,3)/3
    
    # find the maximal miller indizes of mat A which are within of R along the 
    # three reciprocal lattice vectors. 
    h_max = or_fkt.hkl_max(latticeA_rec[0],R)
    k_max = or_fkt.hkl_max(latticeA_rec[1],R)
    l_max = or_fkt.hkl_max(latticeA_rec[2],R)   
    # Calculate all RLP's of A that are in an sphere with a radius of R 
    A_rlp=[]
    h=1
    while h <= h_max:
        
        k = 0 
        while k <= k_max:
            l = 0        
            while l <= l_max:
                g = latticeA_rec[0]*h + latticeA_rec[1]*k + latticeA_rec[2]*l
                if np.sqrt(sum(np.square(g))) <= R: 
                    A_rlp.append([h,k,l,g[0],g[1],g[2]])
                l = l+  1
            k =k  + 1
        
        h = h + 1 
    
    
    alpha = 0
    Valpha=[]            
    ri = r
    rj = r
    #start_calc = time.time()
    while alpha <= alpha_max:
        beta =0
        while beta <= beta_max:
            gamma = 0 
            while gamma <= gamma_max:
                    #print("alpha: " + str(np.rad2deg(alpha))+"°")
                    #print("beta: " + str(np.rad2deg(beta))+"°")
                    latticeB_rec = or_fkt.reziprocal_lattice(or_fkt.rot_z(or_fkt.rot_x(latticeB_orth,alpha),beta))
                    # find the maximal h,k and l values.
                    B_rlp=[]       
                  
                    h=1
                    while h <= h_max:
                        k = 0 
                        while k <= k_max:
                            l = 0        
                            while l <= l_max:
                                g = latticeB_rec[0]*h + latticeB_rec[1]*k + latticeB_rec[2]*l
                                if np.sqrt(sum(np.square(g))) <= R: 
                                     B_rlp.append([h,k,l,g[0],g[1],g[2]])             
                                l = l+ 1
                            k = k + 1
                        h = h + 1     
            
                    # determin all reciprocal lattice points within R for both crystals
                    # calculate the the overlaping volumes of the reciprocal lattice points
                    # and summ it up
                    Vab = 0   
        
                    #start = time.time()
                    for entA in A_rlp:
                        for entB in B_rlp:
                            d = np.sqrt(np.square(entA[3]-entB[3])+np.square(entA[4]-entB[4])+np.square(entA[5]-entB[5]))       
                            if d <= ri+rj:
                                Vij = np.pi/(12*d)*np.square(ri+rj-d)*(np.square(d)+2*d*(ri+rj)-3*np.square(ri-rj))
                                Vab = Vab + Vij
                    #print("Time_inc :"+ str(time.time()-start)) 
                    #print("Time_tot:" + str(time.time()-start_calc))
                   # print("Overlapping Volume :" + str(Vab))
                    #print("")    
                    Valpha.append([Vab/V_RLP,np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)])
                    gamma = gamma +gamma_inc
            beta = beta + beta_inc  
        alpha = alpha + alpha_inc
    
    
    
    #print("Run Time: " + str(time.time()-start_calc))
    #print("")
    with open(settings.path_save+"/V_"+matA.name+"_"+matB.name+"_R"+str(R_scale)+"_r"+str(r_scale)+".dat","w") as dat:
        dat.write("# Interface between {} and {}\n".format(matA.name,matB.name))
        dat.write("# Alpha = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(alpha),np.rad2deg(alpha) ))
        dat.write("# Beta = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(beta),np.rad2deg(beta) ))
        dat.write("# R = {} r = {}\n".format(R_scale,r_scale))    
        for ent in  Valpha:
            dat.write("{:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2]))
        