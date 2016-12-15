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
 
    # calculate the lattice of material A in an orthonormal coordinate system
    # and calculate the reviprcal lattice of the new lattice

    latticeA = np.vstack((matA.a, matA.b,matA.c))
    latticeA_rec = or_fkt.reziprocal_lattice(latticeA)    
    latticeB = np.vstack((matB.a, matB.b, matB.c))
    
    # calculate the lattice of material B in an orthonormal coordinate system
    # the reciprokal lattice of B will be calculated later since it is the 
    # one that will be rotated
    
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
    h=-h_max
    while h <= h_max:
        
        k = -k_max 
        while k <= k_max:
            l = -l_max        
            while l <= l_max:
                g = latticeA_rec[0]*h + latticeA_rec[1]*k + latticeA_rec[2]*l
                gnorm = np.sqrt(sum(np.square(g)))                
                if  gnorm <= R and gnorm > 0: 
                    A_rlp.append([g[0],g[1],g[2]])
                l = l+  1
            k =k  + 1
        h = h + 1 
    
    
    Valpha=[]            
    ri = r
    rj = r
    r2= 2*r 
    pi12 = np.pi/12
    start_calc = time.time()
    alpha = 0
    while alpha <= alpha_max:
        beta =0
        while beta <= beta_max:
            
            start_inc = time.time()            
            gamma = 0 
            #start_time = time.time()
            while gamma <= gamma_max:
                    #print("alpha: " + str(np.rad2deg(alpha))+"°")
                    #print("gamma   : {:.1f}°".format(np.rad2deg(gamma)))
                
                        
                    latticeB_rec = or_fkt.reziprocal_lattice(or_fkt.rot_z(or_fkt.rot_x(or_fkt.rot_z(latticeB,gamma),beta),alpha))
                    # find the maximal h,k and l values.
                    B_rlp=[]       
                    h_max = or_fkt.hkl_max(latticeB_rec[0],R)
                    k_max = or_fkt.hkl_max(latticeB_rec[1],R)
                    l_max = or_fkt.hkl_max(latticeB_rec[2],R)   
                  
                  
                  
                    h=-h_max
                    while h <= h_max:
                        k = -k_max 
                        gh = latticeB_rec[0]*h
                        while k <= k_max:
                            l = -l_max        
                            gk = gh + latticeA_rec[1]*k                            
                            while l <= l_max:
                                g =  gk + latticeB_rec[2]*l
                                gnorm = np.sqrt(sum(np.square(g)))                
                                if  gnorm <= R and gnorm > 0:  
                                     B_rlp.append([g[0],g[1],g[2]])             
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
                            d = (entA[0]-entB[0])*(entA[0]-entB[0])
                            d = d + (entA[1]-entB[1])*(entA[1]-entB[1])
                            d = d + (entA[2]-entB[2])* (entA[2]-entB[2])
                            d = np.sqrt(d)
                            # d = np.sqrt(np.square(entA[0]-entB[0])+np.square(entA[1]-entB[1])+np.square(entA[2]-entB[2]))       
                            
                            
                            if d <= r2 and d > 0:
                                # Vij = np.pi/(12*d)*np.square(ri+rj-d)*(np.square(d)+2*d*(ri+rj)-3*np.square(ri-rj))
                                Vab = Vab + (pi12/d)*(r2-d)*(r2-d) * (d*d + 2*d*r2)
                            elif d == 0:
                                Vab += V_RLP
                    # print("Time_inc :"+ str(time.time()-start)) 
                    # print("Time_tot:" + str(time.time()-start_calc))
                    # print("Overlapping Volume :" + str(Vab))
                    # print("")    
                    Valpha.append([Vab,np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)])
                    gamma = gamma +gamma_inc
                    
            print("R: {}   r: {:.2f}   alpha: {:5.1f}   beta: {:5.1f}   time_tot: {:.2f}   time_inc: {:.2f}".format(R_scale,r_scale,np.rad2deg(alpha),np.rad2deg(beta),time.time()-start_calc,time.time()-start_inc))
            beta = beta + beta_inc  
        alpha = alpha + alpha_inc
    
    
    
    #print("Run Time: " + str(time.time()-start_calc))
    #print("")
    path_save = settings.path_save+"/V_"+matA.name+"_"+matB.name+"_R"+str(R_scale)+"_r"+str(r_scale)+".dat"
    with open(path_save,"w") as dat:
        dat.write("# Interface between {} and {}\n".format(matA.name,matB.name))
        dat.write("# Alpha = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(alpha_max),np.rad2deg(alpha_inc) ))
        dat.write("# Beta = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(beta_max),np.rad2deg(beta_inc) ))
        dat.write("# Gamma = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(gamma_max),np.rad2deg(gamma_inc) ))
        dat.write("# R = {} r = {}\n".format(R_scale,r_scale))    
        for ent in  Valpha:
            dat.write("{:f} {:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2],ent[3]))
    return path_save