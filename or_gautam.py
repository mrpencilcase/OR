# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 11:10:35 2016

@author: Lukas
"""

import numpy as np
from numpy import linalg as la
import or_fkt
import math
import os 
import time

def or_gautam_meth(settings,lattA,unitA,lattB,unitB, HKL):

    """
    Find OR of two materials following the methode proposed by Gautam and Howe.
    Input:    
    Settings:   provides the parameters for the rotaion of material B aswell as informations
                where to save the results.
    lattA/B:    provides the parameters of the lattice of material A and B .
    unitA/B:    contains the atomes and theire positions within the unit cell.
    Output:
    The function itself does not return values but creates a file in which the overlapping 
    intensitys, and the corresponding angles are listed.
    """
    alpha_max = settings.alpha_max
    alpha_inc = settings.alpha_inc
    beta_max = settings.beta_max
    beta_inc = settings.beta_inc
    gamma_max = settings.gamma_max
    gamma_inc = settings.gamma_inc


    """
    MATERIAL A: calculate the CRLPs and theire coresponding intensits 
    """
    # transform lattice to orthonormal basis
    #latticeA_orth = or_fkt.orthon_trans(lattA.a,lattA.b,lattA.c,lattA.alpha,lattA.beta,lattA.gamma) 
    #latticeA_rec  = or_fkt.reziprocal_lattice_Gautam(latticeA_orth)
    # calculate atom positons in new basis
    latticeA_rec = or_fkt.reziprocal_lattice(np.vstack((lattA.a,lattA.b,lattA.c)))
    latticeB = np.vstack((lattB.a,lattB.b,lattB.c))
                                                        
    map_intens = []
    hkl = HKL
    delta0 = 0.2
    intensA = []
    gA = []
    h = -hkl[0]
    rlpA = []
    while h <= hkl[0]:
        gh = latticeA_rec[0] * h       
        k = -hkl[1]    
        while k <= hkl[1]:
            ghk = latticeA_rec[1] * k * gh
            l = -hkl[2]
            while l <= hkl[2]:
                if sum(map(abs,[h,k,l])) > 0: 
                    g = latticeA_rec[2] * l + ghk
                    intensA.append(or_fkt.intensity_lattice_point(unitA,g))
                    gA.append(g)        
                l += 1                
            k += 1
        h += 1
    
    """
    Rotate material B in the set increments and calculate the overlaping 
    intenseties for each iteration
    
    """
    alpha = 0
    latticeB_orth = or_fkt.orthon_trans(lattB.a,lattB.b,lattB.c,lattB.alpha,lattB.beta,lattB.gamma)
    start = time.time()
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
                    intensB=[]                       
                    gB = []
                    h=-hkl[0]
                    while h <= hkl[0]:
                        gh = latticeB_rec[0] * h       
                        k = -hkl[1]    
                        while k <= hkl[1]:
                            ghk = latticeB_rec[1] * k * gh
                            l = -hkl[2]
                            while l <= hkl[2]:
                                if sum(map(abs,[h,k,l])) > 0: 
                                    g = latticeB_rec[2] * l + ghk
                                    intensB.append(or_fkt.intensity_lattice_point(unitA,g))
                                    gB.append(g)
                                l += 1                
                            k += 1
                        h += 1
                    
                    #Find the maximal Intensities in A and B.       
                    ImaxA = max(intensA)
                    ImaxB = max(intensB)

                    # determin all reciprocal lattice points within R for both crystals
                    # calculate the the overlaping volumes of the reciprocal lattice points
                    # and summ it up
                    tot_intens = 0   
                    numb_overl = 0 
                    #start = time.time()
                    for intAi, gAi in zip(intensA,gA):
                        for intBj, gBj in zip(intensB,gB):
                                                   
                            d = la.norm(gAi-gBj)
                            d1 = 1/(la.norm(gAi))
                            p = delta0/(6*d1*(ImaxA+ImaxB))
                            Ri = 6*intAi*p
                            Rj = 6*intBj*p
                            HWi = intAi*p
                            HWj = intBj*p
                            #Only overlapping spheres are from interest
                            if d < Ri+Rj :
                                numb_overl += 1 
                                if d + Rj <= Ri:
                                    I1 = or_fkt.intensity_overlap_tot(intAi,Rj,d,HWi)
                                    I2 = intBj
                                    tot_intens += I1 + I2
                                elif d + Ri <= Rj:
                                    I1 = or_fkt.intensity_overlap_tot(intBj,Ri,d,HWj)
                                    I2 = intAi
                                    tot_intens += I1 + I2
                                else:
                                    I1 = or_fkt.intensity_overlap_part(intAi,Ri,Rj,d,HWi)
                                    I2 = or_fkt.intensity_overlap_part(intBj,Rj,Ri,d,HWj)
                                    tot_intens += I1 + I2

 
                    map_intens.append([tot_intens,np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma), numb_overl])
                    gamma = gamma +gamma_inc
                    
            print("alpha: {:5.1f}   beta: {:5.1f}  time: {:5.1f}".format(np.rad2deg(alpha),np.rad2deg(beta),time.time()-start))
            beta = beta + beta_inc  
        alpha = alpha + alpha_inc
    
    
    
    #print("Run Time: " + str(time.time()-start_calc))
    #print("")
    path_save = settings.path_save+lattA.name+"_"+lattB.name + "_HKL_"+str(hkl[0])+str(hkl[1])+str(hkl[2])+"_d0_"+str(delta0)+ ".dat"
    with open(path_save,"w") as dat:
        dat.write("# Interface between {} and {}\n".format(lattA.name,lattB.name))
        dat.write("# Alpha = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(alpha_max),np.rad2deg(alpha_inc) ))
        dat.write("# Beta = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(beta_max),np.rad2deg(beta_inc) ))
        dat.write("# Gamma = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(gamma_max),np.rad2deg(gamma_inc) ))
        dat.write("# delta0 = {:.2f}\n".format(delta0))
        #dat.write("# R = {} r = {}\n".format(R_scale,r_scale))    
        dat.write("# HKL up to = {} {} {}\n".format(hkl[0],hkl[1],hkl[2]))
        for ent in  map_intens:
            dat.write("{:f} {:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2],ent[3]))

    return path_save
        