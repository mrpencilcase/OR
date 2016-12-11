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

def or_gautam_meth(settings,lattA,unitA,lattB,unitB):

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
    latticeA_orth = or_fkt.orthon_trans(lattA.a,lattA.b,lattA.c,lattA.alpha,lattA.beta,lattA.gamma) 

    latticeA_rec  = or_fkt.reziprocal_lattice_Gautam(latticeA_orth)
    # calculate atom positons in new basis
    map_intens = []
    hkl = [2,2,2]
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
    while alpha <= alpha_max:
        beta =0
        while beta <= beta_max:
            
            start_inc = time.time()            
            gamma = 0 
            #start_time = time.time()
            while gamma <= gamma_max:
                    #print("alpha: " + str(np.rad2deg(alpha))+"°")
                    #print("gamma   : {:.1f}°".format(np.rad2deg(gamma)))

                    latticeB_rec = or_fkt.reziprocal_lattice_Gautam(or_fkt.rot_z(or_fkt.rot_x(or_fkt.rot_z(latticeB_orth,gamma),beta),alpha))
    
    
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
                   
                    #start = time.time()
                    for intAi, gAi in zip(intensA,gA):
                        for intBj, gBj in zip(intensB,gB):
                                                   
                            d = la.norm(gAi-gBj)
                            d1 = 1/(la.norm(gAi))
                            p = delta0/(6*d1(ImaxA+ImaxB))
                            Ri = 6*intAi*p
                            Rj = 6*intBj*p
                            #Only overlapping spheres are from interest
                            if d < Ri+Rj :
                                if d + Rj <= Ri:
                                    I1 = or_fkt.intensity_oberlap_tot()
                                    I2 = 2
                                    tot_intens += I1 + I2
                                elif d + Ri <= Rj:
                                    I1 = or_fkt.intensity_overlap_tot()
                                    I2 = 2
                                    tot_intens += I1 + I2
                                else:
                                    I1 = or_fkt.intensity_overlap_part()
                                    I2 = or_fkt.intensity_overlap_part()
                                    tot_intens += I1 + I2

 
                    map_intens.append([Vab,np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)])
                    gamma = gamma +gamma_inc
                    
            print("R: {}   r: {:.2f}   alpha: {:5.1f}   beta: {:5.1f}   time_tot: {:.2f}   time_inc: {:.2f}".format(R_scale,r_scale,np.rad2deg(alpha),np.rad2deg(beta),time.time()-start_calc,time.time()-start_inc))
            beta = beta + beta_inc  
        alpha = alpha + alpha_inc
    
    
    
    #print("Run Time: " + str(time.time()-start_calc))
    #print("")
    with open(settings.path_save+"/V_"+matA.name+"_"+matB.name+"_R"+str(R_scale)+"_r"+str(r_scale)+".dat","w") as dat:
        dat.write("# Interface between {} and {}\n".format(matA.name,matB.name))
        dat.write("# Alpha = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(alpha_max),np.rad2deg(alpha_inc) ))
        dat.write("# Beta = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(beta_max),np.rad2deg(beta_inc) ))
        dat.write("# Gamma = 0-{:.1f}° in {:.1f}° increments\n".format(np.rad2deg(gamma_max),np.rad2deg(gamma_inc) ))
        dat.write("# R = {} r = {}\n".format(R_scale,r_scale))    
        for ent in  map_intens:
            dat.write("{:f} {:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2],ent[3]))
        