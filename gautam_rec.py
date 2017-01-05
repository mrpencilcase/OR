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

def or_gautam_meth(settings,lattA,unitA,lattB,unitB, hkl):

    """
    Find OR of two materials following the methode proposed by Gautam and Howe.
    Input:    
    Settings:   provides the parameters for the rotaion of material B aswell as informations
                where to save the results.
    lattA/B:    provides the parameters of the lattice of material A and B .
    unitA/B:    contains the atomes and theire positions within the unit cell.
    hkl: 
    Output:
    The function itself does not return values but creates a file in which the overlapping 
    intensitys, and the corresponding angles are listed.
    """
    alpha_start = settings.alpha_start
    alpha_end = settings.alpha_end
    alpha_inc = settings.alpha_inc
    beta_start = settings.beta_start
    beta_end = settings.beta_end
    beta_inc = settings.beta_inc
    gamma_start = settings.gamma_start
    gamma_end = settings.gamma_end
    gamma_inc = settings.gamma_inc


    """
    calculate the CRLPs and theire coresponding intensits 
    """
    #transform lattice to orthonormal basis
    #latticeA_orth = or_fkt.orthon_trans(lattA.a,lattA.b,lattA.c,lattA.alpha,lattA.beta,lattA.gamma) 
    #latticeA_rec  = or_fkt.reziprocal_lattice_Gautam(latticeA_orth)
    #calculate atom positons in new basis
    latticeA = np.vstack((lattA.a,lattA.b,lattA.c))
    latticeA_rec = or_fkt.reziprocal_lattice(latticeA)
    latticeB = np.vstack((lattB.a,lattB.b,lattB.c))
    latticeB_rec = or_fkt.reziprocal_lattice(latticeB)             
    map_intens = []
    hkl_a = []
    hkl_b = []

    """
    Determin the larger reciprocal structure of the two materials and apply the given hkl values to it
    for the smaller structure the maximal hkl values are choosen so that roughly the same reciprpcal 
    space is coverd.
    """
    if la.norm(latticeA_rec[0,:]) <= la.norm(latticeB_rec[0,:]):
       
       hkl_a.append(int(hkl[0] * la.norm(latticeB_rec[0,:]) / la.norm(latticeA_rec[0,:])))
       hkl_a.append(int(hkl[1] * la.norm(latticeB_rec[1,:]) / la.norm(latticeA_rec[1,:])))
       hkl_a.append(int(hkl[2] * la.norm(latticeB_rec[2,:]) / la.norm(latticeA_rec[2,:])))
       
       hkl_b.append(hkl[0])
       hkl_b.append(hkl[1])
       hkl_b.append(hkl[2])
       
         
    elif la.norm(latticeB_rec[0,:]) < la.norm(latticeA_rec[0,:]):

       hkl_b.append(int(hkl[0] * la.norm(latticeA_rec[0,:]) / la.norm(latticeB_rec[0,:])))
       hkl_b.append(int(hkl[1] * la.norm(latticeA_rec[1,:]) / la.norm(latticeB_rec[1,:])))
       hkl_b.append(int(hkl[2] * la.norm(latticeA_rec[2,:]) / la.norm(latticeB_rec[2,:])))

       hkl_a.append(hkl[0])
       hkl_a.append(hkl[1])
       hkl_a.append(hkl[2])
       
    
    delta0 = 0.20
    intensA = []
    gA = []
    rlpA = []
    
    intensB=[]                       
    gB = []
    rlpB = []
    h = -hkl_a[0]
    while h <= hkl_a[0]:
        gh = latticeA_rec[0] * h       
        k = -hkl_a[1]    
        while k <= hkl_a[1]:
            ghk = latticeA_rec[1] * k + gh
            l = -hkl_a[2]
            while l <= hkl_a[2]:
                if sum(map(abs,[h,k,l])) > 0: 
                    g = latticeA_rec[2] * l + ghk
                    #intensA.append(or_fkt.intensity_lattice_point1(unitA,g))
                    intensA.append(or_fkt.intensity_lattice_point2(unitA,[h,k,l],g))
                    gA.append(g)        
                l += 1                
            k += 1
        h += 1

    h=-hkl_b[0]    
    while h <= hkl_b[0]:
        gh = latticeB_rec[0] * h       
        k = -hkl_b[1]    
        while k <= hkl_b[1]:
            ghk = latticeB_rec[1] * k + gh
            l = -hkl_b[2]
            while l <= hkl_b[2]:
                if sum(map(abs,[h,k,l])) > 0: 
                    g = latticeB_rec[2] * l + ghk
                    #intensB.append(or_fkt.intensity_lattice_point1(unitB,g))
                    intensB.append(or_fkt.intensity_lattice_point2(unitB,[h,k,l],g))
                    gB.append(g)
                l += 1                
            k += 1
        h += 1

    #Find the maximal Intensities in A and B. 
    ImaxA = max(intensA)
    ImaxB = max(intensB)
    """
    Rotate material B in the set increments and calculate the overlaping 
    intenseties for each iteration
    
    """

    alpha = alpha_start
    start = time.time()
    while alpha <= alpha_end:
        beta =beta_start
        while beta <= beta_end:
            
            start_inc = time.time()            
            gamma = gamma_start 
            #start_time = time.time()
            while gamma <= gamma_end:
                          
                    gB_rot = or_fkt.rot_z(or_fkt.rot_x(or_fkt.rot_z(gB,gamma),beta),alpha)
                    tot_intens = 0   
                    numb_overl = 0 
                    #start = time.time()
                    for intAi, gAi in zip(intensA,gA):
                        for intBj, gBj in zip(intensB,gB_rot):
                                                   
                            d = la.norm(gAi-gBj)
                            d1 = 1/(la.norm(gAi))
                            p = delta0/(6*d1*(ImaxA+ImaxB))
                            Ri = 6*intAi*p
                            Rj = 6*intBj*p
                            HWi = intAi*p
                            HWj = intBj*p
                            #Only overlapping spheres are from interest
                            if d < Ri + Rj:
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

 
                    map_intens.append([tot_intens,np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)])
                    gamma = gamma +gamma_inc
                    
            print("alpha: {:5.1f}   beta: {:5.1f}  intensity:{:10.1f}  time: {:5.1f}".format(np.rad2deg(alpha),np.rad2deg(beta),tot_intens,time.time()-start))
            beta = beta + beta_inc  
        alpha = alpha + alpha_inc
    
    
    
    #print("Run Time: " + str(time.time()-start_calc))
    #print("")
    #path_save = settings.path_save+lattA.name+"_"+lattB.name + "_HKL_"+str(hkl[0])+str(hkl[1])+str(hkl[2])+"_d0_"+str(delta0)+ ".dat"
    file_name = "{}_{}_HKL_{}{}{}_d0_{}_a{}_{}_b{}_{}".format(lattA.name, lattB.name,
                                                          hkl[0],hkl[1],hkl[2],delta0,
                                                          np.rad2deg(alpha_start),np.rad2deg(alpha_end),
                                                          np.rad2deg(beta_start),np.rad2deg(beta_end))

    path_save = "{}{}.dat".format(settings.path_save,file_name)
    with open(path_save,"w") as dat:
        dat.write("# Interface between {} and {}\n".format(lattA.name,lattB.name))
        dat.write("# Alpha = {:.1f}-{:.1f}� in {:.1f}� increments\n".format(np.rad2deg(alpha_start),np.rad2deg(alpha_end),np.rad2deg(alpha_inc) ))
        dat.write("# Beta = {:.1f}-{:.1f}� in {:.1f}� increments\n".format(np.rad2deg(beta_start),np.rad2deg(beta_end),np.rad2deg(beta_inc) ))
        dat.write("# Gamma = {:.1f}-{:.1f}� in {:.1f}� increments\n".format(np.rad2deg(gamma_start),np.rad2deg(gamma_end),np.rad2deg(gamma_inc) ))
        dat.write("# delta0 = {:.2f}\n".format(delta0))
        #dat.write("# R = {} r = {}\n".format(R_scale,r_scale))    
        dat.write("# HKL up to = {} {} {}\n".format(hkl[0],hkl[1],hkl[2]))
        dat.write("# Calculation Time: {:5.1f} sec \n".format(time.time()-start))
        for ent in  map_intens:
            dat.write("{:f} {:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2],ent[3]))

    return path_save, file_name
        