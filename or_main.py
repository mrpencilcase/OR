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
    delta0 = 0.20

    I = np.eye(3)

    """
    calculate the CRLPs and theire coresponding intensits 
    """

    #transform lattice to orthonormal basis
    #TA = or_fkt.orthon_trans(lattA.a, lattA.b, lattA.c, lattA.alpha, lattA.beta, lattA.gamma) 
    #lattice_A = np.transpose(TA)
    #latticeA_rec  = or_fkt.reziprocal_lattice_Gautam(TA)
    #TB = or_fkt.orthon_trans(lattB.a, lattB.b, lattB.c, lattB.alpha, lattB.beta, lattB.gamma)
    #lattice_B = np.transpose(TA)
    #latticeB_rec = or_fkt.reziprocal_lattice_Gautam(TB)
    #calculate atom positons in new basis
    latticeA = np.vstack((lattA.a,lattA.b,lattA.c))
    latticeB = np.vstack((lattB.a,lattB.b,lattB.c))
    latticeA_rec = or_fkt.reziprocal_lattice(latticeA) 
    latticeB_rec = or_fkt.reziprocal_lattice(latticeB)             
    map_intens = []
   

    """
    Determin the larger reciprocal structure of the two materials and apply the given hkl values to it
    for the smaller structure the maximal hkl values are choosen so that roughly the same reciprpcal 
    space is coverd.
    """
    hkl_a_max, hkl_b_max = or_fkt.get_hkl(latticeA_rec,latticeB_rec,hkl)

    intensA, gA, hkl_a = or_fkt.calc_intensities(hkl_a_max, latticeA_rec, unitA,I)
    intensB, gB, hkl_b = or_fkt.calc_intensities(hkl_b_max, latticeB_rec, unitB,I)
      
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
                    tot_intens,dummy = or_fkt.overlap_lattices(intensA,gA,hkl_a,intensB,gB_rot,hkl_b,ImaxA,ImaxB,delta0)
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
    if not os.path.exists("{}OR".format(settings.path_save)):
            os.makedirs("{}OR".format(settings.path_save))

    path_save = "{}/OR/{}.dat".format(settings.path_save,file_name)
    with open(path_save,"w") as dat:
        dat.write("# Interface between {} and {}\n".format(lattA.name,lattB.name))
        dat.write("# Alpha = {:.1f}-{:.1f}? in {:.1f}? increments\n".format(np.rad2deg(alpha_start),np.rad2deg(alpha_end),np.rad2deg(alpha_inc) ))
        dat.write("# Beta = {:.1f}-{:.1f}? in {:.1f}? increments\n".format(np.rad2deg(beta_start),np.rad2deg(beta_end),np.rad2deg(beta_inc) ))
        dat.write("# Gamma = {:.1f}-{:.1f}? in {:.1f}? increments\n".format(np.rad2deg(gamma_start),np.rad2deg(gamma_end),np.rad2deg(gamma_inc) ))
        dat.write("# delta0 = {:.2f}\n".format(delta0))
        dat.write("# HKL up to = {} {} {}\n".format(hkl[0],hkl[1],hkl[2]))
        dat.write("# Calculation Time: {:5.1f} sec \n".format(time.time()-start))
        
        for ent in  map_intens:
            dat.write("{:f} {:f} {:f} {:f}\n".format(ent[0],ent[1],ent[2],ent[3]))

    return path_save, file_name

def io_gautam_meth(settings,angles,lattA,unitA,lattB,unitB,hkl):
    """
    The function looks into the different sett of angels and returns the individual overlabpping
    intensities of the spots.  
    """
    io_data = []

    I = np.eye(3)


    #TA = or_fkt.orthon_trans(lattA.a, lattA.b, lattA.c, lattA.alpha, lattA.beta, lattA.gamma) 
    #lattice_A = np.transpose(TA)
    #latticeA_rec  = or_fkt.reziprocal_lattice_Gautam(TA)
    #TB = or_fkt.orthon_trans(lattB.a, lattB.b, lattB.c, lattB.alpha, lattB.beta, lattB.gamma)
    #lattice_B = np.transpose(TA)
    #latticeB_rec = or_fkt.reziprocal_lattice_Gautam(TB)
    

    #calculate atom positons in new basis
    latticeA = np.vstack((lattA.a,lattA.b,lattA.c))
    latticeB = np.vstack((lattB.a,lattB.b,lattB.c))
    latticeA_rec = or_fkt.reziprocal_lattice(latticeA) 
    latticeB_rec = or_fkt.reziprocal_lattice(latticeB)  
    
    hkl_a_max, hkl_b_max = or_fkt.get_hkl(latticeA_rec,latticeB_rec,hkl)

    intensA, gA, hkl_a= or_fkt.calc_intensities(hkl_a_max, latticeA_rec, unitA,I)
    intensB, gB, hkl_b = or_fkt.calc_intensities(hkl_b_max, latticeB_rec, unitB,I)
      
    #Find the maximal Intensities in A and B. 
    ImaxA = max(intensA)
    ImaxB = max(intensB)   
    delta0 = 0.20
    for ent in angles:
        gB_rot = or_fkt.rot_z(or_fkt.rot_x(or_fkt.rot_z(gB,np.deg2rad(ent[3])),np.deg2rad(ent[2])),np.deg2rad(ent[1]))
        tot_intens, singel_intens = or_fkt.overlap_lattices(intensA,gA,hkl_a,intensB,gB_rot,hkl_b,ImaxA,ImaxB,delta0)
        io_data.append(singel_intens)
    
    if not os.path.exists("{}/IO".format(settings.path_save)):
            os.makedirs("{}/IO".format(settings.path_save))

    with open("{}/IO/{}.dat".format(settings.path_save,settings.name),"w") as file:
        for io,ang in zip(io_data,angles):
            file.write("#Total Intensity: {}\n".format(ang[0]))
            file.write("#Angles: {} {} {}\n".format(ang[1],ang[2],ang[3]))
            for ent in io:
                #angel_gA_gB = np.rad2deg(or_fkt.angle_between(ent[1],ent[3]))
                file.write("({:2d},{:2d},{:2d})||({:2d},{:2d},{:2d})   {:8.1f}\n".format(ent[2][0],ent[2][1],ent[2][2],ent[4][0],ent[4][1],ent[4][2],ent[0]))
    
    return tot_intens