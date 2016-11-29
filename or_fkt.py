# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 11:38:18 2016

@author: Lukas

Functions used to find the most probable orientational relashionships of an 
hetrocrystal interface. The underlying methode is the Reciprocl Lattice Point
methote (RLP) by Ikuhara and Priouz.

"""
import numpy as np
import string
import cmath
from numpy import linalg as la
from math import sqrt
#function to transorm the coordinats of the lattice in an orhtnormal coordinate
#system. The function requires the lattice parameters a,b,c and the three angle
#alpha, beta and gamma in [rad]
def orthon_trans(a,b,c,alpha,beta,gamma):
    
    t11 = a
    t12 = b*np.cos(gamma)
    t13 = c* np.cos(gamma) 
    t21 = 0
    t22 = b*np.sin(gamma)
    t23 = c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma))
    t31 = 0
    t32 = 0
    t33 = np.power(np.cos(alpha),2) + np.power(np.cos(beta),2) - 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) 
    t33 = t33 /np.power(np.sin(gamma),2)
    t33 = c * np.sqrt(1-t33)
    
    t = np.array([ [t11, t12, t13],
                   [t21, t22, t23],
                   [t31, t32, t33] ])
          
    return t
    
    
# Calculate the reziprokal lattice vectors a*, b* and c* by inverting orhonomral basis 
def reziprocal_lattice_Gautam(lattice):    
    
    rl = np.linalg.pinv(lattice)
    
    return rl

def reziprocal_lattice(lattice_real):
    
    a1 = lattice_real[0,:]
    a2 = lattice_real[1,:]
    a3 = lattice_real[2,:]
    Ve = np.dot(a1,np.cross(a2,a3))

    g1 = 2 * np.pi * (np.cross(a2,a3))/Ve
    g2 = 2 * np.pi * (np.cross(a3,a1))/Ve
    g3 = 2 * np.pi * (np.cross(a1,a2))/Ve

    lattice_rec = np.vstack((g1,g2,g3))

    return lattice_rec
    
# Find the last reziprocal points miller indizes that lies within the the given radius R
# along the reziprokal vector rezVec.
def hkl_max(rezVec,R):
    miller = 1
    g =np.sqrt(np.square(rezVec[0]*miller)+np.square(rezVec[1]*miller)+np.square(rezVec[2]*miller))
    while g <= R:
        miller = miller + 1
        g =np.sqrt(np.square(rezVec[0]*miller)+np.square(rezVec[1]*miller)+np.square(rezVec[2]*miller))
    return miller -1
    
     
# Function to rotate lattice around the [1,0,0] axis
def rot_x(lattice,alpha):
    c = np.cos(alpha)
    s = np.sin(alpha)
    
    rot_m = np.array([ [ c, s, 0] ,
                       [-s, c, 0] ,
                       [  0,0, 1] ])
    
    lattice_rot = np.dot(lattice,rot_m)
    return lattice_rot    
# Function to rotate lattice around the [0,0,1] axis
def rot_z(lattice,alpha):
    c = np.cos(alpha)
    s = np.sin(alpha)
    
    rot_m = np.array([ [ 1, 0, 0] ,
                       [ 0, c, s] ,
                       [ 0,-s, c] ])
    
    lattice_rot = np.dot(lattice,rot_m)
    return lattice_rot
    
    
    
def read_data(path,supercell):
 
    #read structure parameters from the CONTCAR
    CONTCAR=open(path +'/CONTCAR','r')
    CONTCAR.readline()
    aSc=float(CONTCAR.readline())
    a1In=[aSc*float(x) for x in CONTCAR.readline().split()]
    a2In=[aSc*float(x) for x in CONTCAR.readline().split()]
    a3In=[aSc*float(x) for x in CONTCAR.readline().split()]    
    CONTCAR.close()
    # process lattice parameters
    data = lattice()
    data.a = a1In
    data.b = a2In
    data.c = a3In
    data.anorm = la.norm(a1In)/supercell[0]
    data.bnorm = la.norm(a2In)/supercell[1]
    data.cnorm = la.norm(a3In)/supercell[2]
    data.alpha = (np.degrees(np.arccos(np.dot(a1In,a2In)/(la.norm(a1In)*la.norm(a2In)))))
    data.beta = (np.degrees(np.arccos(np.dot(a1In,a3In)/(la.norm(a1In)*la.norm(a3In)))))
    data.gamma = (np.degrees(np.arccos(np.dot(a2In,a3In)/(la.norm(a2In)*la.norm(a3In)))))

    return data
    
# function to select a specific entry form the data set provided by read_data    
def select_data(signature,data):
    for ent in data:
        if signature+":" in ent:
            sdata = [x for x in ent.split()]
            sdata.pop(0)           
            sdata_float = [float(x) for x in sdata ]            
            np.asarray(sdata_float)            
            return sdata_float
            
def def_lattice(name,a,b,c):
    latt = lattice()
    latt.a=a
    latt.b=b
    latt.c=c
    return latt
    
def intensity_lattice_point(unit_cell,g):
    g_norm = la.norm(g)
    
    I1 = 0
    I2 = 0  
    pii = 2j*np.pi      
    for ent in unit_cell:
        fj = at_scat_factr(ent[0],gnorm)
        rj = np.vstack((ent[1],ent[2],ent[3]))
        alpha = rj*g*pii
        I1 += fj * cmath.cos(alpha)
        I2 += fj * cmath.sin(alpha)
    I = np.square(I1) + np.square(I2)
    
    return I
        
def at_scat_factr(ele, dk):
    f = ele 

    return f        
        
class or_setting:
    
        def __init__(self):       
            self.path_save = "no path yet"
            self.method = "no methode selected"            
            self.R_scale = 0
            self.r_scale = 0.0
            self.alpha_max = 0 
            self.alpha_inc = 0 
            self.beta_max = 0 
            self.beta_inc = 0 
            self.gamma_max = 0 
            self.gamma_inc = 0 
            

class lattice:
    def __init__(self):
        self.a = np.zeros(3)
        self.b = np.zeros(3)
        self.c = np.zeros(3)
        self.name = "no name given yet"
        
        def anorm(self):
            return la.norm(self.a)
        
        def bnorm(self):
            return la.norm(self.b)
        
        def cnorm(self):
            return la.norm(self.c)            
            
        def alpha(self):
            return np.degrees(np.arccos(np.dot(self.a,self.b)/(la.norm(self.a)*la.norm(self.b))))

        def beta(self):
            return np.degrees(np.arccos(np.dot(self.a,self.c)/(la.norm(self.a)*la.norm(self.c))))

        def gamma(self):
            return np.degrees(np.arccos(np.dot(self.b,self.c)/(la.norm(self.b)*la.norm(self.c))))


        