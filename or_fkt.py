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
import or_class
from numpy import linalg as la
from math import sqrt

def orthon_trans(a_v,b_v,c_v,alpha,beta,gamma):
    """
    Returns the transformation matrix from the inital crystal basis to 
    a orthogonal one. 

    Input:
    a,b and c: lengths of the unit cell vectors
    alpha, beta and gamma: angles of the unit cell
    
    Output:
    T = |a1|^T  where a1, a2 and a3 are the direct lattice vectors of 
        |a2|    the orhtonormal system. 
        |a3| 
    """    
    a = la.norm(a_v)
    b = la.norm(b_v)
    c = la.norm(c_v)

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
    
def new_at_coord(T, unti_cell_base):
    """
    transforms array of atom positions to a new basis
    uing the transformation matrix T
    
    Input:
    T: transformation matrix as a numpy array
    unti_cell_base: the atom coordinates in form of the custom variable unitcell 

    Output:
    at_coord_new: the atom coordinates after the transformatin in form of the custom variable unitcell
    """
    unit_cell_trans = or_class.unitcell()
    unit_cell_trans.name = unit_cell_base.name

    for atom in unit_cell_case.atoms:
        at_coord_new = 1

   
    return at_coord_new
# Calculate the reziprokal lattice vectors a*, b* and c* by inverting orhonomral basis 
def reziprocal_lattice_Gautam(lattice):    
    
    rl = np.linalg.pinv(lattice)
    
    return rl
# Conventional algorithm to calculate the reciprocal lattice
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
    
     
# Rotate lattice around the [1,0,0] axis
def rot_x(lattice,alpha):
    c = np.cos(alpha)
    s = np.sin(alpha)
    
    rot_m = np.array([ [ c, s, 0] ,
                       [-s, c, 0] ,
                       [  0,0, 1] ])
    
    lattice_rot = np.dot(lattice,rot_m)
    return lattice_rot    
# Rotate lattice around the [0,0,1] axis
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
    CONTCAR=open(path,'r')
    name = CONTCAR.readline()
    aSc=float(CONTCAR.readline())
    a1In=np.asarray([aSc*float(x) for x in CONTCAR.readline().split()])
    a2In=np.asarray([aSc*float(x) for x in CONTCAR.readline().split()])
    a3In=np.asarray([aSc*float(x) for x in CONTCAR.readline().split()])
    atom_type = [x for x in CONTCAR.readline().split()]
    atom_numb = [int(x) for x in CONTCAR.readline().split()]
    CONTCAR.readline()    
    
    i = 1 
    end = 0 
    unit_cell = or_class.unitcell()    
    unit_cell.name = name.rstrip()
    # Get the absolute position of each individual atom in the unit cell. 
    for x,y in zip(atom_type,atom_numb):
        end = end + y        
        while i <= end:
            a_coord_rel = [float(z) for z in CONTCAR.readline().split()]
            a_corrd_abs = a1In *a_coord_rel[0] + a2In *a_coord_rel[1] + a2In *a_coord_rel[2]
            unit_cell.ad_atom(x,a_corrd_abs)
            i = i+1
    CONTCAR.close()
    # process lattice parameters
    data = or_class.lattice()
    data.a = a1In
    data.b = a2In
    data.c = a3In
    data.anorm = la.norm(a1In)/supercell[0]
    data.bnorm = la.norm(a2In)/supercell[1]
    data.cnorm = la.norm(a3In)/supercell[2]
    data.alpha = (np.degrees(np.arccos(np.dot(a1In,a2In)/(la.norm(a1In)*la.norm(a2In)))))
    data.beta = (np.degrees(np.arccos(np.dot(a1In,a3In)/(la.norm(a1In)*la.norm(a3In)))))
    data.gamma = (np.degrees(np.arccos(np.dot(a2In,a3In)/(la.norm(a2In)*la.norm(a3In)))))

    return data, unit_cell
    
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
    latt = or_class.lattice()
    latt.a=a
    latt.b=b
    latt.c=c
    return latt
    
def intensity_lattice_point(unit_cell,g):
    
    g_norm = la.norm(g)   
    I1 = 0
    I2 = 0  
    dummy1 = 0
    dummy2 = 0 
    pi2 = 2*np.pi      

    for ent in unit_cell.atoms:
        fj = at_scat_factr(ent.element,g_norm)
        rj = ent.coord
        alpha = np.dot(rj,g)*pi2
        I1 += fj * cmath.cosh(alpha)
        I2 += fj * cmath.sinh(alpha)
        dummy1 += fj
        dummy2 += fj

    I = np.square(I1) - np.square(I2)
    return I

"""
Function to find the atomic scattering factore of an element(ele) in dependence of 
the reciprocal vector g. The scattering factors are taken from the book
Transmission Electron Microscopy and Diffractometry of Materials by B.Fultz and
H.M.Howe, Springer and can be found in the file atomic_scattering_factor.txt
"""     
def at_scat_factr(ele, dk):
    with open("atomic_scattering_factor.txt","r") as dat:
        scat_factr = list(dat)
        s_inc = [float(x) for x in scat_factr[0].split()]
        for ent in scat_factr: 
            line = [x for x in ent.split()]
            if line[0]  == ele:
                line.pop(0)
                line = [float(x) for x in line]
                s = dk/(np.pi*4)
                i = 1 
                while i < len(line):
                    if s < s_inc[i]:
                        break
                    i += 1
                
                f = (line[i]+line[i-1])/2
                break
    return f


def intensity_overlap_part(I1,Ri,Rj,d,HW):
    """
    Calculate the intensity in a volume of two partially overlapping spheres.
     
    """
    r  = Ri
    ii = (Rj*Rj-d*d)/(4*np.pi*d*HW)*np.log((r*r)/(r*r+HW*HW)) - HW/(4*d*np.pi)*np.log(r*r+HW*HW)+np.arctan(r/HW)/np.pi
    r  = d-Rj
    ii = ii -(Rj*Rj-d*d)/(4*np.pi*d*HW)*np.log((r*r)/(r*r+HW*HW)) - HW/(4*d*np.pi)*np.log(r*r+HW*HW)+np.arctan(r/HW)/np.pi
    ii = ii * I1
      
    return ii

def intensity_overlap_tot(Ij,Ri,d,HWj):
    """
    Calculate the intensity of a sphere that has the ohter sphere fully within itself. This is important 
    when one sphere is significantly larger than the other.
    """
    HWj2 = HWj*HWj
    d2 = d*d
    a = Ri + d
    a2 = a*a
    ij =   d * np.log((a2 + (HWj2)) * d2 /(a2 * (d2 + HWj2))) + Ri*d/a
    ij = (ij + (HWj2 - d2) /HWj * (np.arctan((Ri + d)/ HWj)- np.arctan(d(HWj))))* 2*Ij /(np.pi*HWj)

    return ij