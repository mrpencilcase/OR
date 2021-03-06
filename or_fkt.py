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
    
def new_at_coord(T, unit_cell_base):
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

    for atom in unit_cell_base.atoms:
        at_coord_new = 1
    return at_coord_new

# Calculate the reziprokal lattice vectors a*, b* and c* by inverting orhonomral basis 
def reziprocal_lattice_Gautam(lattice):    
    """
    Calculates the reciprocal lattice of Gautams orthonormal lattice
    by inverting the transformation matrix T. 

    Input: 
    lattice = Transformation matrix T

    Output: 
    rl  = reciprocal lattice
    """
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
    """
    Opens a file that contains information about the materials unit cell
    at the given path and stores the lattice parameter as well as the atom
    positions in seperate variable which are then returned.
    
    Input: 
    path = path to the file including filename and  file extension
    supercell = Information if the data in the file represents a singel unit cell
                or multiple cells. The variable is expectet as a list with three entries
                representing the number of unit cell in x,y and z direction [nx,ny,nz].
                
    Output:
    lattice = Variable of the class Lattice which contains the geometrical
              data of the unit cell.  
    unit_cell = Variable of the class UnitCell wich contains the positons of the
                singel atoms within the unit cell.
    """

    # Read the lattice data

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
    unit_cell = or_class.UnitCell()    
    unit_cell.name = name.rstrip()
    unit_cell.elements = atom_type
    # Get the absolute position of each individual atom in the unit cell. 
    for x,y in zip(atom_type,atom_numb):
        end = end + y        
        while i <= end:
            a_coord_rel = [float(z) for z in CONTCAR.readline().split()]
            a_corrd_abs = a1In *a_coord_rel[0] + a2In *a_coord_rel[1] + a3In *a_coord_rel[2]
            #unit_cell.ad_atom(x,a_corrd_abs)
            unit_cell.ad_atom(x,a_coord_rel)
            i = i+1
    CONTCAR.close()
    # process lattice parameters
    lattice = or_class.Lattice()
    lattice.a = a1In
    lattice.b = a2In
    lattice.c = a3In
    lattice.anorm = la.norm(a1In)/supercell[0]
    lattice.bnorm = la.norm(a2In)/supercell[1]
    lattice.cnorm = la.norm(a3In)/supercell[2]
    lattice.alpha = (np.degrees(np.arccos(np.dot(a1In,a2In)/(la.norm(a1In)*la.norm(a2In)))))
    lattice.beta = (np.degrees(np.arccos(np.dot(a1In,a3In)/(la.norm(a1In)*la.norm(a3In)))))
    lattice.gamma = (np.degrees(np.arccos(np.dot(a2In,a3In)/(la.norm(a2In)*la.norm(a3In)))))

    return lattice, unit_cell
    
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
    latt = or_class.Lattice()
    latt.a=a
    latt.b=b
    latt.c=c
    return latt
    
def intensity_lattice_point1(unit_cell,g):
    
    g_norm = la.norm(g)   
    I1 = 0
    I2 = 0  
    pi2 = 2*np.pi      
    i = 1
    for ent in unit_cell.elements:
        fj = at_scat_factr(ent,g_norm)
        for at in unit_cell.atoms:
            if ent == at.element:
                print(i)
                i += 1
                rj = at.coord
                alpha = np.dot(rj,g)
                alpha = alpha*pi2
                I1 += fj * cmath.cos(alpha)
                I2 += fj * cmath.sin(alpha)

    I = (np.square(I1) - np.square(I2)).real
    return I

def intensity_lattice_point2(unit_cell,hkl,g,T):
    
    g_norm = la.norm(g)   
    I1 = 0
    I2 = 0  
    pi2 = 2*np.pi      
    i = 1
    for ent in unit_cell.elements:  
        fj = at_scat_factr(ent,g_norm)      
        for at in unit_cell.atoms:
            if ent == at.element:    
                rj = np.dot(at.coord,T)
                alpha = rj[0]*hkl[0] + rj[1]*hkl[1]+ rj[2]*hkl[2]
                alpha = alpha*pi2
                I1 = I1 + fj * cmath.cos(alpha)
                I2 = I2 + fj * cmath.sin(alpha)

    I = np.square(I1) - np.square(I2)
    return I.real

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
        s = dk/(np.pi*4)
        for ent in scat_factr: 
            line = [x for x in ent.split()]
            if line[0]  == ele:
                line.pop(0)
                line = [float(x) for x in line]             
                i = 1 
                while i < len(line):
                    if s < s_inc[i]:
                        break
                    i += 1               
                f = (line[i]+line[i-1])/2
                break
    return f


def intensity_overlap_part(Ii,Ri,Rj,d,HWi):
    """
    Calculate the intensity in a volume of two partially overlapping spheres.
     
    """
    r  = Ri
    ii =     (Rj*Rj-d*d)/(4*np.pi*d*HWi)*np.log((r*r)/(r*r+HWi*HWi)) - HWi/(4*d*np.pi)*np.log(r*r+HWi*HWi)+np.arctan(r/HWi)/np.pi
    r  = d-Rj
    ii = ii -(Rj*Rj-d*d)/(4*np.pi*d*HWi)*np.log((r*r)/(r*r+HWi*HWi)) - HWi/(4*d*np.pi)*np.log(r*r+HWi*HWi)+np.arctan(r/HWi)/np.pi
    ii = ii * Ii
      
    return ii

def intensity_overlap_tot(Ij,Ri,d,HWj):
    """
    Calculate the intensity of a sphere that has the ohter sphere fully within itself. This is important 
    when one sphere is significantly larger than the other.
    Inpunt:
    Ij = Intensity of bigger sphere
    Ri = Radius of smaller sphere
    d = Distance between the spheres
    HWj = Half width at half maximum of Intensity in j 
    """
    HWj2 = HWj*HWj
    d2 = d*d
    a = Ri + d
    a2 = a*a
    ij =   d * np.log((a2 + (HWj2)) * d2 /(a2 * (d2 + HWj2))) + Ri*d/a
    ij = (ij + (HWj2 - d2) /HWj * (np.arctan((Ri + d)/ HWj)- np.arctan(d/HWj)))* 2*Ij /(np.pi*HWj)

    return ij

def find_max_inten(path, count):
  
    raw_data = []
    return_data = []

    with open(path, "r") as dat:
        data = list(dat)

        for ent in data:
            if "#" not in ent:   
                line = [x for x in ent.split()]          
                intens = (float(line[0]))
                alpha = (float(line[1]))
                beta = (float(line[2]))
                gamma = (float(line[3]))
                raw_data.append([intens,alpha,beta,gamma])

    data_sorted = sorted(raw_data, key=getKey)

    for i in range(count):
        return_data.append(data_sorted[-(i+1)])
    
    return return_data

def getKey(item):
    return item[0]

def get_hkl(latticeA,latticeB,hkl):

    hkl_a = []
    hkl_b = []

    if la.norm(latticeA[0,:]) <= la.norm(latticeB[0,:]):
       
        hkl_a.append(int(hkl[0] * la.norm(latticeB[0,:]) / la.norm(latticeA[0,:])))
        hkl_a.append(int(hkl[1] * la.norm(latticeB[1,:]) / la.norm(latticeA[1,:])))
        hkl_a.append(int(hkl[2] * la.norm(latticeB[2,:]) / la.norm(latticeA[2,:])))
       
        hkl_b.append(hkl[0])
        hkl_b.append(hkl[1])
        hkl_b.append(hkl[2])
       
         
    elif la.norm(latticeB[0,:]) < la.norm(latticeA[0,:]):

       hkl_b.append(int(hkl[0] * la.norm(latticeA[0,:]) / la.norm(latticeB[0,:])))
       hkl_b.append(int(hkl[1] * la.norm(latticeA[1,:]) / la.norm(latticeB[1,:])))
       hkl_b.append(int(hkl[2] * la.norm(latticeA[2,:]) / la.norm(latticeB[2,:])))

       hkl_a.append(hkl[0])
       hkl_a.append(hkl[1])
       hkl_a.append(hkl[2])

    return hkl_a, hkl_b

def calc_intensities(hkl, lattice,unit_cell,T):
    g = [] 
    intens =  []
    hkl_v = []

    h = -hkl[0]

    while h <= hkl[0]:
        gh = lattice[0] * h       
        k = -hkl[1]    
        while k <= hkl[1]:
            ghk = lattice[1] * k + gh
            l = -hkl[2]
            while l <= hkl[2]:
                if sum(map(abs,[h,k,l])) > 0: 
                    ghkl = lattice[2] * l + ghk
                    intens.append(intensity_lattice_point2(unit_cell,[h,k,l],ghkl,T))
                    g.append(ghkl)   
                    hkl_v.append([h,k,l])     
                l += 1                
            k += 1
        h += 1

    return intens, g, hkl_v

def overlap_lattices(intensA,gA,hkl_a,intensB,gB,hkl_b,ImaxA,ImaxB,delta0):
    tot_intens = 0  
    singel_intens = [] 
    numb_overl = 0 

    for intAi, gAi, hkla in zip(intensA,gA,hkl_a):
        for intBj, gBj, hklb in zip(intensB,gB,hkl_b):
                                                   
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
                    I1 = intensity_overlap_tot(intAi,Rj,d,HWi)
                    I2 = intBj
                    singel_intens.append([I1+I2, gAi,hkla,gBj,hklb])
                    tot_intens += I1 + I2
                elif d + Ri <= Rj:
                    I1 = intensity_overlap_tot(intBj,Ri,d,HWj)
                    I2 = intAi
                    singel_intens.append([I1+I2, gAi,hkla,gBj,hklb])
                    tot_intens += I1 + I2
                else:
                    I1 = intensity_overlap_part(intAi,Ri,Rj,d,HWi)
                    I2 = intensity_overlap_part(intBj,Rj,Ri,d,HWj)
                    singel_intens.append([I1+I2, gAi, hkla , gBj, hklb])
                    tot_intens += I1 + I2
    return tot_intens , singel_intens

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def get_planes(Mat,hkl):
    
    
    plane_atoms = []
    R = 2 # radius around the atom positions

    xs = 1
    ys = 1
    zs = 1
    A = 1
    B = 1 
    C = 1 
    
    for ent in Mat.unitcell.atoms:
    
        t = A * xs + B * ys + C * zs
        t = t/(A*A + B*B + C*C)
        xc = xs + A * t
        yc = ys + B * t
        zc = zs + C * t

        d = np.sqrt(np.square(xs-xc) + np.square(ys-yc) +np.square(zs-zc))

        if d <= R:
            plane_atoms.append([xc,yc,zc])
