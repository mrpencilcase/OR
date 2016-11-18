# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 11:38:18 2016

@author: Lukas

Functions used to find the most probable orientational relashionships of an 
hetrocrystal interface. The underlying methode is the Reciprocl Lattice Point
methote (RLP) by Ikuhara and Priouz.

"""
import numpy as np

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
def reziprocal_lattice(lattice):    
    
    rl = np.linalg.pinv(lattice)
    
    return rl

# Find the last reziprocal points miller indizes that lies within the the given radius R
# along the reziprokal vector rezVec.
def hkl_max(rezVec,R):
    miller = 1
    g =np.sqrt(np.square(rezVec[0]*miller)+np.square(rezVec[1]*miller)+np.square(rezVec[2]*miller))
    while g <= R:
        miller = miller + 1
        g =np.sqrt(np.square(rezVec[0]*miller)+np.square(rezVec[1]*miller)+np.square(rezVec[2]*miller))
    return miller -1