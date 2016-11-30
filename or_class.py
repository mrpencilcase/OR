# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:21:24 2016

@author: Lukas
classes used in the OR calculations
"""
import numpy as np 
from numpy import linalg as la
        
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



class unitcell:
    def __init__(self):
        self.name = "no name given"        
        self.atoms =  []
    def ad_atom(self,ele,coordinates):
        adatom = atom()
        adatom.element=ele
        adatom.coord=coordinates
        self.atoms.append(adatom)

class atom:
    def __init__(self):
        self.element = "no element"
        self.coord = [0,0,0]
        