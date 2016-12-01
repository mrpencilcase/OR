# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 16:32:24 2016

@author: Lukas
"""

import numpy as np
import string
import cmath
import or_class
from numpy import linalg as la
from math import sqrt
import matplotlib.pyplot as plt


x = np.linspace(-100,100, num=1000)
hw = 10
intens = 1000
y = hw /(np.pi*(np.power(x,2)+np.square(hw))) *intens

i6hw = hw /(np.pi*(np.power(6*hw,2)+np.square(hw))) *intens
icent = hw /(np.pi*np.square(hw)) *intens 
print(i6hw)

print(i6hw/icent)

plt.plot(x,y)
plt.show