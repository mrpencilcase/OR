import os
import numpy as np

import or_fkt
import or_class

def fit_lattice_planes(MatA, hkl_A, MatB, hkl_B):

    planeA = get_plane(MatA, hkl_A)
    planeB = get_plane(MatB, hkl_B)
    
    