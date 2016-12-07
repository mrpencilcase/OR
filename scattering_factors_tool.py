# -*- coding: utf-8 -*-
"""
Created on Fri Dec 07 15:55:35 2016

@author: Lukas
short script to format the text from the ocr programm which was used to convert the scans of the 
atomic scattering factors.

"""
import numpy as np
#path of the textfile
def read_sf():  
    path = "C:\Users\Lukas\Documents\Diplomarbeit\Python_Scripts\OR/test.txt"
    rows = 36
    columms = 25
    #scat_factr = np.zeros((lines,rows))
    scat_factr = []
    with open(path,"r") as raw_text:
        sf = list(raw_text)
    i = 0
    j = 0
    line = []
    for ent in sf:
        #scat_factr[i,j] = ent.rstrip()
        line.append(str(ent.rstrip()))
        if i == rows-1:
            j = j+1
            i = 0
            scat_factr.append(line)
            line = []
        i = i +1 

    scat_factr_tr = map(list,zip(*scat_factr))
    with open("C:\Users\Lukas\Documents\Diplomarbeit\Python_Scripts\OR/testout.txt","w") as out:
        for row in scat_factr:
            for ent in row: 
               out.write("{}".format(ent))
            out.write("\n")
