import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ex


path = os.path.dirname(os.path.realpath(__file__))

for ent in os.listdir(path):
    if os.path.isfile(os.path.join(path,ent)) == TRUE:
        if ".py"  in ent and ".pyx"  not in ent:
            shutil.copyfile(os.path.join(path,ent),os.path.join(path,ent+"x"))            

    
ext_modules=[

            Extension("or_start2", ["or_start.pyx"]),
            Extension("or_main", ["or_main.pyx"]),
            Extension("or_fkt", ["or_fkt.pyx"]),
            Extension("or_class", ["or_class.pyx"]),
            Extension("or_plot_intensity", ["or_plot_intensity.pyx"])
            ]
setup(
    name = "OR_calc",
    cmdclass = {"build_ext" : build_ext},
    ext_modules = ext_modules,
    )
