import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ex

for ent in os.listdir():
    if os.path.isfile(ent) == TRUE:



    
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
