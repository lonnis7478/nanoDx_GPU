import os
from os.path import join as pjoin
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext = Extension('cudaclassifier',
                sources=['cuda_wrapper.pyx'],
                extra_objects=['kernel'],
                libraries=['cudart', 'curand', 'cuda'],
                language='c++',
                include_dirs=["/usr/local/cuda-12/include" , "/home/sander/Documents/nanoDxGPU/nanoDx_GPU/scripts",  np.get_include()],
                library_dirs=["/usr/local/cuda-12/lib64"],
                extra_compile_args=['-fopenmp'],
                extra_link_args=[]
                )

setup(
    name='calibratedClassifer',
    ext_modules=[ext],
    cmdclass={'build_ext': build_ext},
)