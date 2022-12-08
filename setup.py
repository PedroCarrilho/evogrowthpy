from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import os
from distutils.sysconfig import get_python_lib

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

pack_name='evogrowthpy'

extensions=Extension(name=pack_name,
                           sources=["evogrowthpy.pyx","Evo.cpp"],
                           include_dirs=["/usr/local/include","../"],
                           libraries=["gsl","gslcblas"],
                           extra_link_args=['-L/usr/local/lib'],
                           language="c++",
                           extra_compile_args=['-std=c++11']
                           )


setup(name=pack_name,
      version="1.0.0",
      author="Pedro Carrilho",
      author_email="pedromgcarrilho@gmail.com",
      cmdclass={'build_ext': build_ext},
      ext_modules = cythonize(extensions,language_level = 3),
      packages=[pack_name],
      package_dir={pack_name: ''}
      )
