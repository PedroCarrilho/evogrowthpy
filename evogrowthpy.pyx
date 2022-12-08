
import cython
cimport cython
from libc.stdlib cimport free
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import sys as _sys
import numpy as np
from scipy.interpolate import CubicSpline as _CubicSpline
import warnings as _warnings

# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8

#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION



# Import the C++ functions

cdef extern from "Evo.h":
    void get_growth(double A, double omega0, double omegacb, double h, double w0, double wa, double xi, double Om_rc, double * res, double accuracy);


def get_growth_wrap(omega0,hubble,par, accuracy=1e-3,om_cb=0):

    if om_cb:
      omega_cb=om_cb
    else:
      omega_cb=omega0

    z_real=par[0]
    w0=par[1]
    wa=par[2]
    if len(par)>=4:
        xi=par[3]
    else:
        xi=0
    if len(par)==5:
        Om_rc=par[4]
    else:
        Om_rc=0

    res=np.zeros(2)

    cdef double* DF = <double *> PyMem_Malloc(2 * sizeof(double))

    get_growth(1./(1.+z_real),omega0,omega_cb,hubble,w0,wa,xi,Om_rc,DF,accuracy)

    res[0]=DF[0]
    res[1]=DF[1]

    PyMem_Free(DF)

    return res
