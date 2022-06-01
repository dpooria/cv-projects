#! /usr/bin/env python3

import pathlib
import ctypes
import numpy as np
from sklearn.datasets import load_boston

C_p = ctypes.c_void_p
C_i = ctypes.c_int
C_d = ctypes.c_double
X = np.load("Hoda_data.npz")["img"]
Y = np.load("Hoda_data.npz")["target"]
n = X.shape[0]
m = X.shape[1] * X.shape[2]
c = Y.shape[1]
x = np.array(X, dtype=C_i)
y = np.array(Y, dtype=C_i)

libpath = pathlib.Path().absolute() / "libread_data.so"
lib = ctypes.CDLL(libpath)
lib.read_hoda.argtypes = [C_p, C_p, C_i, C_i, C_i]
lib.read_hoda(C_p(x.ctypes.data),
              C_p(y.ctypes.data), C_i(n), C_i(m), C_i(c))

X_boston, y_boston = load_boston(return_X_y=True)
X_boston = np.array(X_boston, dtype=C_d)
y_boston = np.array(y_boston, dtype=C_d)
n_boston = len(y_boston)
m_boston = X_boston.shape[1]
lib.read_boston.argtypes = [C_p, C_p, C_i, C_i]
lib.read_boston(C_p(X_boston.ctypes.data), 
        C_p(y_boston.ctypes.data), C_i(n_boston), C_i(m_boston))
print("data converted successfully")
