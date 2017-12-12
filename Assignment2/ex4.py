# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np


def ForwardSubstitution(A, b):
    rows = A.shape[0]
    y = np.zeros(rows)
    for i in range(rows):
        s = 0.
        for j in range(i):
            s = s + A[i, j] * y[j]
        y[i] = (b[i] - s) / A[i, i]
    return y


def BackSubstitution(A, b):
    rows = A.shape[0]
    x = np.zeros(rows)
    for i in reversed(range(rows)):
        s = 0
        for j in range(i + 1, rows):
            s = s + A[i, j] * x[j]
        x[i] = (b[i] - s) / A[i, i]
    return x


def LU(A):
    n, m = A.shape
    assert (n == m)
    U = np.copy(A)
    L = np.identity(n)
    for k in range(n-1):
        for j in range(k+1, n):
            L[j, k] = U[j, k] / U[k, k]
            for p in range(k, n):
                U[j, p] = U[j, p] - L[j, k] * U[k, p]
    return {'L': L, 'U': U}


A = np.array([[1., 1., 2.], [-1., 0., 2.], [3., 2., -1.]])
B = np.array([1., -3., 8.])
ret = LU(A)
L = ret['L']
U = ret['U']
y = ForwardSubstitution(L, B)
x = BackSubstitution(U, y)
