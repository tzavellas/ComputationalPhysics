# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np


def swap(r, p, E):
    n = A.shape
    for q in range(r, n):
        b = E[r][q]
        E[r][q] = E[p][q]
        E[p][q] = b


def BackSubstitution(x, A):
    n = A.shape[0]
    x[n-1] = A[n-1][n]/A[n-1][n-1]
    for i in reversed(range(n-1)):
        sum = 0
        for j in range(i+1, n):
            sum = sum + A[i][j] * x[j]
        x[i] = (A[i][n] - sum) / A[i][i]


def Pivoting(r, E, n):
    PivotFound = False
    for p in range(r, n - 1):
        if np.isclose(E[p][r], 0):  # Keep looking for non zero pivot
            continue
        else:  # if pivot is found
            PivotFound = True
            if p != r:  # Only swap if p>r
                swap(r, p, E)  # swap p and r rows
            break
    return PivotFound


def GaussElimination(x, A, b):
    isSingular = False
    n = b.size
    b = b.reshape(n, 1)
    E = np.append(A, b, axis=1)  # Append b as extra column
    for r in range(n - 1):
        if Pivoting(r, E, n) is False:
            isSingular = True
            print("Matrix is singular")
            break
        for i in range(r + 1, n):
            m = - E[i][r]/E[r][r]
            E[i][r] = 0
            for j in range(r + 1, n + 1):
                E[i][j] = E[i][j] + m * E[r][j]
    if np.isclose(E[n-1][n-1], 0):
        print("There is no unique solution")
    if (isSingular is False):
        BackSubstitution(x, E)
    return {'E': E, 'x': x, 'isSingular': isSingular}


A = np.array([[1, 1, 2], [-1, 0, 2], [3, 2, -1]])
B = np.array([1, -3, 8])
x = np.zeros(3)

g = GaussElimination(x, A, B)
