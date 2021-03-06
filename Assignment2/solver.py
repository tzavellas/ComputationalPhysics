# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np


def swap(r, p, E, cols=None):
    if cols is None:
        b = np.copy(E[r, :])
        E[r, :] = E[p, :]
        E[p, :] = np.copy(b)
    else:
        b = np.copy(E[:, r])
        E[:, r] = E[:, p]
        E[:, p] = np.copy(b)


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


def foundNonZeroPivot(r, E):
    rows = E.shape[0]
    PivotFound = False
    for p in range(r, rows - 1):
        if np.isclose(E[p, r], 0):  # Keep looking for non-zero pivot
            continue
        else:  # if pivot is found
            PivotFound = True
            if p > r:  # Only swap if p>r
                swap(r, p, E)
            break
    return PivotFound


def partialPivot(r, A):
    rows = A.shape[0]  # Number of rows
    Amax = np.abs(A[r, r])
    rmax = r
    for p in range(r, rows):
        Apr = np.abs(A[p, r])
        if Apr > Amax:
            Amax = Apr
            rmax = p
    if rmax != r:
        swap(r, rmax, A)
    return


def completePivot(r, A):
    rows, cols = A.shape
    cols = cols - 1  # ignore the last column of the augmented matrix
    Amax = np.abs(A[r, r])
    rmax = r
    cmax = r
    for i in range(r, rows):
        for j in range(r, cols):
            Aij = np.abs(A[i, j])
            if Aij > Amax:
                Amax = Aij
                rmax = i
                cmax = j
    if (rmax != r) and (cmax != r):
        swap(r, rmax, A)
        swap(r, cmax, A, True)
    return


def GaussElimination(A, b, pivot=None):
    isSingular = False
    rows = A.shape[0]
    b = b.reshape(rows, 1)
    E = np.append(A, b, axis=1)  # Append b as extra column
    for r in range(rows - 1):
        if pivot == 'partial':
            partialPivot(r, E)
        elif pivot == 'complete':
            completePivot(r, E)
        else:  # Simple pivot is required to avoid division by 0
            isSingular = not foundNonZeroPivot(r, E)
            if isSingular:
                break
        for i in range(r + 1, rows):
            if np.isclose(E[i, r], 0):  # skip line if pivot is already 0
                continue
            m = - E[i, r]/E[r, r]
            E[i][r] = 0
            for j in range(r + 1, rows + 1):
                E[i, j] = E[i, j] + m * E[r, j]
    if isSingular:
        print("Matrix is singular")
    elif np.isclose(E[rows-1, rows-1], 0):
        print("There is no unique solution")
    else:
        y = E[:, rows]
        E = np.delete(E, [rows], axis=1)
        x = BackSubstitution(E, y)
    return {'E': E, 'x': x, 'isSingular': isSingular}


def LU(A):
    rows, cols = A.shape
    assert (rows == cols)
    U = np.copy(A)
    L = np.identity(rows)
    for k in range(rows-1):
        for j in range(k+1, rows):
            L[j, k] = U[j, k] / U[k, k]
            for p in range(k, rows):
                U[j, p] = U[j, p] - L[j, k] * U[k, p]
    return {'L': L, 'U': U}


def inverseU(U):
    rows = U.shape[0]
    Uinv = np.zeros(U.shape)
    unit = np.identity(rows)
    for j in range(rows):
        b = unit[:, j]
        Uinv[:, j] = BackSubstitution(U, b)
    return Uinv


def IterativeMethod(A, d, tol, tau=0.1, N=None):
    rows = A.shape[0]
    Dinv = np.diag(1./A.diagonal())
    CL = -np.tril(A, -1)
    CU = -np.triu(A, 1)

    L = np.dot(Dinv, CL)
    U = np.dot(Dinv, CU)

    xk = d
    M = np.identity(rows) - U
    Minv = inverseU(M)
    C = np.dot(Dinv, d)
    C = tau * np.dot(Minv, C)
    T = (1.-tau) * np.identity(rows) + tau * np.dot(Minv, L)
    i = 1
    while (True if N is None else (i < N)):
        x = np.dot(T, xk) + C
        dx = x - xk  # Check if scheme converges
        if np.sqrt(np.dot(dx, dx)) < tol:
            break
        i = i + 1
        xk = x  # Prepare for next iteration
    return {'x': x, 'dx': dx, 'iterations': i}
