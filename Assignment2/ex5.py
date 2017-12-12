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


def BackSubstitution(x, A):
    rows = A.shape[0]
    x[rows-1] = A[rows-1, rows]/A[rows-1, rows-1]
    for i in reversed(range(rows-1)):
        sum = 0
        for j in range(i+1, rows):
            sum = sum + A[i, j] * x[j]
        x[i] = (A[i, rows] - sum) / A[i, i]


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


def GaussElimination(x, A, b, pivot=None):
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
        BackSubstitution(x, E)
    return {'E': E, 'x': x, 'isSingular': isSingular}


def determinant(A):
    rows = A.shape[0]
    det = 1
    for i in range(rows):
        det = det * A[i, i]
    return det


#A = np.array([[1, 1, 2], [-1, 0, 2], [3, 2, -1]])
#B = np.array([1, -3, 8])
A = np.array([[1., 1., 2.], [1., -1., 0.], [-2., 2., 4.]])
B = np.array([7., 1., 2.])
x = np.zeros(3)

result1 = GaussElimination(x, A, B, 'complete')
det1 = determinant(result1['E'])
print('Complete Pivot Determinant is: ', det1)
print('Complete Pivot Solution is: ', result1['x'])

result2 = GaussElimination(x, A, B, 'partial')
det2 = determinant(result2['E'])
print('Partial Pivot Determinant is: ', det2)
print('Partial Pivot Solution is: ', result2['x'])

result3 = GaussElimination(x, A, B)
det3 = determinant(result3['E'])
print('Simple Pivot Determinant is: ', det3)
print('Simple Pivot Solution is: ', result3['x'])



A = np.array([[-1., 2., -1.], [2., -1., 0.], [1., 7., -3.]])
B = np.array([0., 1., 5.])
res = GaussElimination(x, A, B, 'partial')
