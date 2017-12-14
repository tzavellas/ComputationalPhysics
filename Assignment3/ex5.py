# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import solver


def determinant(A):
    rows = A.shape[0]
    det = 1.  # Det of a triangular matrix is equal
    for i in range(rows):  # to the product of the
        det = det * A[i, i]  # diagonal entries
    return det


A = np.array([[1., 1., 2.], [1., -1., 0.], [-2., 2., 4.]])
B = np.array([7., 1., 2.])
x = np.zeros(3)

result1 = solver.GaussElimination(A, B, 'partial')
det1 = determinant(result1['E'])
print('Partial Pivot Determinant is: ', det1)
print('Partial Pivot Solution is: ', result1['x'])

result2 = solver.GaussElimination(A, B, 'complete')
det2 = determinant(result2['E'])
print('Complete Pivot Determinant is: ', det2)
print('Complete Pivot Solution is: ', result2['x'])
