# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import systems


def determinant(A):
    rows = A.shape[0]
    det = 1.
    for i in range(rows):
        det = det * A[i, i]
    return det


A = np.array([[1., 1., 2.], [1., -1., 0.], [-2., 2., 4.]])
B = np.array([7., 1., 2.])
x = np.zeros(3)

result1 = systems.GaussElimination(A, B, 'complete')
det1 = determinant(result1['E'])
print('Complete Pivot Determinant is: ', det1)
print('Complete Pivot Solution is: ', result1['x'])

result2 = systems.GaussElimination(A, B, 'partial')
det2 = determinant(result2['E'])
print('Partial Pivot Determinant is: ', det2)
print('Partial Pivot Solution is: ', result2['x'])

result3 = systems.GaussElimination(A, B)
det3 = determinant(result3['E'])
print('Simple Pivot Determinant is: ', det3)
print('Simple Pivot Solution is: ', result3['x'])
