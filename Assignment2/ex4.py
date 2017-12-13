# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import solver


A = np.array([[1., 1., 2.],
              [-1., 0., 2.],
              [3., 2., -1.]])
b = np.array([1., -3., 8.])
ret = solver.LU(A)
L = ret['L']
U = ret['U']
print('L\n', L)
print('\nU\n', U)
y = solver.ForwardSubstitution(L, b)
print('\ny=', y)
x = solver.BackSubstitution(U, y)
print('x=', x)
