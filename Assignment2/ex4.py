# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import systems


A = np.array([[1., 1., 2.],
              [-1., 0., 2.],
              [3., 2., -1.]])
b = np.array([1., -3., 8.])
ret = systems.LU(A)
L = ret['L']
U = ret['U']
print('L\n', L)
print('\nU\n', U)
y = systems.ForwardSubstitution(L, b)
print('\ny=', y)
x = systems.BackSubstitution(U, y)
print('x=', x)
