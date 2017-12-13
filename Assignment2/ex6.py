# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import systems


#A = np.array([[2., -1., 0., 0.], [-1., 2., -1., 0.], [0., -1., 2., -1.], [0., 0., -1., 2.]])
#d = np.array([1., 0., 0., 1.])
A = np.array([[9., 1, 0], [4, -7, 2], [0, 3, 8]])
d = np.array([5., 6, 2])
taus = np.arange(.1, 2., .1)
for tau in taus:
    result = systems.IterativeMethod(A, d, .5e-6, tau, 1000)
    x = result['x']
    iterations = result['iterations']
    print('tau=', tau, ': ', iterations, ', x=', x)
    print('A*x-d=', np.dot(A, x) - d)

#result = systems.GaussElimination(A, d)
#x = result['x']
#print('GaussElimination: x=', x)
#print('A*x-d=', np.dot(A, x) - d)
