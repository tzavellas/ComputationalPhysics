# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import solver


A = np.array([[2., -1., 0., 0.],
              [-1., 2., -1., 0.],
              [0., -1., 2., -1.],
              [0., 0., -1., 2.]])
d = np.array([1., 0., 0., 1.])

taus = np.arange(.1, 2., .1)
for tau in taus:
    result = solver.IterativeMethod(A, d, .5e-6, tau, 1000)
    x = result['x']
    iterations = result['iterations']
    print('tau=', tau, ': ', iterations, ', x=', x)
    delta = np.dot(A, x) - d
    print('Norm (Ax-d) =', np.sqrt(np.dot(delta, delta)))
