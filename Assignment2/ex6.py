# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 20:02:37 2017

@author: tzave
"""

import numpy as np
import systems


#A = np.array([[2., 3., 0., 0.], [1., 5., 8., 0.], [0., 7., 6., 1.], [0., 0., 9., 7.]])
#d = np.array([1., 2., 3., 4.])
A = np.array([[2., -1., 0., 0.], [-1., 2., -1., 0.], [0., -1., 2., -1.], [0., 0., -1., 2.]])
d = np.array([1., 0., 0., 1.])
#A = np.array([[9., 1, 0], [4, -7, 2], [0, 3, 8]])
#d = np.array([5., 6, 2])
taus = np.arange(.1, 2., .1)
for tau in taus:
    result = systems.IterativeMethod(A, d, .5e-6, tau)
    print('tau=', tau, ': ', result['iterations'], ', x=', result['x'])

result = systems.GaussElimination(A, d)
print('GaussElimination: x=', result['x'])