# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 17:51:54 2017

@author: tzave
"""

import numpy as np
import systems


A = np.array([[-1., 2., -1.], [2., -1., 0.], [1., 7., -3.]])
B = np.array([0., 1., 5.])

#result1 = systems.GaussElimination(A, B, 'complete')
#print('Complete Pivot Solution is: ', result1['x'])

result2 = systems.GaussElimination(A, B, 'partial')
print('Partial Pivot Solution is: ', result2['x'])

#result3 = systems.GaussElimination(A, B)
#print('Simple Pivot Solution is: ', result3['x'])
