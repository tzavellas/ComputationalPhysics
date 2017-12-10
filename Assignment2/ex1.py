# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 20:28:24 2017

@author: tzave
"""

import numpy as np
import matplotlib.pyplot as plt


def fixed_point(x0, g, tol, N = None):
	i = 1
	p0 = x0
	while (True if N is None else (i<N) ):
		p = g(p0)
		if np.abs(p-p0)< tol :
			break
		i = i + 1
		p0 = p
	return {'p': p, 'iterations': i}


def f(x):
	return np.power(x, 2) + 10 * np.cos(x)


def g1(x):
	return np.sqrt(-10 * np.cos(x))


def dg1(x):
	return 10*np.sin(x) / (2 * np.sqrt(-10*np.cos(x)))


def g2(x):
	return np.arccos( - np.power(x, 2) / 10 )


def dg2(x):
	return 2*x/np.sqrt(100 - np.power(x,4))


x0 = 2.6
roots = np.zeros(4)
iterations = np.zeros(2)

result = fixed_point(x0, g1, 1e-4)
roots[0] = result['p']
roots[1] = - roots[0]
iterations[0] = result['iterations']
print('Root found after ', iterations[0], ' iterations')

result = fixed_point(x0, g2, 1e-4)
roots[2] = result['p']
roots[3] = - roots[2]
iterations[1] = result['iterations']
print('Root found after ', iterations[1], ' iterations')

x = np.arange(-4.5, 4.5, 0.1)
plt.close('all')
plt.figure(1)
plt.plot(x, f(x), roots[0], f(roots[0]), '*r', roots[1], f(roots[1]), '*r', roots[2], f(roots[2]), '*r', roots[3], f(roots[3]), '*r')
plt.grid()

x = np.arange(1.95, 3.17, 0.01)
plt.figure(2)
plt.plot(x, dg1(x), label = 'dg1')
plt.plot(x, dg2(x), label = 'dg2')
plt.plot(x, np.ones(x.shape), label = 'y = 1')
plt.plot(x, -np.ones(x.shape), label = 'y = -1')
plt.ylim(-1.5, 1.5)
plt.legend()
plt.grid()