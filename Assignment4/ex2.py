# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 22:23:05 2017

@author: tzave
"""

import numpy as np
import matplotlib.pyplot as plt
import diffEquation as de


def f(t, y):
    g = 10.
    A = np.array([[0, 1, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 0, 0]])
    b = np.array([0, 0, 0, -g])
    return np.dot(A, y) + b


v0 = 10.  # Initial conditions
phi = np.pi / 4
v0x = v0 * np.cos(phi)
v0y = v0 * np.sin(phi)
Y0 = np.array([0, v0x, 0, v0y])
h = 0.05

eulerSolution = de.solve(Y0, h, f, de.eulerStep)
rK4Solution = de.solve(Y0, h, f, de.rungeKutta4)
x1 = eulerSolution['x']
y1 = eulerSolution['y']
t1 = eulerSolution['t']
x2 = rK4Solution['x']
y2 = rK4Solution['y']
t2 = rK4Solution['t']

g = 10
x3 = v0x * t2
y3 = v0y * t2 - 1/2 * g * np.power(t2, 2.)

plt.close('all')
plt.figure(1)
plt.plot(t1, x1, '--', label='Euler')
plt.plot(t2, x2, '.', label='RungeKutta4')
plt.plot(t2, x3, label='Analytical')
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(t1, y1, '--', label='Euler')
plt.plot(t2, y2, '.', label='RungeKutta4')
plt.plot(t2, y3, label='Analytical')
plt.grid()
plt.legend()
