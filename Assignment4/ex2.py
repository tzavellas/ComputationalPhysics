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


def run(Y0, f, method, h, figNum):
    t = np.array([])
    for hi in h:
        solution = de.solve(Y0, f, method, hi)
        print("h=%.3f, iterations=%d"% (hi,
                                        solution['iterations']))
        x = solution['x']
        y = solution['y']
        t = solution['t']
        plt.figure(figNum)
        plt.plot(t, x, '--', label="h=%.3f" % hi)
        plt.figure(figNum + 1)
        plt.plot(t, y, '--', label="h=%.3f" % hi)
    g = 10.
    v0x = Y0[1]
    v0y = Y0[3]
    x = v0x * t
    y = v0y * t - 1/2 * g * np.power(t, 2.)
    plt.figure(figNum)
    plt.plot(t, x, label="Analytical")
    plt.grid()
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.figure(figNum+1)
    plt.plot(t, y, label="Analytical")
    plt.grid()
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('y(t)')


v0 = 10.  # Initial conditions
phi = np.pi / 4
v0x = v0 * np.cos(phi)
v0y = v0 * np.sin(phi)
Y0 = np.array([0, v0x, 0, v0y])  # initial vector
h = np.array([0.005, 0.01, 0.05, 0.1])  # time steps
t = np.array([])
plt.close('all')

print('Euler')
run(Y0, f, de.eulerStep, h, 1)  # run with Euler
print('\nRunge-Kutta 4th')
run(Y0, f, de.rungeKutta4, h, 3)  # run with Runge-Kutta 4th
