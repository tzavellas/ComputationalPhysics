# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:01:02 2017

@author: Anastasios Tzavellas
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x, a, b):
    """
    Line equation, used 
    for curve fitting algorithm
    
    Args:
        x: x value
         a: line coefficient
         b: line constant term
    
    Returns:
        The y coordinate of the point
    """
    return a*x + b

def period(L):
    """
    Calculates the theoretical period of a
    pendulum with length L. It assumes that
    the actual value of g is 9.8m/s2
    
    Args:
        L: The length of the pendulum
    
    Returns:
        The period of the pendulum in sec
    """
    g = 9.8
    return 2*np.pi*np.sqrt(L/g)

def gError(L, T, dL, dT):
    """
    Calculates the error of g 
    estimation given the length, the 
    period and their respective errors
    using the error propagation rule
    
    Args:
        L: A vector of length values
        T: A vector of period values
        dL: The error in length measurement
        dT: The error in period measurement
    
    Returns:
        A vector with g error values
    """
    dg = np.power(2*np.pi, 2) * np.sqrt( np.power(dL/(T*T), 2) + np.power(2*L*dT/np.power(T, 3), 2) )
    return dg

def experiment(L, T, dL, dT, dLsystm = 0):
    """
    Performs a g-measurement experiment
    
    Args:
        L: A vector of length measurements of the pendulum
        T: A vector of period measurements of the pendulum
        dL: The error in length measurement
        dT: The error in period measurement
        dLsystm: Systematic error of length measurement, default value 0
    
    Returns:
        A dictionary with the mean values of g, 
        the g-error values and the measured period 
        values, for each length
    """
    L = L + dLsystm             # Add systematic error, if it exists
    g = np.power(2*np.pi, 2) * L / np.power(T, 2) # Indirect g measurement from
                                                  # length and period
    dg = gError(L, T, dL, dT)   # g measurement error
    gMean = np.sum(g)/g.size    # Mean value of g measurements
    dgMean = np.sqrt(np.sum(dg*dg))/dg.size # Error of mean value of g
    return {'g':gMean, 'dg':dgMean}

def fit(experiment, L, T, dLsystm = 0):
    """
    Performs Least Square Fit on the given experiment
    
    Args:
        experiment: The experiment to perform LSF
        L: A vector of length measurements of the pendulum
        T: A vector of period measurements of the pendulum
        dLsystm: Systematic error of length measurement, default value 0
    
    Returns:
        A dictionary with the LSF value of g, the
        LSF coefficients, and the values used for the fit
    """
    x = np.power(T, 2)
    y = L + dLsystm
    result = curve_fit(line, x, y)  # y = A + Bx
    A = result[0][1]
    B = result[0][0]
    dBA = np.sqrt(np.diag(result[1]))
    g = np.power(2*np.pi, 2) * B    # Coefficient A gives g: A = (2*pi)^2 / g
    dg = np.power(2*np.pi, 2)*dBA[0]# Error of g is using error propagation rule
    return {'g':g, 'dg':dg, 'A':A, 'B':B, 'x':x, 'y':y}


gTheory = 9.8
L = np.array([0.2, 0.4, 0.8, 1.0])
dL = 0.01
dT = 0.1
noise = np.random.standard_normal(L.size) * dT
T = period(L) + noise               # Period measurements with dt=0.1s accuracy

print("Experiment without systematic error")
experiment1 = experiment(L, T, dL, dT)          # Perform experiment 1
print("Mean Value Method")
print("-----------------")
print("g = {:.2f} +- {:.2f}".format(experiment1['g'], experiment1['dg']))
print("\nLeast Squares Fit Method")
print("------------------------")
lsq1 = fit(experiment1, L, T)                   # Perform LSF for experiment 1
print("g = {:.2f} +- {:.2f}".format(lsq1['g'], lsq1['dg']))
xn = lsq1['x']
yn = np.polyval([lsq1['B'], lsq1['A']], xn)
plt.plot(xn, yn, 'r', label='$\delta L_{system} = 0$') # Plot least square line
plt.errorbar(lsq1['x'], lsq1['y'], xerr=dL, yerr=dT, fmt='r.') # Plot measurements
plt.xlabel('$T^2[sec^2]$')
plt.ylabel('$L[m]$')

dLsystm = 0.05
T = period(L+dLsystm) + noise
T = period(L+dLsystm) + np.random.standard_normal(L.size) * dT
print("\n\nExperiment with systematic error={:.2f}".format(dLsystm))
experiment2 = experiment(L, T, dL, dT, dLsystm) # Perform experiment 2,dL = 0.05
print("Mean Value Method")
print("-----------------")
print("g = {:.2f} +- {:.2f}".format(experiment2['g'], experiment2['dg']))
print("\nLeast Squares Fit Method")
print("-------------------------")
lsq2 = fit(experiment2, L, T, dLsystm)
print("g = {:.2f} +- {:.2f}".format(lsq2['g'], lsq2['dg']))
xn = lsq2['x']
yn = np.polyval([lsq2['B'], lsq2['A']], xn)
plt.plot(xn, yn, 'b', label='$\delta L_{system} = 0.05$') # Plot least square line
plt.errorbar(lsq2['x'], lsq2['y'], xerr=dL, yerr=dT, fmt='b*') # Plot measurements
plt.legend()
plt.show()