# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 18:57:05 2017

@author: Anastasios Tzavellas
"""

import numpy as np

# ------------------ Function Definitions ------------------

def partialSums(Lmin, Lmax, N):
    """
    Calculates the sum values of function r^2 * (1 + exp(r^2))
    as well as the sum of their squares, at interval [Lmin, Lmax]
    
    Args:
        Lmin: The lower bound of the interval
        Lmax: The upper bound of the interval
        N: The number of random samples
    
    Returns:
        A dictionary with the sum of the
        values and the sum of their squares
    """
    r = np.random.uniform(Lmin, Lmax, N)                # Generate N random samples from U(Lmin, Lmax)
    f = np.power(r, 2) * (1 + np.exp(np.power(r, 2)))   # Evaluate integrant at these points
    meanSum = np.sum(f)                                 # Store their sum in meanSum - Will be used to compute mean value
    sigmaSum = np.sum(f*f)                              # Store the sum of their squares in sigmaSum - Will be used to compute variance
    return {'meanSum':meanSum, 'sigmaSum':sigmaSum}

def StratifiedMonteCarlo(Lmin, Lmax, N, strata):
    """
    Performs Stratified Monte Carlo Integration
    or Crude Monte Carlo integration, if strata=1 
    at interval [Lmin, Lmax]
    
    Args:
        Lmin: The lower bound of the interval
        Lmax: The upper bound of the interval
        N: The number of random samples
        strata: The number of strata
    
    Returns:
        A dictionary with the integral value
        and the error value of the integrals
    """
    r = np.linspace(Lmin, Lmax, strata + 1)             # Split integration interval in strata number of sub-intervals
    Ni = N // strata                                    # N/strata number of samples will be generated within that interval
    integral = np.zeros(strata)                         # Initialize integral and error vectors that store the result of each interval
    error = np.zeros(strata)
    for i in range(strata):
        ret = partialSums(r[i], r[i+1], Ni)             # Calculate meanSum and sigmaSum
        mean = ret['meanSum'] / Ni                      # mean holds the mean value of the samples
        integral[i] =  mean * (r[i+1]-r[i])             # integral[i] holds the integral of the sub-interval
        error[i] =  np.sqrt( ret['sigmaSum']/Ni - np.power(mean, 2)) * (r[i+1]-r[i]) / np.sqrt(Ni) # error[i] holds the error of the sub-interval
                                                        # The total integral is the sum of the integrals of every sub-interval I = I1 + ... + Ik
                                                        # The total error is calculated by the error propagation rule DI = sqrt( DI1^2 + .. + DIk^2)
    return {'integral': np.sum(integral), 'error': np.sqrt(np.sum(error * error))}

# ------------------ Start of the program ------------------
interval = np.array([0, 4]) = 4         # interval of integration
N = 100000                              # Number of Samples

strata = 1
result = StratifiedMonteCarlo(0, L, N, strata)      # Crude Monte Carlo, only 1 stratum
print('\nCrude Monte Carlo ({:d} samples)'.format(N))
print('momentInertia = {:.2e} +- {:4.2e}'.format(result['integral'], result['error']))

strata = 4
result = StratifiedMonteCarlo(0, L, N, strata)      # Stratified Monte Carlo wtih 4 strata
print('\nStratified Monte Carlo ({:d} strata)'.format(strata))
print('momentInertia = {:.2e} +- {:4.2e}'.format(result['integral'], result['error']))

strata = 16
result = StratifiedMonteCarlo(0, L, N, strata)      # Stratified Monte Carlo with 16 strata
print('\nStratified Monte Carlo ({:d} strata)'.format(strata))
print('momentInertia = {:.2e} +- {:4.2e}'.format(result['integral'], result['error']))

print("\nAnalytic Evaluation for Comparison")       # Analytic result for compilation
M = 1.71975e7
print('mass = {:.2e}'.format(M))