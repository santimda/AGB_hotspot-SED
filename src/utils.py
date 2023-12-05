#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.integrate import trapz

from constants import *


def nBB(nu, norm, T):
    '''Black-body spectrum: B_nu = 2*h*nu**3/c**2 * 1/(exp(h*nu/k/T)-1)
    Input:
    norm = normalization constant (related to radius)
    nu = frequency in 1e15 Hz (IT MUST BE CONSISTENT WITH THE nsc PARAMETER!!!)
    T = temperature in Kelvin
    '''
    nu2 = nu * nsc  
    return norm * 2. * h * (nu2**3/c**2) / ( np.exp(h*nu2/k_B/T) - 1. )


def get_norm_nBB(R, T):
    nu0 = (k_B * T) / h # Characteristic frequency for BB emission
    nu_min =  nu0 / 300.
    nu_max = nu0 * 1e2
    nu = np.geomspace( nu_min, nu_max, 300 )

    L_norm = nBB(nu/nsc, 1.0, T)
    # Get the normalization just that integral L_E(E) dE = 4 pi R^2 sigma T^4 --> norm = ()/integral 
    integral = trapz(L_norm, nu) 
    norm = (4.0 * np.pi * R**2 * sigma_SB * T**4) / integral
    return norm 


def S_BB(nu, R, T, D):
    '''Returns the SED in mJy'''
    norm = get_norm_nBB(R, T)
    L = nBB(nu, norm, T)
    S = L / (4.0 * np.pi * D**2) / mJy
    return S


def Snu_to_TB(S_nu, nu, Omega_s, Omega_b):
    '''Converts the fluxes S_nu at the frequencies nu to brightness temperature T_B. 
    Omega_s and Omega_b are the solid angles of the source and beam, respectively. 

    Input:
    [S_nu] = mJy 
    [nu] = Hz 
    [Omega_s], [Omega_b] = stereoradian
    
    Output:
    [T_B] = K
    '''
    corr = [ min(1, Omega_s/Omega_i) for Omega_i in Omega_b]      # If the source size is smaller than the beam, the T is "diluted"
    return c**2 * (S_nu*mJy) / (2 * k_B * np.power(nu,2) * Omega_s) * corr 


def Sobs_to_TB(S_nu, nu, Omega_b):
    '''Converts the fluxes S_nu at the frequencies nu to brightness temperature T_B. 
    Omega_b are the solid angles of the beam at each frequency. Similar to Snu_to_TB
    but without knowledge of the true size of the source.

    Input:
    [S_nu] = mJy 
    [nu] = Hz 
    [Omega_b] = stereoradian
    
    Output:
    [T_B] = K
    '''
    return c**2 * (S_nu*mJy) / (2 * k_B * np.power(nu,2) * Omega_b) 


def TB_thick(B, nu, Omega_s, Omega_b):
    '''Returns the brightness temperature for an optically thick synchrotron source 
    using Eq. 5.91 here: https://www.cv.nrao.edu/~sransom/web/Ch5.html.
    Note that this *does not* include the beam factor for unresolved sources!
    
    Input:
    [B] = G 
    [nu] = Hz 
    [Omega_s], [Omega_b] = stereoradian

    Output:
    [T_B] = K
    '''
    corr = [ min(1, Omega_s/Omega_i) for Omega_i in Omega_b]      # If the source size is smaller than the beam, the T is "diluted"
    return 1.18e6 * np.sqrt(nu / B) * corr 


def Omega_beam(theta1, theta2):
    '''Beam solid angle in terms of the theta_HPBW (theta1 x theta2)
    Use the equation here: https://science.nrao.edu/facilities/vla/proposing/TBconv'''
    return pi * theta1 * theta2 / (4 * np.log(2))

