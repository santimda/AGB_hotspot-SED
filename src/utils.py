#!/usr/bin/env python
# coding: utf-8

import numpy as np

from constants import *


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

def calc_T_hs(T_star, T_beam, x):
    # The temperature of the hot spot has to be consistent with the assumed stellar temperature 
    # and the average temperature measured in the beam (T_beam). When x=x_max, T_hot = T_beam
    return T_star + (T_beam - T_star)/x**2


