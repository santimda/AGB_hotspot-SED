#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.integrate import trapz
from scipy import interpolate

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
    '''Returns the SED in Jy'''
    norm = get_norm_nBB(R, T)
    L = nBB(nu, norm, T)
    S = L / (4.0 * np.pi * D**2) / Jy
    return S


def S_ff(nu, R, D, Te, n, H=None):
    '''Free-free spectrum from an homogeneous cylindrical region with radius and height equal to R.
    It assumes n_i = n_e = n.
    Formulae in Olnon 1975 ( https://articles.adsabs.harvard.edu/pdf/1975A%26A....39..217O )
    
    Input:
    nu = frequency in Hz
    R = source size in cm
    D = source distance in cm
    T_e = electron temperature in Kelvin
    n = gas density in cm^-3
    H = linear depth of the source in cm (default: 0.1*R)
    
    Output:
    S_nu = flux density in Jy
    '''
    
    H = H if H is not None else 0.1 * R

    #Transform units
    nu_GHz = nu / 1e9 # Note that this is usually a vector of size m! 
    m = len(nu)
    
    fs = 8.235e-2 * Te**(-1.35) * nu_GHz**(-2.1) / pc # Eq.4 is in [cm^6 pc^-1]
#     cte_Snu = 2 * np.pi * k_B * Te * R**2 / (c**2 * D**2) # First part in Eq. 11
    ps = [2 * n**2 * H * f for f in fs] # Eq. 10
    Snu = []
    for i in range(m):
        opac = 1. - np.exp(-ps[i]) if ps[i] > 1e-3 else ps[i] 
#         S = cte_Snu * nu[i]**2 * opac * Jy # Eq. 11 Olnon 1975 for sphere (in Jy)
        cte_Snu = nBB(nu[i]/nsc, np.pi*R**2/(D**2), Te)
        S = cte_Snu * opac / Jy # Eq. 11 Olnon 1975 for sphere (in Jy)
        Snu.append(S)
    
    return Snu # in Jy


def CircumstellarDustExtinction(Snu_start,nu,Mdot_dust,vExp,Rin,Rout,kappaFile, dustFolder="../dust_opacities/"):
    ''' Routine to correct an SED for dust extinction.

    Input: 
    Snu_start: unabsorbed SED
    vExp: velocity in km/s (!)

    Output:
    Snu_out: absorption-corrected SED''' 

    Rin = Rin*au
    Rout = Rout*au
    kappaIn = np.loadtxt(dustFolder+kappaFile,usecols=(0,1)).transpose()
    kappaIn[0] = (299792.458/kappaIn[0])*1E9    
    f = interpolate.interp1d(kappaIn[0], kappaIn[1])
    kappa = f(nu)
    
    Mdot_dust = Mdot_dust * M_sun_yr
    #vExp = vExp * 1E3  # Convert from km/s to m/s
    vExp = vExp * 1E5  # Convert from km/s to cm/s
    
    N_dust = Mdot_dust/(2*np.pi*vExp)*((1.0/Rin) - (1.0/Rout)) #in g/cm^2
    tau  = kappa * N_dust
    Snu_out = Snu_start * np.exp(-tau) 
    
    return Snu_out


def S_ff_abs(nu, R, D, Te, n, H, Mdot_dust, vExp, Rin, Rout, kappaFile):
    '''Free-free spectrum corrected for dust absorption'''
    Snu_ff = S_ff(nu, R, D, Te, n, H)
    Snu_abs = CircumstellarDustExtinction(Snu_ff, nu, Mdot_dust, vExp, Rin, Rout, kappaFile) 
    return Snu_abs