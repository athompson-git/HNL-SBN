import numpy as np
from numpy import power, sqrt, pi

import sys
sys.path.append("../")
from alplib.constants import *

r_ratio_dat = np.genfromtxt("data/r_ratio_by_sqrtS-GeV.txt")

mu = [2.3, 1275.0, 173210]
md = [4.8, 95.0, 4180.0]

def r_ratio(sqrtS):
    # input sqrtS in GeV
    # from 1801.04847
    return np.interp(sqrtS, r_ratio_dat[:,0], r_ratio_dat[:,1], left=0.0)

def decay_width_zprime_nunu(mZp, gBL=1.0):
    return 3*(gBL**2 * mZp / 24 / pi)  # 3x generations

def decay_width_zprime_dilepton(mZp, ml, gBL=1.0):
    return np.nan_to_num((gBL**2 * mZp / 12 / pi) * sqrt(1 - 4*power(ml/mZp, 2)) * (1 + 2*power(ml/mZp, 2)))

def decay_width_zprime_hadrons(mZp, gBL=1.0):
    width_virtual_gamma_mumu = decay_width_zprime_dilepton(mZp, M_MU, gBL)
    prefactor = np.sum([np.heaviside(mZp - 2*mu[i], 0.0)*decay_width_zprime_dilepton(mZp, mu[i], 1.0) \
                           + np.heaviside(mZp - 2*md[i], 0.0)*decay_width_zprime_dilepton(mZp, md[i], 1.0) for i in range(3)]) \
                / (np.sum([power(2/3, 2)*np.heaviside(mZp - 2*mu[i], 0.0)*decay_width_zprime_dilepton(mZp, mu[i], 1.0) \
                           + power(1/3, 2)*np.heaviside(mZp - 2*md[i], 0.0)*decay_width_zprime_dilepton(mZp, md[i], 1.0) for i in range(3)])) / 9
    threshold = mZp > 770.0
    return threshold*np.nan_to_num(prefactor * width_virtual_gamma_mumu * r_ratio(mZp*1e-3))

def decay_width_zprime_2HNL(mZp, mN, gBL=1.0):
    # x2 for symmetry factor
    return np.nan_to_num((gBL**2 * mZp / 24 / pi) * sqrt(1 - 4*power(mN/mZp, 2)) * (1 - power(mN/mZp, 2)))

def total_width_zprime(mZp, mN, gBL=1.0):
    return decay_width_zprime_nunu(mZp, gBL) + decay_width_zprime_hadrons(mZp, gBL) \
        + decay_width_zprime_dilepton(mZp, M_E, gBL) + decay_width_zprime_dilepton(mZp, M_MU, gBL) \
        + decay_width_zprime_dilepton(mZp, M_TAU, gBL) + decay_width_zprime_2HNL(mZp, mN, gBL)

def br_zprime_2HNL(mZp, mN):
    gamma_total = total_width_zprime(mZp, mN, 1.0)
    return decay_width_zprime_2HNL(mZp, mN, gBL=1.0) / gamma_total
