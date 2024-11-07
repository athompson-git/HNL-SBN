import sys
sys.path.append("../")

from alplib.constants import *
from alplib.fmath import *


theta_12 = 33.45*2*pi/360 #rad
theta_13 = 8.60*2*pi/360; #rad
theta_23 = 52.1*2*pi/360; #rad
   
delta_cp = -3*pi/2 # CP violation angle


#PMNS Matrix---------------------------------------------------------------
U_23 = np.array([[1, 0, 0],[0, np.cos(theta_23), np.sin(theta_23)],[ 0, -np.sin(theta_23),  np.cos(theta_23)]])
U_13 = np.array([[np.cos(theta_13), 0, np.sin(theta_13)*np.exp(-1j*delta_cp)], [0, 1, 0], [-np.sin(theta_13)*np.exp(1j*delta_cp), 0, np.cos(theta_13)]])
U_12 = np.array([[np.cos(theta_12), np.sin(theta_12), 0],[-np.sin(theta_12), np.cos(theta_12), 0], [0, 0, 1]])

#U_PMNS = (U_23.dot(U_13)).dot(U_12)
U_PMNS = U_23@U_13@U_12
U_PMNS_conj = np.conjugate(U_PMNS)
U_PMNS_adj = np.transpose(U_PMNS_conj)

lepton_masses = [M_E, M_MU, M_TAU]


# Transition magnetic moment functions
def mag_moment_dirac(mnu):
    return 3*E_QED*G_F*mnu / (8*sqrt(2)*pi**2)

def transition_mag_moment_dirac(mi, mj, i=0, j=1):
    return ((3*E_QED*G_F*mi)/(32*sqrt(2)*pi**2)) * (1 + mj/mi) \
        * np.sum([U_PMNS[l,j]*U_PMNS_conj[l,i]*(lepton_masses[l]/M_W)**2 for l in range(2)])

def transition_mag_moment_sterile(mN, Ualpha):
    # i = mN
    # j = mnu << mN
    return ((3*E_QED*G_F*mN)/(32*sqrt(2)*pi**2)) \
        * np.sum([np.sum([np.real(U_PMNS[l,j])*Ualpha*(lepton_masses[l]/M_W)**2 for l in range(2)])
           for j in range(2)])


### Individual decay widths

# N -> nu gamma
def decay_width_hnl_nu_gamma(mN, Ualpha):
    # 2x for Majorana
    return 2*9*ALPHA * G_F**2 * mN**5 * Ualpha**2 / (512*pi**2)




# N -> pi^0 nu
def decay_width_hnl_pi_nu(mN, Ualpha):
    # 2x for Majorana
    return 2*np.nan_to_num(heaviside(mN - M_PI0, 0.0) * \
        power(Ualpha * G_F * F_PI * (1 - (M_PI0/mN)**2), 2) * mN**3 / (32*pi))




# N -> eta^0 nu
F_ETA = (cos(THETA_8) * ETA_F_8 / sqrt(3)) + (sin(THETA_0) * ETA_F_0 / sqrt(6))
def decay_width_hnl_eta_nu(mN, Ualpha):
    # 2x for Majorana
    return 2*np.nan_to_num(heaviside(mN - M_ETA, 0.0) * \
        power(Ualpha * G_F * F_ETA * (1 - (M_ETA/mN)**2), 2) * mN**3 / (32*pi))




# N -> eta^prime nu
F_ETA_PRIME = (sin(THETA_8) * ETA_F_8 / sqrt(3)) - (cos(THETA_0) * ETA_F_0 / sqrt(6))
def decay_width_hnl_etaprime_nu(mN, Ualpha):
    # 2x for Majorana
    return 2*np.nan_to_num(heaviside(mN - M_ETA_PRIME, 0.0) * \
        power(Ualpha * G_F * F_ETA_PRIME * (1 - (M_ETA_PRIME/mN)**2), 2) * mN**3 / (32*pi))




# N -> nu nu nu
def decay_width_hnl_3nu(mN, Ue, Umu, Utau):
    # x2 for CP-conjugated state
    # 2x for Majorana
    return 2*np.nan_to_num(G_F**2 * mN**5 * (Ue**2 + Umu**2 + Utau**2) / (192*pi**3))




# N -> nu l^+_alpha l^-_beta
def decay_width_hnl_null_mixed(mN, Ualpha, alpha_flavor=0, beta_flavor=0):
    if alpha_flavor not in [0,1,2] or beta_flavor not in [0,1,2]:
        raise Exception("nu_flavor and l_flavor must be in [0,1,2]!")
    if alpha_flavor == beta_flavor:
        return 0.0
    
    ml_a = lepton_masses[alpha_flavor]
    ml_b = lepton_masses[beta_flavor]
    xl = max(ml_a, ml_b) / mN
    xl2 = np.power(xl, 2)
    xl4 = np.power(xl, 4)
    xl6 = np.power(xl, 6)
    xl8 = np.power(xl, 8)

    prefactor = G_F**2 * mN**5 * Ualpha**2 / (192*pi**3) #* np.heaviside(mN - ml_a - ml_b, 0.0)

    # 2x for Majorana
    return 2 * prefactor * np.clip(np.nan_to_num(1 - 8*xl2 + 8*xl6 - xl8 - 12*xl4*log(xl2)),
                                   a_min=0.0, a_max=np.inf)




# N -> nu l+ l-
def decay_width_hnl_null(mN, Ualpha, nu_flavor=0, l_flavor=0):
    if nu_flavor not in [0,1,2] or l_flavor not in [0,1,2]:
        raise Exception("nu_flavor and l_flavor must be in [0,1,2]!")
    
    ml = lepton_masses[l_flavor]
    xl = ml / mN
    xl2 = np.power(xl, 2)
    xl4 = np.power(xl, 4)
    xl6 = np.power(xl, 6)
    log_arg_num = np.clip((1 - 3*xl2 - (1-xl2)*sqrt(1 - 4*xl2)),a_min=1.1e-16,a_max=np.inf)
    log_arg_denom = ((1+sqrt(1-4*xl2))*xl2)
    logterm = log(log_arg_num) - log(log_arg_denom)
    c1 = 0.25*(1- 4*SSW + 8*SSW**2)
    c2 = 0.5*SSW*(2*SSW - 1)
    c3 = 0.25*(1 + 4*SSW + 8*SSW**2)
    c4 = 0.5*SSW*(2*SSW + 1)

    prefactor = np.heaviside(mN - 2*ml, 0.0) * G_F**2 * mN**5 * Ualpha**2 / (192*pi**3)

    delta = nu_flavor == l_flavor

    # 2x for CP conjugated states
    # 2x for Majorana
    return 4*np.nan_to_num(prefactor * ( (c1*(1-delta) + c3*delta) \
                        * ((1 - 14*xl2 - 2*xl4 - 12*xl6)*sqrt(1-4*xl2) + 12*xl4*(xl4-1)*logterm) \
                        + 4*(c2*(1-delta) + c4*delta) * (xl2*(2 + 10*xl2 - 12*xl4)*sqrt(1-4*xl2) + 6*xl4*(1-2*xl2+2*xl4)*logterm)))




# N -> pi^+- l^-+_alpha
def decay_width_hnl_pi_lep(mN, Ualpha, l_flavor=0):
    ml = lepton_masses[l_flavor]
    prefactor = np.heaviside(mN - M_PI - ml, 0.0) * power(G_F * F_PI * Ualpha * V_UD, 2) * mN**3 / (16*pi)
    
    # x2 for CP-conjugated final states
    # 2x for Majorana
    return 4*np.nan_to_num(prefactor * power((1-power(ml/mN, 2))**2 - power(M_PI/mN, 2)*(1 + power(ml/mN, 2)), 2) \
                        * sqrt((1-power((M_PI - ml)/mN, 2))*(1-power((M_PI + ml)/mN, 2))))




# N -> P^+- l^-+_alpha: general charged pseudoscalar meson P
def decay_width_hnl_P_lep(mN, Ualpha, l_flavor=0, fP=F_PI, m_meson=M_PI, ckm=V_UD):
    ml = lepton_masses[l_flavor]
    xP = m_meson/mN
    xl = ml/mN
    prefactor = power(G_F * fP * Ualpha * ckm, 2) * mN**3 / (16*pi) * np.heaviside(mN - m_meson - ml, 0.0)
    
    # x2 for CP-conjugated final states
    # 2x for Majorana
    kallen = sqrt(kallen_alplib(1.0, xP**2, xl**2))

    return prefactor * np.clip(np.nan_to_num(kallen * (1 - xP**2 - xl**2 * (2 + xP**2 - xl**2))),
                               a_min=0.0, a_max=np.inf)




# N -> K^+- l^-+_alpha: general charged pseudoscalar meson P
def decay_width_hnl_K_lep(mN, Ualpha, l_flavor=0):
    return decay_width_hnl_P_lep(mN, Ualpha, l_flavor, F_K, M_K, V_US)




# N -> D^+- l^-+_alpha: general charged pseudoscalar meson P
def decay_width_hnl_D_lep(mN, Ualpha, l_flavor=0):
    return decay_width_hnl_P_lep(mN, Ualpha, l_flavor, F_D, M_D_MESON, V_CD)




# N -> D_s^+- l^-+_alpha: general charged pseudoscalar meson P
def decay_width_hnl_Ds_lep(mN, Ualpha, l_flavor=0):
    return decay_width_hnl_P_lep(mN, Ualpha, l_flavor, F_DS, M_DS_MESON, V_CS)




# N -> V^+- l^-+_alpha
charged_vector_meson_dict = {
    'rho': {'f': F_RHO, 'm': M_RHO_PLUS, 'ckm': V_UD},
    'kstar': {'f':F_KSTAR, 'm': M_KSTAR, 'ckm': V_US},
}
def decay_width_hnl_V_lep(mN, Ualpha, l_flavor=0, meson_species='rho'):
    ml = lepton_masses[l_flavor]
    m_meson = charged_vector_meson_dict[meson_species]['m']
    fV = charged_vector_meson_dict[meson_species]['f']
    ckm = charged_vector_meson_dict[meson_species]['ckm']
    xV = m_meson/mN
    xl = ml/mN
    prefactor = np.heaviside(mN - m_meson - ml, 0.0) \
        * power(G_F * fV * Ualpha * ckm / m_meson, 2) * mN**3 / (16*pi)
    
    kallen = sqrt(kallen_alplib(1.0, xV**2, xl**2))

    return prefactor * np.clip(np.nan_to_num(kallen * ((1 - xV**2)*(1 + 2*xV**2) + xl**2 * (xV**2 + xl**2 - 2))),
                               a_min=0.0, a_max=np.inf)




# N -> V nu
neutral_vector_meson_dict = {
    'rho': {'f': F_RHO, 'm': M_RHO0, 'gv': 1 - 2*SSW},
    'omega': {'f':F_OMEGA, 'm': M_OMEGA, 'gv':  -2*SSW/3},
    'phi': {'f': F_PHI, 'm': M_PHI, 'gv': -sqrt(2)*(0.5 - 2*SSW/3)}
}
def decay_width_hnl_V_nu(mN, Ualpha, meson_species='rho'):
    m_meson = neutral_vector_meson_dict[meson_species]['m']
    fV = neutral_vector_meson_dict[meson_species]['f']
    gV = neutral_vector_meson_dict[meson_species]['gv']
    xV = m_meson/mN
    prefactor = np.heaviside(mN - m_meson, 0.0) \
        * power(G_F * gV * fV * Ualpha / m_meson, 2) * mN**3 / (32*pi)
    
    return prefactor * np.clip(np.nan_to_num((1 + 2*xV**2) * power(1 - xV**2, 2)),
                               a_min=0.0, a_max=np.inf)




# hadronic BRs (taken from data)
ehadr_data = np.genfromtxt("data/hnl_br_e_hadronic_by_mass-MeV.txt")
muhadr_data = np.genfromtxt("data/hnl_br_mu_hadronic_by_mass-MeV.txt")
nuhadr_data = np.genfromtxt("data/hnl_br_nu_hadronic_by_mass-MeV.txt")

def decay_width_e_hadr(mN, Ualpha):
    return power(Ualpha, 2) * np.interp(mN, ehadr_data[:,0], ehadr_data[:,1], left=0.0)




def decay_width_mu_hadr(mN, Ualpha):
    return power(Ualpha, 2) * np.interp(mN, muhadr_data[:,0], muhadr_data[:,1], left=0.0)




def decay_width_nu_hadr(mN, Ualpha):
    return power(Ualpha, 2) * np.interp(mN, nuhadr_data[:,0], nuhadr_data[:,1], left=0.0)




### Summed decay widths
# N -> visible (+ invisible)
def decay_width_visible(mN, Ualpha, flavors=[0,1,2]):
    delta_e = 0 in flavors
    delta_mu = 1 in flavors
    delta_tau = 2 in flavors
    # N -> P nu
    gamma_neutral_mesons = len(flavors)*decay_width_hnl_pi_nu(mN, Ualpha) \
                            + len(flavors)*decay_width_hnl_eta_nu(mN, Ualpha) \
                            + len(flavors)*decay_width_hnl_etaprime_nu(mN, Ualpha)
    # N -> P^+- l^-+
    gamma_charged_mesons = delta_e*decay_width_hnl_pi_lep(mN, Ualpha, l_flavor=0) \
                            + delta_mu*decay_width_hnl_pi_lep(mN, Ualpha, l_flavor=1) \
                            + delta_tau*decay_width_hnl_pi_lep(mN, Ualpha, l_flavor=2) \
                            + delta_e*decay_width_hnl_K_lep(mN, Ualpha, l_flavor=0) \
                            + delta_mu*decay_width_hnl_K_lep(mN, Ualpha, l_flavor=1) \
                            + delta_tau*decay_width_hnl_K_lep(mN, Ualpha, l_flavor=2) \
                            + delta_e*decay_width_hnl_D_lep(mN, Ualpha, l_flavor=0) \
                            + delta_mu*decay_width_hnl_D_lep(mN, Ualpha, l_flavor=1) \
                            + delta_tau*decay_width_hnl_D_lep(mN, Ualpha, l_flavor=2) \
                            + delta_e*decay_width_hnl_Ds_lep(mN, Ualpha, l_flavor=0) \
                            + delta_mu*decay_width_hnl_Ds_lep(mN, Ualpha, l_flavor=1) \
                            + delta_tau*decay_width_hnl_Ds_lep(mN, Ualpha, l_flavor=2)
    # N -> V nu
    gamma_vector_nu = len(flavors)*decay_width_hnl_V_nu(mN, Ualpha, 'rho') \
                        + len(flavors)*decay_width_hnl_V_nu(mN, Ualpha, 'omega') \
                        + len(flavors)*decay_width_hnl_V_nu(mN, Ualpha, 'phi')
    # N -> V^+- l^-+
    gamma_charged_vector_lep = delta_e*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=0, meson_species='rho') \
                    + delta_mu*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=1, meson_species='rho') \
                    + delta_tau*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=2, meson_species='rho') \
                    + delta_e*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=0, meson_species='kstar') \
                    + delta_mu*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=1, meson_species='kstar') \
                    + delta_tau*decay_width_hnl_V_lep(mN, Ualpha, l_flavor=2, meson_species='kstar')
    # l+ l-
    gamma_ll = delta_e*decay_width_hnl_null(mN, Ualpha, nu_flavor=0, l_flavor=0) + delta_mu*decay_width_hnl_null(mN, Ualpha, nu_flavor=1, l_flavor=0) \
            + delta_e*decay_width_hnl_null(mN, Ualpha, nu_flavor=0, l_flavor=1) + delta_mu*decay_width_hnl_null(mN, Ualpha, nu_flavor=1, l_flavor=1) \
            + delta_e*decay_width_hnl_null(mN, Ualpha, nu_flavor=0, l_flavor=2) + delta_mu*decay_width_hnl_null(mN, Ualpha, nu_flavor=1, l_flavor=2) \
            + delta_tau*decay_width_hnl_null(mN, Ualpha, nu_flavor=2, l_flavor=0) + delta_tau*decay_width_hnl_null(mN, Ualpha, nu_flavor=2, l_flavor=1) \
            + delta_tau*decay_width_hnl_null(mN, Ualpha, nu_flavor=2, l_flavor=2)
    
    gamma_ll_mixed = 0.0
    for alpha in flavors:
        gamma_ll_mixed += decay_width_hnl_null_mixed(mN, Ualpha, alpha_flavor=alpha, beta_flavor=0)
        gamma_ll_mixed += decay_width_hnl_null_mixed(mN, Ualpha, alpha_flavor=alpha, beta_flavor=1)
        gamma_ll_mixed += decay_width_hnl_null_mixed(mN, Ualpha, alpha_flavor=alpha, beta_flavor=2)
    
    # hadronic
    #gamma_nu_hadr = decay_width_nu_hadr(mN, Ualpha)
    #gamma_e_hadr = decay_width_e_hadr(mN, Ualpha)
    #gamma_mu_hadr = decay_width_mu_hadr(mN, Ualpha)

    total_gamma = gamma_charged_mesons + gamma_neutral_mesons + gamma_vector_nu \
            + gamma_charged_vector_lep + gamma_ll + gamma_ll_mixed #\
            #+ gamma_nu_hadr + gamma_e_hadr + gamma_mu_hadr
    
    return total_gamma




# N -> anything for a SINGLE alpha coupling
def total_decay_width_hnl(mN, Ualpha, alpha=0):
    delta_e = alpha==0
    delta_mu = alpha==1
    delta_tau = alpha==2
    return decay_width_visible(mN, Ualpha, flavors=[alpha]) + 3*decay_width_hnl_3nu(mN, Ualpha*delta_e, Ualpha*delta_mu, Ualpha*delta_tau)



