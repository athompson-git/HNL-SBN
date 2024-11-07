import sys
sys.path.append("../")

from alplib.fmath import *
from alplib.constants import *

import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

frho = [0.616, 0.223, -0.339]
fomega = [1.011, -0.881, 0.369]

M_RHO = 770.0
M_OMEGA = 782.6
M_OMEGA_1420 = 1420.0

m_vmesons = [M_RHO, M_OMEGA, M_OMEGA_1420]
widths_mesons = [147.4, 8.68, 290.0]
Lambda_cutoff = 1500.0  # [1, 2] GeV, central value = 1.5 GeV

def Hfunc(z, pT, mD):
    return power(pT, 2) + power(z*M_P, 2) + (1-z)*mD**2

def fppD(p2):
    return power(Lambda_cutoff, 4) / (power(Lambda_cutoff, 4) + power(p2 - M_P**2, 2))

def fVsq(kV2, p2):
    fV = np.sum([frho[i] * m_vmesons[i]**2 \
                   / (m_vmesons[i]**2 - kV2 - 1j * m_vmesons[i] * widths_mesons[i]) \
                    for i in range(3)])
    return np.real(fV * np.conjugate(fV)) * fppD(p2)**2

def wV(z, pT, mV, g_V):
    return (g_V**2 / (8*pi**2)) * (1/Hfunc(z, pT, mV)) * fVsq(mV**2, M_P**2 - Hfunc(z, pT, mV)/z) \
                    * (z - z*(1-z)*(2*M_P**2 + mV**2)/Hfunc(z, pT, mV) \
                        + Hfunc(z, pT, mV)/(2*z*mV**2))

def wS(z, pT, mS, g_V):
    return (g_V**2 / (8*pi**2)) * (1/Hfunc(z, pT, mS)) * fVsq(mS**2, M_P**2 - Hfunc(z, pT, mS)/z) \
                    * (z + z*(1-z)*(4*M_P**2 - mS**2)/Hfunc(z, pT, mS))

def sigma_NSD(s):
    # s in GeV
    # returns in mb
    return 1.76 + 19.8*power(s, 0.057)

def dsigma(z, pT, Ep, mV, g_V, rep="vector"):
    # returns dsigma / dz dpT in mb / MeV
    s = (2*M_P**2 + 2*Ep*M_P)*1e-6
    phase_space = np.heaviside(0.2 - Hfunc(z, pT, mV) / (4*z*power((1-z)*Ep,2)), 0.0)

    # convert wV from MeV^-1 to mb
    if rep == "vector":
        return np.clip(2*pT*wV(z, pT, mV, g_V) * sigma_NSD(s) \
                                    * phase_space, a_min=0.0, a_max=np.inf)
    elif rep == "scalar":
        return np.clip(2*pT*wS(z, pT, mV, g_V) * sigma_NSD(s) \
                                    * phase_space, a_min=0.0, a_max=np.inf)
    else:
        raise Exception("rep = {} not found in [scalar, vector]".format(rep))

def sigma_pbrem(Ep, mV, g_V, n_samples=10000):
    z_rnd = 10**np.random.uniform(-3,0,n_samples)
    pt_rnd = Ep*power(10, np.random.uniform(-6,np.log10(0.2),n_samples))
    return (np.log(10)*np.log(10)*(1+3)*(np.log10(0.2) + 6)) \
        * np.sum(z_rnd * pt_rnd * dsigma(z_rnd, pt_rnd, Ep, mV, g_V)) / n_samples


def pbrem_mc_momentum_and_weights(Ep, mV, g_V, n_samples=10000):
    z_rnd = 10**np.random.uniform(-3,0,n_samples)
    pt_rnd = Ep*power(10, np.random.uniform(-6,np.log10(0.2),n_samples))
    mc_vol = (np.log(10)*np.log(10)*(1+3)*(np.log10(0.2) + 6)) / n_samples
    weights = mc_vol* z_rnd * pt_rnd * dsigma(z_rnd, pt_rnd, Ep, mV, g_V)
    return z_rnd, pt_rnd, weights


def plot_dsiga_2d_dist():
    Ep0 = 120.0e3
    z_pts = np.logspace(-1.5, 0, 100) # np.linspace(0.0001, 0.1, 100)
    pt_pts = np.logspace(-4, 0, 100) # Ep0*np.linspace(0.0001, 0.1, 100)
    Z, PT = np.meshgrid(z_pts, pt_pts)
    SIGMA = np.zeros_like(Z)
    for i in range(z_pts.shape[0]):
        for j in range(100):
            SIGMA[i, j] = dsigma(Z[i,j], Ep0*PT[i,j], Ep0, 10000.0, 1.0)

    sigma_min = min(SIGMA.flatten())
    sigma_max = max(SIGMA.flatten())
    print(sigma_min, sigma_max)
    levels_sigma = np.linspace(sigma_min, sigma_max, 50)

    fig, ax = plt.subplots()
    cs = ax.contourf(Z, PT, SIGMA, locator=ticker.LogLocator())
    cbar = fig.colorbar(cs, label=r"$d\sigma/(dz dp_T)$")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r"$z$", fontsize=16)
    plt.ylabel(r"$p_T / E_p$", fontsize=16)
    plt.title(r"$m_V = 10$ GeV", loc="right", fontsize=16)
    plt.show()
    plt.close()


def plot_total_xs():

    """
    ep_list = np.linspace(5000.0, 120000.0, 20)

    mV0 = 1000.0
    sigma_list = np.array([sigma_pbrem(ep, mV0, 1.0, n_samples=1000000) for ep in ep_list])

    plt.plot(ep_list, sigma_list)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r"$E_p$ [MeV]")
    plt.ylabel(r"$\sigma$ [mb]")
    plt.show()
    """

    mv_list = np.linspace(0.0, 2000.0, 400)
    sigma_list_mV = np.array([sigma_pbrem(120.0e3, mV, 0.3, n_samples=200000) for mV in mv_list])

    print(sigma_list_mV)

    plt.plot(mv_list, sigma_list_mV)
    plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel(r"$m_V$ [MeV]")
    plt.ylabel(r"$\sigma$ [mb]")
    plt.show()



def plot_kinematics():
    p_proton = 120.0e3
    z_100MeV, pt_100MeV, wgts_100MeV = pbrem_mc_momentum_and_weights(p_proton, 100.0, 1.0, n_samples=100000)
    z_500MeV, pt_500MeV, wgts_500MeV = pbrem_mc_momentum_and_weights(p_proton, 500.0, 1.0, n_samples=100000)
    z_1GeV, pt_1GeV, wgts_1GeV = pbrem_mc_momentum_and_weights(p_proton, 1000.0, 1.0, n_samples=100000)
    z_5GeV, pt_5GeV, wgts_5GeV = pbrem_mc_momentum_and_weights(p_proton, 5000.0, 1.0, n_samples=100000)

    ptotal_100MeV = sqrt((p_proton*z_100MeV)**2 + pt_100MeV**2)
    ptotal_500MeV = sqrt((p_proton*z_500MeV)**2 + pt_500MeV**2)
    ptotal_1GeV = sqrt((p_proton*z_1GeV)**2 + pt_1GeV**2)
    ptotal_5GeV = sqrt((p_proton*z_5GeV)**2 + pt_5GeV**2)

    theta_100MeV = np.arccos(p_proton*z_100MeV/ptotal_100MeV)
    theta_500MeV = np.arccos(p_proton*z_500MeV/ptotal_500MeV)
    theta_1GeV = np.arccos(p_proton*z_1GeV/ptotal_1GeV)
    theta_5GeV = np.arccos(p_proton*z_5GeV/ptotal_5GeV)

    # plot thetas
    rad2deg = 180.0/np.pi
    theta_bins = np.linspace(0.0, 90.0, 500)
    plt.hist(theta_100MeV*rad2deg, weights=wgts_100MeV, bins=theta_bins, histtype='step', label=r"$m_V = 100$ MeV")
    plt.hist(theta_500MeV*rad2deg, weights=wgts_500MeV, bins=theta_bins, histtype='step', label=r"$m_V = 500$ MeV")
    plt.hist(theta_1GeV*rad2deg, weights=wgts_1GeV, bins=theta_bins, histtype='step', label=r"$m_V = 1$ GeV")
    plt.hist(theta_5GeV*rad2deg, weights=wgts_5GeV, bins=theta_bins, histtype='step', label=r"$m_V = 5$ GeV")
    plt.yscale('log')
    plt.xlabel(r"$\theta_z$ [deg]", fontsize=14)
    plt.ylabel(r"Vector Flux Counts", fontsize=14)
    plt.xlim((0.0, 90.0))
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    theta_bins_log = np.logspace(-3, np.log10(90.0), 100)
    plt.hist(theta_100MeV*rad2deg, weights=wgts_100MeV, bins=theta_bins_log, histtype='step', label=r"$m_V = 100$ MeV")
    plt.hist(theta_500MeV*rad2deg, weights=wgts_500MeV, bins=theta_bins_log, histtype='step', label=r"$m_V = 500$ MeV")
    plt.hist(theta_1GeV*rad2deg, weights=wgts_1GeV, bins=theta_bins_log, histtype='step', label=r"$m_V = 1$ GeV")
    plt.hist(theta_5GeV*rad2deg, weights=wgts_5GeV, bins=theta_bins_log, histtype='step', label=r"$m_V = 5$ GeV")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r"$\theta_z$ [deg]", fontsize=14)
    plt.ylabel(r"Vector Flux Counts", fontsize=14)
    plt.xlim((1e-3, 90.0))
    plt.legend()
    plt.tight_layout()
    plt.show()

    # plot momenta
    p_bins = np.logspace(0, np.log10(p_proton), 100)
    plt.hist(ptotal_100MeV, weights=wgts_100MeV, bins=p_bins, histtype='step', label=r"$m_V = 100$ MeV")
    plt.hist(ptotal_500MeV, weights=wgts_500MeV, bins=p_bins, histtype='step', label=r"$m_V = 500$ MeV")
    plt.hist(ptotal_1GeV, weights=wgts_1GeV, bins=p_bins, histtype='step', label=r"$m_V = 1$ GeV")
    plt.hist(ptotal_5GeV, weights=wgts_5GeV, bins=p_bins, histtype='step', label=r"$m_V = 5$ GeV")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r"$p_V$ [MeV]", fontsize=14)
    plt.ylabel(r"Vector Flux Counts", fontsize=14)
    plt.xlim((1.0, 120000.0))
    plt.legend()
    plt.tight_layout()
    plt.show()



def main():
    plot_total_xs()
    #plot_kinematics()



if __name__ == "__main__":
    main()