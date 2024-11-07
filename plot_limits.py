import numpy as np
from numpy import power, sqrt

import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from matplotlib.lines import Line2D

from scipy.signal import savgol_filter

# import existing limits

# muons
muon_limits = {
    "ATLAS1": np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/ATLAS_1_mu_data.csv", delimiter=","),
    "ATLAS2" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/ATLAS_2_mu_data.csv", delimiter=","),
    "ATLAS3" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/ATLAS_3_mu_data.csv", delimiter=","),
    "bebc" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/BEBC_data.csv", delimiter=","),
    "belle" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/Belle_mu_data.csv", delimiter=","),
    "ccfr" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CCFR_data.csv", delimiter=","),
    "cdhs" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CDHS_data.csv", delimiter=","),
    "charm" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CHARM_mu_data.csv", delimiter=","),
    "CMS" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CMS_mu_data.csv", delimiter=","),
    "delphi" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/DELPHI_data.csv", delimiter=","),
    "e949" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/E949_data.csv", delimiter=","),
    "ewpd" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/EWPD_mu_data.csv", delimiter=","),
    "fmmf" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/FMMF_data.csv", delimiter=","),
    "Higgs" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/Higgs_mu_data.csv", delimiter=","),
    "ic_dc" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/IC_DC_data.csv", delimiter=","),
    "IceCube" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/IceCube_data.csv", delimiter=","),
    "KEK" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/KEK_data.csv", delimiter=","),
    "L3" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/L3_mu_data.csv", delimiter=","),
    "LHCB" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/LHCb_data.csv", delimiter=","),
    "MesonDecays_LNV" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MesonDecays_LNV_mu_data.csv", delimiter=","),
    "MicroBooNE_2" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MicroBooNE_2_data.csv", delimiter=","),
    "MicroBooNE" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MicroBooNE_data.csv", delimiter=","),
    "MiniBooNE" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MiniBooNE_data.csv", delimiter=","),
    "MINOS" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MINOS_data.csv", delimiter=","),
    "MuSpectrum" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/MesonDecays_LNV_mu_data.csv", delimiter=","),
    "NA3" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/NA3_mu_data.csv", delimiter=","),
    "na62" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/NA62_mu_data.csv", delimiter=","),
    "nova" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/NOvA_data.csv", delimiter=","),
    "nutev" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/NuTeV_data.csv", delimiter=","),
    "pienu" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/PIENU_mu_data.csv", delimiter=","),
    "ps191" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/PS191_mu_data.csv", delimiter=","),
    "psi" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/PSI_data.csv", delimiter=","),
    "SBN" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/SBN_data.csv", delimiter=","),
    "SuperK" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/SuperK_mu_data.csv", delimiter=","),
    "SuperKdecay" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/SuperK_mu_decay_data.csv", delimiter=","),
    "T2K" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/T2K_mu_data.csv", delimiter=","),
    "Xray" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/Xray_data.csv", delimiter=","),
    "Xray2" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/Xray_2_data.csv", delimiter=","),
}

muon_limits_astro = {
    "bbn" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/BBN_data.csv", delimiter=","),
    "KoppSN2" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CHARM_mu_data.csv", delimiter=","),
    "KoppSN" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/CHARM_mu_data.csv", delimiter=","),
    "Planck" : np.genfromtxt("data/limits_muon_mixingSq_vs_hnlmass/Planck_data.csv", delimiter=","),
}

# tau mixing
tau_limits = {
    "bdecays": np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/B_decays_data.csv", delimiter=","),
    "bebc" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/BEBC_tau_data.csv", delimiter=","),
    "charm" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/CHARM_tau_data.csv", delimiter=","),
    "charm2" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/CHARM_tau_2_data.csv", delimiter=","),
    "CMS" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/CMS_tau_data.csv", delimiter=","),
    "delphi" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/DELPHI_data.csv", delimiter=","),
    "ewpd" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/EWPD_tau_data.csv", delimiter=","),
    "ic_dc" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/IC_DC_data.csv", delimiter=","),
    "L3" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/L3_tau_data.csv", delimiter=","),
    "nomad" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/NOMAD_data.csv", delimiter=","),
    "nova" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/NOvA_data.csv", delimiter=","),
    "SuperK" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/SuperK_tau_data.csv", delimiter=","),
    "SuperKdecay" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/SuperK_tau_decay_data.csv", delimiter=","),
    "T2K" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/T2K_data.csv", delimiter=","),
    "T2Kdecay" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/T2K_tau_decay_data.csv", delimiter=","),
    "univ" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/tau_universality_data.csv", delimiter=","),
    "xray2" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/Xray_2_data.csv", delimiter=","),
    "xray" : np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/Xray_data.csv", delimiter=","),
    "argoneut": np.genfromtxt("data/limits_tau_mixingSq_vs_hnlmass/ArgoNeuT_data.txt")
}

# future limits
dune_future = np.genfromtxt("data/DUNE_mu_data.csv", delimiter=",")
dune_future_2 = np.genfromtxt("data/DUNE_Indirect_data.csv", delimiter=",")

dq_18POT_Umu = np.genfromtxt("data/darkQuest_1e18_projection_Umu2_vs_mNGeV.txt")
dq_20POT_Umu = np.genfromtxt("data/darkQuest_1e20_projection_Umu2_vs_mNGeV.txt")
dq_18POT_Utau = np.genfromtxt("data/darkQuest_1e18_projection_Utau2_vs_mNGeV.txt")
dq_20POT_Utau = np.genfromtxt("data/darkQuest_1e20_projection_Utau2_vs_mNGeV.txt")


def seesaw_upper(mN):
    return 0.23e-6 / mN

def seesaw_lower(mN):
    return sqrt(7.5e-5) * 1e-6 / mN


def get_Usq_upper_lower_limits(file, color='royalblue', cls=[3.0, 10.0, 100.0],
                                ls=['solid', 'dashed', 'dotted']):
    # takes in mass, coupling, chi2
    dat = np.genfromtxt(file)
    m_unique = np.unique(dat[:,0])
    mass_array = dat[:,0]
    coupling_array = dat[:,1]**2
    chi2_array = dat[:,2]

    fine_mass_array = np.logspace(np.log10(min(m_unique)), np.log10(max(m_unique)), 1000)

    for i, t in enumerate(cls):
        coupling_2s_lower = []
        mass_list_i = []

        for m in m_unique:
            this_coupling_arr = coupling_array[mass_array == m]
            this_chi2_arr = chi2_array[mass_array == m]
            this_chi2_arr -= min(this_chi2_arr)

            searching_2s_lower = True

            for j, g in enumerate(this_coupling_arr):
                if (this_chi2_arr[j] > t) and (searching_2s_lower):
                    coupling_2s_lower.append(g)
                    mass_list_i.append(m)
                    searching_2s_lower = False
                continue
        
        coupling_interp = lambda m: np.interp(m, mass_list_i, coupling_2s_lower, left=1.0e-1, right=1.0e-1)

        coupling_interp = savgol_filter(coupling_interp(fine_mass_array), 51, 1)

        plt.plot(fine_mass_array, coupling_interp, color=color, ls=ls[i])


def get_Usq_contour_from_file(file, color='royalblue', cls=[3.0, 10.0, 100.0], ls=['solid', 'dashed', 'dotted']):
    dat = np.genfromtxt(file)

    MASS, COUPLING = np.meshgrid(np.unique(dat[:,0]),np.unique(dat[:,1]))
    CHI2 = np.reshape(dat[:,2], (np.unique(dat[:,0]).shape[0],np.unique(dat[:,1]).shape[0])).transpose()

    UALPHA2 = COUPLING**2

    return plt.contour(MASS, UALPHA2, CHI2, levels=cls, colors=color, linestyles=ls)


def plot_muon_limits(mass_ratio=2.1):

    for lim in muon_limits:
        dat = muon_limits[lim]
        Usq = power(10,dat[:,1])
        mN = 1e3 * power(10, dat[:,0])
        plt.fill_between(mN, Usq, y2=1.0, alpha=1.0, color='silver')
    
    # cosmology
    dat = muon_limits_astro["bbn"]
    plt.plot()


    if mass_ratio == 2.1:
        c1 = get_Usq_upper_lower_limits("scans/hnl_DUNE_scan_mass-ratio-2_ebrem_umu_highstats.txt", color='royalblue')
        c2 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-2_pbrem_umu.txt", color='indianred')
        c3 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-2_meson_umu.txt", color='g')
    elif mass_ratio == 5:
        c3 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-5_pbrem_umu.txt", color='indianred')
        c1 = get_Usq_upper_lower_limits("scans/hnl_DUNE_scan_mass-ratio-5_ebrem_umu.txt", color='royalblue')
        c3 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-5_meson_umu.txt", color='g')


    # plot seesaw line
    masses = np.logspace(0, 4, 100)
    plt.fill_between(masses, y1=seesaw_lower(masses), y2=seesaw_upper(masses), color='gold', alpha=0.3)

    # future dune
    plt.plot(1e3*power(10,dune_future[:,0]), power(10, dune_future[:,1]), ls='dotted', color='teal', linewidth=1.0)
    plt.text(43.0, 1.2e-8, "DUNE-ND\n(Mixing Only)", rotation=-35.0, color='teal')

    plt.text(30.0, 1e-9, "Seesaw", rotation=-14.0, fontsize=12)
    plt.text(1500.0, 5e-4, "Laboratory Limits\n(Mixing-only)", fontsize=12)

    line_pbrem = Line2D([0], [0], label=r'$p N \to p N Z^\prime$', color='indianred')
    line_ebrem = Line2D([0], [0], label=r'$e^\pm N \to e^\pm N Z^\prime$, $e^+ e^- \to Z^\prime$', color='royalblue')
    line_pion = Line2D([0], [0], label=r'$\pi^0 (\eta^0) \to \gamma Z^\prime$', color='g')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([line_ebrem, line_pbrem, line_pion])
    
    plt.yscale('log')
    plt.xscale('log')
    if mass_ratio == 2.1:
        plt.legend(handles=handles, loc="lower right", framealpha=1, fontsize=10, title="DUNE-ND")
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 2.1$", loc="right", fontsize=14)
    elif mass_ratio == 5:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 5$", loc="right", fontsize=14)
    plt.ylabel(r"$|U_\mu|^2$", fontsize=14)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((10.0, 1e4))
    plt.ylim((1e-12, 1.0e-2))
    plt.tight_layout()
    plt.show()


def plot_muon_limits_sbn(mass_ratio=2.1):

    for lim in muon_limits:
        dat = muon_limits[lim]
        Usq = power(10,dat[:,1])
        mN = 1e3 * power(10, dat[:,0])
        if lim == "MicroBooNE":
            plt.plot(mN, Usq, color='teal', ls='dotted')
        if lim == "MicroBooNE_2":
            plt.plot(mN, Usq, color='teal', ls='dotted')

        plt.fill_between(mN, Usq, y2=1.0, alpha=1.0, color='silver')


    if mass_ratio == 2.1:
        c1 = get_Usq_contour_from_file("scans/hnl_SBND_scan_mass-ratio-2_pbrem_umu.txt", color='indianred', cls=[3.0])
        c3 = get_Usq_contour_from_file("scans/hnl_MUB_scan_mass-ratio-2_pbrem_umu.txt", color='cadetblue', cls=[3.0])
        c2 = get_Usq_contour_from_file("scans/hnl_ICARUS_scan_mass-ratio-2_pbrem_umu.txt", color='orange', cls=[3.0])
        c4 = get_Usq_contour_from_file("scans/hnl_SBND-dump_scan_mass-ratio-2_pbrem_umu.txt",
                                       color='indianred', cls=[3.0*3.0], ls=['dashed'])
        c5 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-2_pbrem_umu.txt", color='g', cls=[3.0])
        c6 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-2_pbrem_umu.txt", color='g', cls=[3.0e-2], ls=['dashed'])
    elif mass_ratio == 5:
        c1 = get_Usq_contour_from_file("scans/hnl_SBND_scan_mass-ratio-5_pbrem_umu.txt", color='indianred', cls=[3.0])
        c3 = get_Usq_contour_from_file("scans/hnl_MUB_scan_mass-ratio-5_pbrem_umu.txt", color='cadetblue', cls=[3.0])
        c2 = get_Usq_contour_from_file("scans/hnl_ICARUS_scan_mass-ratio-5_pbrem_umu.txt", color='orange', cls=[3.0])
        c4 = get_Usq_contour_from_file("scans/hnl_SBND-dump_scan_mass-ratio-5_pbrem_umu.txt",
                                       color='indianred', cls=[3.0*3.0], ls=['dashed'])
        c5 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-5_pbrem_umu.txt", color='g', cls=[3.0])
        c6 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-5_pbrem_umu.txt", color='g', cls=[3.0e-2], ls=['dashed'])
    
    # plot DQ projections
    plt.plot(1e3*dq_18POT_Umu[:,0], dq_18POT_Umu[:,1], linewidth=1.0, color='mediumpurple', ls='dashdot')
    plt.plot(1e3*dq_20POT_Umu[:,0], dq_20POT_Umu[:,1], linewidth=1.0, color='mediumpurple', ls='dotted')

    # plot seesaw line
    masses = np.logspace(0, 4, 100)
    plt.fill_between(masses, y1=seesaw_lower(masses), y2=seesaw_upper(masses), color='gold', alpha=0.3)

    plt.text(1500.0, 1e-11, "Seesaw", rotation=-14.0, fontsize=12)
    plt.text(1750.0, 5e-4, "Laboratory Limits\n(Mixing-only)", fontsize=12)

    # text for microboone
    plt.text(50, 1.0e-7, "MicroBooNE\n(Mixing Only)", rotation=-45.0, color='teal')

    # text for DQ
    plt.text(600.0, 3e-7, "DQ $10^{18}$\n(Mixing Only)", rotation=-25.0, color="mediumpurple")
    plt.text(2000.0, 7e-9, "DQ $10^{20}$\n(Mixing Only)", rotation=90.0, color="mediumpurple")

    line_sbnd = Line2D([0], [0], label=r'SBND', color='indianred')
    line_sbnd_dump = Line2D([0], [0], label=r'SBND Dump Mode', color='indianred', ls='dashed')
    line_mub = Line2D([0], [0], label=r'MicroBooNE', color='cadetblue')
    line_icarus = Line2D([0], [0], label=r'ICARUS-BNB', color='orange')
    line_dq = Line2D([0], [0], label=r'DarkQuest ($10^{18}$ POT)', color='g')
    line_dq_20 = Line2D([0], [0], label=r'DarkQuest ($10^{20}$ POT)', color='g', ls='dashed')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([line_sbnd, line_mub, line_icarus, line_sbnd_dump, line_dq, line_dq_20])
    plt.legend(handles=handles, loc="lower left", framealpha=1, fontsize=10)

    plt.yscale('log')
    plt.xscale('log')
    if mass_ratio == 2.1:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 2.1$", loc="right", fontsize=14)
    elif mass_ratio == 5:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 5$", loc="right", fontsize=14)
    plt.ylabel(r"$|U_\mu|^2$", fontsize=14)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((10.0, 1e4))
    plt.ylim((1e-12, 1.0e-2))
    plt.tight_layout()
    plt.show()


def plot_tau_limits(mass_ratio=2.1):

    for lim in tau_limits:
        dat = tau_limits[lim]
        Usq = power(10,dat[:,1])
        mN = 1e3 * power(10, dat[:,0])

        if lim == "argoneut":
            mN = 1e3*dat[:,0]
            Usq = dat[:,1]

        plt.fill_between(mN, Usq, y2=1.0, alpha=1.0, color='silver')


    if mass_ratio == 2.1:
        c1 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-2_pbrem_utau.txt", color='indianred')
        c2 = get_Usq_upper_lower_limits("scans/hnl_DUNE_scan_mass-ratio-2_ebrem_utau.txt", color='royalblue')
        c3 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-2_meson_utau.txt", color='g')
    elif mass_ratio == 5:
        c1 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-5_pbrem_utau.txt", color='indianred')
        c2 = get_Usq_upper_lower_limits("scans/hnl_DUNE_scan_mass-ratio-5_ebrem_utau.txt", color='royalblue')
        c3 = get_Usq_contour_from_file("scans/hnl_DUNE_scan_mass-ratio-5_meson_utau.txt", color='g')


    # plot seesaw line
    masses = np.logspace(0, 4, 100)
    plt.fill_between(masses, y1=seesaw_lower(masses), y2=seesaw_upper(masses), color='gold', alpha=0.3)

    plt.text(30.0, 1e-9, "Seesaw", rotation=-14.0, fontsize=12)
    plt.text(1500.0, 5e-4, "Laboratory Limits\n(Mixing-only)", fontsize=12)

    line_pbrem = Line2D([0], [0], label=r'$p N \to p N Z^\prime$', color='indianred')
    line_ebrem = Line2D([0], [0], label=r'$e^\pm N \to e^\pm N Z^\prime$, $e^+ e^- \to Z^\prime$', color='royalblue')
    line_pion = Line2D([0], [0], label=r'$\pi^0 (\eta^0) \to \gamma Z^\prime$', color='g')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([line_ebrem, line_pbrem, line_pion])
    #plt.legend(handles=handles, loc="lower right", framealpha=1, fontsize=10, title="DUNE-ND")

    plt.yscale('log')
    plt.xscale('log')
    if mass_ratio == 2.1:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 2.1$", loc="right", fontsize=14)
    elif mass_ratio == 5:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 5$", loc="right", fontsize=14)

    plt.ylabel(r"$|U_\tau|^2$", fontsize=14)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((10.0, 1e4))
    plt.ylim((1e-12, 1.0e-2))
    plt.tight_layout()
    plt.show()


def plot_tau_limits_sbn(mass_ratio=2.1):

    for lim in tau_limits:
        dat = tau_limits[lim]
        Usq = power(10,dat[:,1])
        mN = 1e3 * power(10, dat[:,0])
        plt.fill_between(mN, Usq, y2=1.0, alpha=1.0, color='silver')


    if mass_ratio == 2.1:
        c1 = get_Usq_contour_from_file("scans/hnl_SBND_scan_mass-ratio-2_pbrem_utau.txt", color='indianred', cls=[3.0])
        c3 = get_Usq_contour_from_file("scans/hnl_MUB_scan_mass-ratio-2_pbrem_utau.txt", color='cadetblue', cls=[3.0])
        c2 = get_Usq_contour_from_file("scans/hnl_ICARUS_scan_mass-ratio-2_pbrem_utau.txt", color='orange', cls=[3.0])
        c4 = get_Usq_contour_from_file("scans/hnl_SBND-dump_scan_mass-ratio-2_pbrem_utau.txt",
                                       color='indianred', cls=[3.0*3.0], ls=['dashed'])
        c5 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-2_pbrem_utau.txt", color='g', cls=[3.0])
        c6 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-2_pbrem_utau.txt", color='g', cls=[3.0e-2], ls=['dashed'])
    elif mass_ratio == 5:
        c1 = get_Usq_contour_from_file("scans/hnl_SBND_scan_mass-ratio-5_pbrem_utau.txt", color='indianred', cls=[3.0])
        c3 = get_Usq_contour_from_file("scans/hnl_MUB_scan_mass-ratio-5_pbrem_utau.txt", color='cadetblue', cls=[3.0])
        c2 = get_Usq_contour_from_file("scans/hnl_ICARUS_scan_mass-ratio-5_pbrem_utau.txt", color='orange', cls=[3.0])
        c4 = get_Usq_contour_from_file("scans/hnl_SBND-dump_scan_mass-ratio-5_pbrem_utau.txt",
                                       color='indianred', cls=[3.0*3.0], ls=['dashed'])
        c5 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-5_pbrem_utau.txt", color='g', cls=[3.0])
        c6 = get_Usq_contour_from_file("scans/hnl_DQ_scan_mass-ratio-5_pbrem_utau.txt", color='g', cls=[3.0e-2], ls=['dashed'])

    # plot DQ projections
    plt.plot(1e3*dq_18POT_Utau[:,0], dq_18POT_Utau[:,1], linewidth=1.0, color='mediumpurple', ls='dashdot')
    plt.plot(1e3*dq_20POT_Utau[:,0], dq_20POT_Utau[:,1], linewidth=1.0, color='mediumpurple', ls='dotted')

    # plot seesaw line
    masses = np.logspace(0, 4, 100)
    plt.fill_between(masses, y1=seesaw_lower(masses), y2=seesaw_upper(masses), color='gold', alpha=0.3)

    plt.text(1500.0, 1e-11, "Seesaw", rotation=-14.0, fontsize=12)
    plt.text(1750.0, 5e-4, "Laboratory Limits\n(Mixing-only)", fontsize=12)

    # text for DQ
    if mass_ratio == 2.1:
        plt.text(200.0, 1e-4, "DQ $10^{18}$\n(Mixing Only)", rotation=-32.0, color="mediumpurple")
        plt.text(200.0, 2.5e-6, "DQ $10^{20}$\n(Mixing Only)", rotation=-32.0, color="mediumpurple")
    else:
        plt.text(474.0, 1e-5, "DQ $10^{18}$\n(Mixing Only)", rotation=-30.0, color="mediumpurple")
        plt.text(1750.0, 5e-7, "DQ $10^{20}$\n(Mixing Only)", rotation=90.0, color="mediumpurple")

    line_sbnd = Line2D([0], [0], label=r'SBND', color='indianred')
    line_sbnd_dump = Line2D([0], [0], label=r'SBND Dump Mode', color='indianred', ls='dashed')
    line_mub = Line2D([0], [0], label=r'MicroBooNE', color='cadetblue')
    line_icarus = Line2D([0], [0], label=r'ICARUS-BNB', color='orange')
    line_dq = Line2D([0], [0], label=r'DarkQuest ($10^{18}$ POT)', color='g')
    line_dq_20 = Line2D([0], [0], label=r'DarkQuest ($10^{20}$ POT)', color='g', ls='dashed')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([line_sbnd, line_mub, line_icarus, line_sbnd_dump, line_dq, line_dq_20])
    plt.legend(handles=handles, loc="lower left", framealpha=1, fontsize=10)

    plt.yscale('log')
    plt.xscale('log')
    if mass_ratio == 2.1:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 2.1$", loc="right", fontsize=14)
    elif mass_ratio == 5:
        plt.title(r"$g_{B-L} = 10^{-4}$, $m_{Z^\prime}/m_N = 5$", loc="right", fontsize=14)
    plt.ylabel(r"$|U_\tau|^2$", fontsize=14)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim((10.0, 1e4))
    plt.ylim((1e-12, 1.0e-2))
    plt.tight_layout()
    plt.show()


def main():
    plot_muon_limits(mass_ratio=2.1)
    plot_tau_limits(mass_ratio=2.1)
    plot_muon_limits(mass_ratio=5)
    plot_tau_limits(mass_ratio=5)
    
    #plot_muon_limits_sbn(mass_ratio=2.1)
    #plot_tau_limits_sbn(mass_ratio=2.1)
    #plot_muon_limits_sbn(mass_ratio=5)
    #plot_tau_limits_sbn(mass_ratio=5)


if __name__ == "__main__":
    main()