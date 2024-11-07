from hnl_decays import *
import matplotlib.pyplot as plt
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def plot_limited_brs():
    mn_list = np.logspace(1, np.log10(1400), 200)

    widths_3nu = 3*np.nan_to_num(decay_width_hnl_3nu(mn_list, Ue=1.0, Umu=1.0, Utau=1.0))
    widths_eta = np.nan_to_num(decay_width_hnl_eta_nu(mn_list, Ualpha=1.0))
    widths_pi = 3*np.nan_to_num(decay_width_hnl_pi_nu(mn_list, Ualpha=1.0))
    widths_nuee = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=0))
    widths_numumu = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=1))
    widths_nutautau = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=2) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=2) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=2))
    widths_pi_e = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=0))
    widths_pi_mu = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=1))
    widths_pi_tau = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=2))

    width_emu = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=0, beta_flavor=1)

    widths_hnl_nu_gamma = np.nan_to_num(decay_width_hnl_nu_gamma(mn_list, Ualpha=1.0))



    total_width = widths_3nu + widths_pi + widths_eta + widths_nuee \
        + widths_numumu + widths_nutautau + widths_hnl_nu_gamma + widths_pi_e \
        + widths_pi_mu + widths_pi_tau #+ width_emu


    plt.plot(mn_list, widths_3nu/total_width, label=r"$N \to 3 \nu$")
    plt.plot(mn_list, widths_pi/total_width, label=r"$N \to \nu \pi^0$")
    plt.plot(mn_list, widths_pi_e/total_width, label=r"$N \to \pi^\pm e^\mp$")
    plt.plot(mn_list, widths_pi_mu/total_width, label=r"$N \to \pi^\pm \mu^\mp$")
    plt.plot(mn_list, widths_pi_tau/total_width, label=r"$N \to \pi^\pm \tau^\mp$")
    plt.plot(mn_list, widths_eta/total_width, label=r"$N \to \nu \eta$")
    plt.plot(mn_list, widths_nuee/total_width, label=r"$N \to \nu e^+ e^-$")
    plt.plot(mn_list, widths_numumu/total_width, label=r"$N \to \nu \mu^+ \mu^-$")
    plt.plot(mn_list, widths_nutautau/total_width, label=r"$N \to \nu \tau^+ \tau^-$")
    plt.plot(mn_list, widths_hnl_nu_gamma/total_width, label=r"$N \to \nu \gamma$")
    #plt.plot(mn_list, width_emu/total_width, label=r"$N \to \nu e \mu$")
    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim((0.0, 1400.0))
    plt.ylim((3e-4, 1.1))
    plt.xlabel(r"$m_N$ [MeV]", fontsize=14)
    plt.ylabel(r"$\Gamma_i / \Gamma_{total}$", fontsize=14)
    plt.title(r"$|U_e| = |U_\mu| = |U_\tau| = 1.0$", loc="right", fontsize=14)
    plt.legend(fontsize=12, loc="lower right", ncol=2)
    plt.tight_layout()
    plt.show()
    plt.close()



def plot_all_brs():
    mn_list = np.logspace(1, np.log10(2200), 500)

    Ualpha=1.0

    widths_3nu = 3*np.nan_to_num(decay_width_hnl_3nu(mn_list, Ue=1.0, Umu=1.0, Utau=1.0))
    widths_eta = 3*np.nan_to_num(decay_width_hnl_eta_nu(mn_list, Ualpha=1.0))
    widths_etaprime = 3*np.nan_to_num(decay_width_hnl_etaprime_nu(mn_list, Ualpha=1.0))
    widths_pi = 3*np.nan_to_num(decay_width_hnl_pi_nu(mn_list, Ualpha=1.0))
    widths_nuee = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=0))
    widths_numumu = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=1))

    widths_pi_e = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=0))
    widths_pi_mu = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=1))

    width_emu = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=1, beta_flavor=0)
    width_etau = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=2, beta_flavor=0)

    widths_hnl_nu_gamma = np.nan_to_num(decay_width_hnl_nu_gamma(mn_list, Ualpha=1.0))

    widths_K_e = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=0)
    widths_K_mu = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=1)
    widths_D_e = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=0)
    widths_D_mu = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=1)
    widths_Ds_e = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=0)
    widths_Ds_mu = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=1)

    # N -> V nu
    widths_rho = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'rho')
    widths_omega = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'omega')
    widths_phi = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'phi')

    # N -> V^+- l^-+
    widths_rho_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='rho')
    widths_rho_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='rho')
    widths_kstar_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='kstar')
    widths_kstar_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='kstar')

    widths_nu_hadr = decay_width_nu_hadr(mn_list, Ualpha)
    widths_e_hadr = decay_width_e_hadr(mn_list, Ualpha)
    widths_mu_hadr = decay_width_mu_hadr(mn_list, Ualpha)
    
    # TOTAL
    total_width = decay_width_visible(mn_list, Ualpha, flavors=[0,1,2]) + 3*decay_width_hnl_3nu(mn_list, Ualpha, Ualpha, Ualpha)

    plt.plot(mn_list, widths_nu_hadr, color='darkgray', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_e_hadr, color='darkgray', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_mu_hadr, color='darkgray', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_3nu/total_width, color='k', linewidth=1.0)

    plt.plot(mn_list, widths_rho/total_width, color='magenta', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_omega/total_width, color='magenta', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_phi/total_width, color='magenta', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_rho_e/total_width, color='goldenrod', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_rho_mu/total_width, color='goldenrod', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_kstar_e/total_width, color='navy', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_kstar_mu/total_width, color='navy', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, width_emu/total_width, color='rosybrown', ls='solid', linewidth=1.0)
    plt.plot(mn_list, width_etau/total_width, color='rosybrown', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, widths_pi/total_width, color='crimson', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_eta/total_width, color='crimson', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_etaprime/total_width, color='crimson', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_pi_e/total_width, color='lightseagreen', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_pi_mu/total_width, color='lightseagreen', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_pi_tau/total_width, label=r"$N \to \pi^\pm \tau^\mp$", color='lightseagreen', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_K_e/total_width, color='sienna', linewidth=1.0)
    plt.plot(mn_list, widths_K_mu/total_width, color='sienna', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_D_e/total_width, color='turquoise', linewidth=1.0)
    plt.plot(mn_list, widths_D_mu/total_width, color='turquoise', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_Ds_e/total_width, color='khaki', linewidth=1.0)
    plt.plot(mn_list, widths_Ds_mu/total_width, color='khaki', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, widths_nuee/total_width, color='royalblue', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_numumu/total_width, color='royalblue', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_nutautau/total_width, label=r"$N \to \nu \tau^+ \tau^-$", color='royalblue', ls='dotted')
    
    plt.plot(mn_list, widths_hnl_nu_gamma/total_width, color='forestgreen', linewidth=1.0)
    #plt.plot(mn_list, width_emu/total_width, label=r"$N \to \nu e \mu$")
    
    # Plot labels
    fsize = 12
    plt.text(45.0, 0.4, r"$3 \nu$", rotation=5.0, color='k', fontsize=fsize)
    plt.text(693.0, 0.1, r"$\nu \pi^0$", rotation=-15.0, color='crimson', fontsize=fsize)
    plt.text(610.0, 0.011, r"$\nu \eta$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1407.0, 0.007, r"$\nu \eta^\prime$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1740.0, 0.02, r"$\nu \rho^0$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1338.0, 8e-4, r"$\nu \omega$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1131.0, 0.006, r"$\nu \phi$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1050.0, 0.14, r"$\pi^\pm e^\mp$", rotation=-12.0, color='lightseagreen', fontsize=fsize)
    plt.text(400.0, 0.117, r"$\pi^\pm \mu^\mp$", rotation=0.0, color='lightseagreen', fontsize=fsize)
    plt.text(468.0, 0.04, r"$\nu e^+ e^-$", rotation=10.0, color='royalblue', fontsize=fsize)
    plt.text(390.0, 0.006, r"$\nu \mu^+ \mu^-$", rotation=0.0, color='royalblue', fontsize=fsize)
    plt.text(391.0, 0.002, r"$\nu \gamma$", rotation=10.0, color='forestgreen', fontsize=fsize)
    plt.text(670.0, 0.0046, r"$\rho^\pm e^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(950.0, 0.01, r"$\rho^\pm \mu^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(340.0, 6e-5, r"$K^\pm e^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(605.0, 7e-5, r"$K^\pm \mu^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(1773.0, 2e-3, r"${K^*}^\pm e^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1028.0, 3e-5, r"${K^*}^\pm \mu^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1811.0, 2.3e-5, r"$D^\pm e^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    plt.text(1800.0, 1.3e-5, r"$D^\pm \mu^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    
    plt.text(1910.0, 4e-4, r"$D_s^\pm e^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)
    plt.text(2000.0, 1e-4, r"$D_s^\pm \mu^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)

    plt.text(1720.0, 0.17, r"$\nu$+hadr.", rotation=0.0, color='dimgray', fontsize=fsize)
    plt.text(1127.0, 1e-4, r"$e^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)
    plt.text(1270.0, 1e-4, r"$\mu^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)

    plt.text(303.0, 9e-3, r"$e^\pm \mu^\mp$", rotation=40.0, color='rosybrown', fontsize=fsize)
    plt.text(2060.0, 1.4e-5, r"$e^\pm \tau^\mp$", rotation=40.0, color='rosybrown', fontsize=fsize)


    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim((0.0, 2200.0))
    plt.ylim((1e-5, 1.1))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=16)
    plt.ylabel(r"$\Gamma_i / \Gamma_{\rm total}$", fontsize=16)
    plt.title(r"$|U_e| = |U_\mu| = |U_\tau| = 1.0$", loc="right", fontsize=14)
    #plt.legend(fontsize=12, loc="lower right", ncol=2)
    plt.tight_layout()
    plt.show()
    plt.close()



def plot_all_brs_long():
    mn_list = np.logspace(-1, np.log10(2200), 1200)

    Ualpha=1.0

    widths_3nu = 3*np.nan_to_num(decay_width_hnl_3nu(mn_list, Ue=1.0, Umu=1.0, Utau=1.0))
    widths_eta = 3*np.nan_to_num(decay_width_hnl_eta_nu(mn_list, Ualpha=1.0))
    widths_etaprime = 3*np.nan_to_num(decay_width_hnl_etaprime_nu(mn_list, Ualpha=1.0))
    widths_pi = 3*np.nan_to_num(decay_width_hnl_pi_nu(mn_list, Ualpha=1.0))
    widths_nuee = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=0))
    widths_numumu = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=1))

    widths_pi_e = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=0))
    widths_pi_mu = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=1))

    width_emu = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=1, beta_flavor=0)
    width_etau = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=2, beta_flavor=0)

    widths_hnl_nu_gamma = np.nan_to_num(decay_width_hnl_nu_gamma(mn_list, Ualpha=1.0))

    widths_K_e = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=0)
    widths_K_mu = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=1)
    widths_D_e = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=0)
    widths_D_mu = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=1)
    widths_Ds_e = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=0)
    widths_Ds_mu = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=1)

    # N -> V nu
    widths_rho = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'rho')
    widths_omega = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'omega')
    widths_phi = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'phi')

    # N -> V^+- l^-+
    widths_rho_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='rho')
    widths_rho_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='rho')
    widths_kstar_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='kstar')
    widths_kstar_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='kstar')
    
    # TOTAL
    total_width = decay_width_visible(mn_list, Ualpha, flavors=[0,1,2]) + 3*decay_width_hnl_3nu(mn_list, Ualpha, Ualpha, Ualpha)

    # taken from digitizing 2007.03701, approximating total width without multi hadron final states (3nu dominates!)
    widths_nu_hadr = decay_width_nu_hadr(mn_list, Ualpha) * total_width
    widths_e_hadr = decay_width_e_hadr(mn_list, Ualpha) * total_width
    widths_mu_hadr = decay_width_mu_hadr(mn_list, Ualpha) * total_width

    total_width += widths_nu_hadr + widths_e_hadr + widths_mu_hadr  # add back in the hadronic widths.

    plt.figure(figsize=(14,5))

    plt.plot(mn_list, widths_nu_hadr/total_width, color='darkgray', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_e_hadr/total_width, color='darkgray', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_mu_hadr/total_width, color='darkgray', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_3nu/total_width, color='k', linewidth=1.0)

    plt.plot(mn_list, widths_rho/total_width, color='magenta', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_omega/total_width, color='magenta', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_phi/total_width, color='magenta', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_rho_e/total_width, color='goldenrod', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_rho_mu/total_width, color='goldenrod', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_kstar_e/total_width, color='navy', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_kstar_mu/total_width, color='navy', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, width_emu/total_width, color='rosybrown', ls='solid', linewidth=1.0)
    plt.plot(mn_list, width_etau/total_width, color='rosybrown', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, widths_pi/total_width, color='crimson', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_eta/total_width, color='crimson', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_etaprime/total_width, color='crimson', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_pi_e/total_width, color='lightseagreen', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_pi_mu/total_width, color='lightseagreen', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_pi_tau/total_width, label=r"$N \to \pi^\pm \tau^\mp$", color='lightseagreen', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, widths_K_e/total_width, color='sienna', linewidth=1.0)
    plt.plot(mn_list, widths_K_mu/total_width, color='sienna', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_D_e/total_width, color='turquoise', linewidth=1.0)
    plt.plot(mn_list, widths_D_mu/total_width, color='turquoise', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_Ds_e/total_width, color='khaki', linewidth=1.0)
    plt.plot(mn_list, widths_Ds_mu/total_width, color='khaki', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, widths_nuee/total_width, color='royalblue', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_numumu/total_width, color='royalblue', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_nutautau/total_width, label=r"$N \to \nu \tau^+ \tau^-$", color='royalblue', ls='dotted')
    
    plt.plot(mn_list, widths_hnl_nu_gamma/total_width, color='forestgreen', linewidth=1.0)
    #plt.plot(mn_list, width_emu/total_width, label=r"$N \to \nu e \mu$")
    
    # Plot labels
    fsize = 12
    plt.text(45.0, 0.4, r"$3 \nu$", rotation=0.0, color='k', fontsize=fsize)
    plt.text(693.0, 0.1, r"$\nu \pi^0$", rotation=-15.0, color='crimson', fontsize=fsize)
    plt.text(610.0, 0.011, r"$\nu \eta$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1407.0, 0.005, r"$\nu \eta^\prime$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1740.0, 0.02, r"$\nu \rho^0$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(812.0, 3e-4, r"$\nu \omega$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1131.0, 0.006, r"$\nu \phi$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1050.0, 0.14, r"$\pi^\pm e^\mp$", rotation=-12.0, color='lightseagreen', fontsize=fsize)
    plt.text(400.0, 0.117, r"$\pi^\pm \mu^\mp$", rotation=0.0, color='lightseagreen', fontsize=fsize)
    plt.text(45.0, 0.09, r"$\nu e^+ e^-$", rotation=0.0, color='royalblue', fontsize=fsize)
    plt.text(390.0, 0.006, r"$\nu \mu^+ \mu^-$", rotation=0.0, color='royalblue', fontsize=fsize)
    plt.text(45.0, 0.004, r"$\nu \gamma$", rotation=0.0, color='forestgreen', fontsize=fsize)
    plt.text(750.0, 0.0046, r"$\rho^\pm e^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(900.0, 4.7e-3, r"$\rho^\pm \mu^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(450.0, 6e-5, r"$K^\pm e^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(605.0, 7e-5, r"$K^\pm \mu^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(1773.0, 1.5e-3, r"${K^*}^\pm e^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1028.0, 3e-5, r"${K^*}^\pm \mu^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1911.0, 2.3e-5, r"$D^\pm e^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    plt.text(2020.0, 1.3e-5, r"$D^\pm \mu^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    
    plt.text(2000.0, 4e-4, r"$D_s^\pm e^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)
    plt.text(2101.0, 2e-4, r"$D_s^\pm \mu^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)

    plt.text(1720.0, 0.11, r"$\nu$+hadr.", rotation=0.0, color='dimgray', fontsize=fsize)
    plt.text(1160.0, 1e-4, r"$e^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)
    plt.text(1270.0, 1e-4, r"$\mu^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)

    plt.text(303.0, 9e-3, r"$\nu e^\pm \mu^\mp$", rotation=25.0, color='rosybrown', fontsize=fsize)
    plt.text(2100.0, 2.2e-5, r"$\nu e^\pm \tau^\mp$", rotation=40.0, color='rosybrown', fontsize=fsize)


    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim((0.0, 2200.0))
    plt.ylim((1e-5, 1.1))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=16)
    plt.ylabel(r"Branching Ratio $\Gamma(N \to X) / \Gamma_{\rm N, total}$", fontsize=16)
    plt.title(r"$|U_e| = |U_\mu| = |U_\tau| = 1.0$", loc="right", fontsize=14)
    #plt.legend(fontsize=12, loc="lower right", ncol=2)
    plt.tight_layout()
    plt.show()
    plt.close()



def plot_relevant_brs_long():
    mn_list = np.logspace(-1, np.log10(2200), 1200)

    Ualpha=1.0

    widths_3nu = 3*np.nan_to_num(decay_width_hnl_3nu(mn_list, Ue=1.0, Umu=1.0, Utau=1.0))
    widths_eta = 3*np.nan_to_num(decay_width_hnl_eta_nu(mn_list, Ualpha=1.0))
    widths_etaprime = 3*np.nan_to_num(decay_width_hnl_etaprime_nu(mn_list, Ualpha=1.0))
    widths_pi = 3*np.nan_to_num(decay_width_hnl_pi_nu(mn_list, Ualpha=1.0))
    widths_nuee = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=0) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=0))
    widths_numumu = np.nan_to_num(decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=0, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=1, l_flavor=1) \
                    + decay_width_hnl_null(mn_list, Ualpha=1.0, nu_flavor=2, l_flavor=1))

    widths_pi_e = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=0))
    widths_pi_mu = np.nan_to_num(decay_width_hnl_pi_lep(mn_list, Ualpha=1.0, l_flavor=1))

    width_emu = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=1, beta_flavor=0)
    width_etau = decay_width_hnl_null_mixed(mn_list, Ualpha=1.0, alpha_flavor=2, beta_flavor=0)

    widths_hnl_nu_gamma = np.nan_to_num(decay_width_hnl_nu_gamma(mn_list, Ualpha=1.0))

    widths_K_e = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=0)
    widths_K_mu = decay_width_hnl_K_lep(mn_list, Ualpha, l_flavor=1)
    widths_D_e = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=0)
    widths_D_mu = decay_width_hnl_D_lep(mn_list, Ualpha, l_flavor=1)
    widths_Ds_e = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=0)
    widths_Ds_mu = decay_width_hnl_Ds_lep(mn_list, Ualpha, l_flavor=1)

    # N -> V nu
    widths_rho = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'rho')
    widths_omega = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'omega')
    widths_phi = 3*decay_width_hnl_V_nu(mn_list, Ualpha, 'phi')

    # N -> V^+- l^-+
    widths_rho_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='rho')
    widths_rho_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='rho')
    widths_kstar_e = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=0, meson_species='kstar')
    widths_kstar_mu = decay_width_hnl_V_lep(mn_list, Ualpha, l_flavor=1, meson_species='kstar')

    widths_nu_hadr = decay_width_nu_hadr(mn_list, Ualpha)
    widths_e_hadr = decay_width_e_hadr(mn_list, Ualpha)
    widths_mu_hadr = decay_width_mu_hadr(mn_list, Ualpha)
    
    # TOTAL
    total_width = decay_width_visible(mn_list, Ualpha, flavors=[0,1,2]) + 3*decay_width_hnl_3nu(mn_list, Ualpha, Ualpha, Ualpha)

    plt.figure(figsize=(14,5))

    br_thresh = 0.05*total_width

    plt.plot(mn_list, widths_nu_hadr*(widths_nu_hadr > br_thresh/total_width), color='darkgray', ls='solid', linewidth=1.0)
    plt.plot(mn_list, widths_e_hadr*(widths_e_hadr > br_thresh/total_width), color='darkgray', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, widths_mu_hadr*(widths_mu_hadr > br_thresh/total_width), color='darkgray', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, (widths_3nu > br_thresh)*widths_3nu/total_width, color='k', linewidth=1.0)

    plt.plot(mn_list, (widths_rho > br_thresh)*widths_rho/total_width, color='magenta', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_omega > br_thresh)*widths_omega/total_width, color='magenta', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, (widths_phi > br_thresh)*widths_phi/total_width, color='magenta', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, (widths_rho_e > br_thresh)*widths_rho_e/total_width, color='goldenrod', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_rho_mu > br_thresh)*widths_rho_mu/total_width, color='goldenrod', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, (widths_kstar_e > br_thresh)*widths_kstar_e/total_width, color='navy', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_kstar_mu > br_thresh)*widths_kstar_mu/total_width, color='navy', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, (width_emu > br_thresh)*width_emu/total_width, color='rosybrown', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (width_etau > br_thresh)*width_etau/total_width, color='rosybrown', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, (widths_pi > br_thresh)*widths_pi/total_width, color='crimson', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_eta > br_thresh)*widths_eta/total_width, color='crimson', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, (widths_etaprime > br_thresh)*widths_etaprime/total_width, color='crimson', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, (widths_pi_e > br_thresh)*widths_pi_e/total_width, color='lightseagreen', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_pi_mu > br_thresh)*widths_pi_mu/total_width, color='lightseagreen', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_pi_tau/total_width, label=r"$N \to \pi^\pm \tau^\mp$", color='lightseagreen', ls='dotted', linewidth=1.0)

    plt.plot(mn_list, (widths_K_e > br_thresh)*widths_K_e/total_width, color='sienna', linewidth=1.0)
    plt.plot(mn_list, (widths_K_mu > br_thresh)*widths_K_mu/total_width, color='sienna', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, (widths_D_e > br_thresh)*widths_D_e/total_width, color='turquoise', linewidth=1.0)
    plt.plot(mn_list, (widths_D_mu > br_thresh)*widths_D_mu/total_width, color='turquoise', ls='dashed', linewidth=1.0)
    plt.plot(mn_list, (widths_Ds_e > br_thresh)*widths_Ds_e/total_width, color='khaki', linewidth=1.0)
    plt.plot(mn_list, (widths_Ds_mu > br_thresh)*widths_Ds_mu/total_width, color='khaki', ls='dashed', linewidth=1.0)

    plt.plot(mn_list, (widths_nuee > br_thresh)*widths_nuee/total_width, color='royalblue', ls='solid', linewidth=1.0)
    plt.plot(mn_list, (widths_numumu > br_thresh)*widths_numumu/total_width, color='royalblue', ls='dashed', linewidth=1.0)
    #plt.plot(mn_list, widths_nutautau/total_width, label=r"$N \to \nu \tau^+ \tau^-$", color='royalblue', ls='dotted')
    
    plt.plot(mn_list, (widths_hnl_nu_gamma > br_thresh)*widths_hnl_nu_gamma/total_width, color='forestgreen', linewidth=1.0)
    #plt.plot(mn_list, width_emu/total_width, label=r"$N \to \nu e \mu$")
    
    # Plot labels
    fsize = 12
    plt.text(45.0, 0.4, r"$3 \nu$", rotation=5.0, color='k', fontsize=fsize)
    plt.text(693.0, 0.1, r"$\nu \pi^0$", rotation=-15.0, color='crimson', fontsize=fsize)
    plt.text(610.0, 0.011, r"$\nu \eta$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1407.0, 0.007, r"$\nu \eta^\prime$", rotation=0.0, color='crimson', fontsize=fsize)
    plt.text(1740.0, 0.02, r"$\nu \rho^0$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1338.0, 8e-4, r"$\nu \omega$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1131.0, 0.006, r"$\nu \phi$", rotation=0.0, color='magenta', fontsize=fsize)
    plt.text(1050.0, 0.14, r"$\pi^\pm e^\mp$", rotation=-12.0, color='lightseagreen', fontsize=fsize)
    plt.text(400.0, 0.117, r"$\pi^\pm \mu^\mp$", rotation=0.0, color='lightseagreen', fontsize=fsize)
    plt.text(468.0, 0.04, r"$\nu e^+ e^-$", rotation=10.0, color='royalblue', fontsize=fsize)
    plt.text(390.0, 0.006, r"$\nu \mu^+ \mu^-$", rotation=0.0, color='royalblue', fontsize=fsize)
    plt.text(391.0, 0.002, r"$\nu \gamma$", rotation=10.0, color='forestgreen', fontsize=fsize)
    plt.text(670.0, 0.0046, r"$\rho^\pm e^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(950.0, 0.01, r"$\rho^\pm \mu^\mp$", rotation=0.0, color='goldenrod', fontsize=fsize)
    plt.text(450.0, 6e-5, r"$K^\pm e^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(605.0, 7e-5, r"$K^\pm \mu^\mp$", rotation=0.0, color='sienna', fontsize=fsize)
    plt.text(1773.0, 2e-3, r"${K^*}^\pm e^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1028.0, 3e-5, r"${K^*}^\pm \mu^\mp$", rotation=0.0, color='navy', fontsize=fsize)
    plt.text(1911.0, 2.3e-5, r"$D^\pm e^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    plt.text(2020.0, 1.3e-5, r"$D^\pm \mu^\mp$", rotation=0.0, color='turquoise', fontsize=fsize)
    
    plt.text(2000.0, 4e-4, r"$D_s^\pm e^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)
    plt.text(2101.0, 2e-4, r"$D_s^\pm \mu^\mp$", rotation=0.0, color='darkkhaki', fontsize=fsize)

    plt.text(1720.0, 0.17, r"$\nu$+hadr.", rotation=0.0, color='dimgray', fontsize=fsize)
    plt.text(1160.0, 1e-4, r"$e^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)
    plt.text(1270.0, 1e-4, r"$\mu^\pm$+hadr.", rotation=90.0, color='dimgray', fontsize=fsize)

    plt.text(303.0, 9e-3, r"$e^\pm \mu^\mp$", rotation=25.0, color='rosybrown', fontsize=fsize)
    plt.text(2140.0, 4e-5, r"$e^\pm \tau^\mp$", rotation=40.0, color='rosybrown', fontsize=fsize)


    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim((0.0, 2200.0))
    plt.ylim((1e-5, 1.1))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=16)
    plt.ylabel(r"Branching Ratio $\Gamma(N \to X) / \Gamma_{\rm N, total}$", fontsize=16)
    plt.title(r"$|U_e| = |U_\mu| = |U_\tau| = 1.0$", loc="right", fontsize=14)
    #plt.legend(fontsize=12, loc="lower right", ncol=2)
    plt.tight_layout()
    plt.show()
    plt.close()




def plot_visible_decay_width():
    mn_list = np.logspace(1, np.log10(2200), 500)

    Ualpha=1.0

    visible_width_0 = decay_width_visible(mn_list, Ualpha, flavors=[0])
    visible_width_1 = decay_width_visible(mn_list, Ualpha, flavors=[1])
    visible_width_2 = decay_width_visible(mn_list, Ualpha, flavors=[2])

    total_width_0 = total_decay_width_hnl(mn_list, Ualpha, alpha=0)
    total_width_1 = total_decay_width_hnl(mn_list, Ualpha, alpha=1)
    total_width_2 = total_decay_width_hnl(mn_list, Ualpha, alpha=2)

    plt.plot(mn_list, visible_width_0/total_width_0, color='royalblue', ls='solid', linewidth=1.0, label=r"$|U_e| = 1.0$")
    plt.plot(mn_list, visible_width_1/total_width_1, color='crimson', ls='solid', linewidth=1.0, label=r"$|U_\mu| = 1.0$")
    plt.plot(mn_list, visible_width_2/total_width_2, color='goldenrod', ls='solid', linewidth=1.0, label=r"$|U_\tau| = 1.0$")

    plt.yscale('log')
    #plt.xscale('log')
    plt.xlim((0.0, 2200.0))
    plt.ylim((7e-2, 1.1))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(r"$m_N$ [MeV]", fontsize=16)
    plt.ylabel(r"$\Gamma_{\rm visible} / \Gamma_{\rm total}$", fontsize=16)
    plt.legend(fontsize=12, loc="lower right")
    plt.tight_layout()
    plt.show()
    plt.close()




def main():
    #plot_all_brs()

    plot_all_brs_long()

    #plot_relevant_brs_long()

    #plot_visible_decay_width()


if __name__ == "__main__":
    main()



