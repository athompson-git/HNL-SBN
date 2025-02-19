import sys
sys.path.append("../")

from alplib.fluxes import *
from hnl_decays import *
from hnl_flux import *

pi0_flux_per_1e5POT = np.genfromtxt("data/dune_target_pi0_4vectors_1e5POT.txt")
eta_flux_per_1e5POT = np.genfromtxt("data/dune_target_eta0_4vectors_1e5POT.txt")

# trim to 20 mrad
pi0_theta = np.arccos(pi0_flux_per_1e5POT[:,2] / np.sqrt(pi0_flux_per_1e5POT[:,0]**2 + pi0_flux_per_1e5POT[:,1]**2 + pi0_flux_per_1e5POT[:,2]**2))
eta_theta = np.arccos(eta_flux_per_1e5POT[:,2] / np.sqrt(eta_flux_per_1e5POT[:,0]**2 + eta_flux_per_1e5POT[:,1]**2 + eta_flux_per_1e5POT[:,2]**2))

pi0_flux_per_1e5POT = 1e3*pi0_flux_per_1e5POT[pi0_theta < 2e-2]
eta_flux_per_1e5POT = 1e3*eta_flux_per_1e5POT[eta_theta < 2e-2]

# add the mass energy back in
pi0_p2 = pi0_flux_per_1e5POT[:,0]**2 + pi0_flux_per_1e5POT[:,1]**2 + pi0_flux_per_1e5POT[:,2]**2
pi0_T2 = pi0_flux_per_1e5POT[:,3]**2
pi0_mass = (pi0_p2 - pi0_T2)/(2*pi0_flux_per_1e5POT[:,3])
pi0_flux_per_1e5POT[:,3] += pi0_mass

eta_p2 = eta_flux_per_1e5POT[:,0]**2 + eta_flux_per_1e5POT[:,1]**2 + eta_flux_per_1e5POT[:,2]**2
eta_T2 = eta_flux_per_1e5POT[:,3]**2
eta_mass = (eta_p2 - eta_T2)/(2*eta_flux_per_1e5POT[:,3])
eta_flux_per_1e5POT[:,3] += eta_mass


def flux_generator(m_HNL, m_Zp, gBL, Ualpha, ns=10000):
    flux = FluxHNLFromElectronPositron(axion_mass=m_Zp, mixing_angle=Ualpha, hnl_mass=m_HNL, coupling_BL=gBL, n_samples=ns)

    #resonant_flux.simulate()
    print("simulating res v...")
    flux.simulate()
    flux.decay_to_hnl()
    flux.propagate(Ualpha)

    return np.array(flux.hnl_energy), np.array(flux.hnl_angle), np.array(flux.hnl_flux)


def parameter_scan(mzp_to_mN_ratio=2.1, out_file="scans/hnl_DUNE_scan_mass-ratio-2.1_brem.txt", flavor=1):

    min_pbrem_HNL_mass = 250.0/mzp_to_mN_ratio

    #mN_array = np.logspace(np.log10(min_pbrem_HNL_mass), 4, 60)
    #mN_array = np.logspace(1, np.log10(M_ETA/mzp_to_mN_ratio), 30)  # for pion/eta
    #ualpha_array = np.logspace(-7, 0, 50)

    # for ebrem
    ualpha_array = np.logspace(-4, 0, 40)
    mN_array = np.logspace(np.log10(50.0), np.log10(300.0), 40)

    file = open(out_file, "w")
    file.close()
    
    # meson flux is simulated from 1e5 pot
    flux = FluxHNLFromNeutralMeson(meson_flux=pi0_flux_per_1e5POT, flux_weight=1.0e-5, coupling_BL=1e-4,
                                   mixing_flavor=flavor, n_samples=1, apply_angle_cut=True, meson_species="Pion")
    flux_eta = FluxHNLFromNeutralMeson(meson_flux=eta_flux_per_1e5POT, flux_weight=1.0e-5, coupling_BL=1e-4,
                                   mixing_flavor=flavor, n_samples=1, apply_angle_cut=True, meson_species="Eta")
    
    #flux = FluxHNLFromElectronPositron(mixing_angle=1.0, hnl_mass=mN_array[0],
     #                                  coupling_BL=1e-4, n_samples=500000, mixing_flavor=flavor, max_track_length=10.0)
    
    #flux = FluxHNLFromProtonBrem(proton_energy=120000.0, mixing_flavor=flavor,
    #                             mixing_angle=1.0, hnl_mass=mN_array[0], coupling_BL=1e-4, n_samples=100000)
    
    for mN in mN_array:
        flux.set_new_params(zprime_mass=mzp_to_mN_ratio*mN, hnl_mass=mN)
        print("------------- Setting new mass mN = {} and simulating".format(mN))
        flux.simulate()
        flux.decay_to_hnl()

        flux_eta.set_new_params(zprime_mass=mzp_to_mN_ratio*mN, hnl_mass=mN)
        print("Simulating eta...")
        flux_eta.simulate()
        flux_eta.decay_to_hnl()

        print("Sum of HNLS = {}".format(np.sum(flux.hnl_flux)*DUNE_POT_TOTAL))

        for Ualpha in ualpha_array:
            flux.propagate(Ualpha, verbose=True)
            flux_eta.propagate(Ualpha, verbose=True)

            # flux is normalized to per POT, so scale to total POT
            # 100% efficiency
            wgts = DUNE_EFFICIENCY * DUNE_POT_TOTAL * flux.decay_axion_weight
            wgts_eta = DUNE_EFFICIENCY * DUNE_POT_TOTAL * flux_eta.decay_axion_weight

            counts = np.sum(wgts) + np.sum(wgts_eta)

            print("Found {} counts at mN = {} with Ualpha={}".format(counts, mN, Ualpha))

            file = open(out_file, "a")
            file.write("{} {} {}\n".format(str(mN), str(Ualpha), str(counts)))
            file.close()





def main():
    parameter_scan(mzp_to_mN_ratio=2.1, out_file="scans/hnl_DUNE_scan_mass-ratio-2_meson_umu_v2.txt", flavor=1)
    parameter_scan(mzp_to_mN_ratio=2.1, out_file="scans/hnl_DUNE_scan_mass-ratio-2_meson_utau_v2.txt", flavor=2)
    #parameter_scan(mzp_to_mN_ratio=2.1, out_file="scans/hnl_DUNE_scan_mass-ratio-2_ebrem_ue.txt", flavor=0)

    #parameter_scan(mzp_to_mN_ratio=5, out_file="scans/hnl_DUNE_scan_mass-ratio-5_ebrem_ue_highstats.txt", flavor=0)
    parameter_scan(mzp_to_mN_ratio=5, out_file="scans/hnl_DUNE_scan_mass-ratio-5_meson_umu_v2.txt", flavor=1)
    parameter_scan(mzp_to_mN_ratio=5, out_file="scans/hnl_DUNE_scan_mass-ratio-5_meson_utau_v2.txt", flavor=2)

if __name__ == "__main__":
    main()