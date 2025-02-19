import sys
sys.path.append("../")

from alplib.fluxes import *
from hnl_decays import *
from hnl_flux import *
from other_constants import *


def flux_generator(m_HNL, m_Zp, gBL, Ualpha, ns=10000):
    flux = FluxHNLFromElectronPositron(axion_mass=m_Zp, mixing_angle=Ualpha, hnl_mass=m_HNL, coupling_BL=gBL, n_samples=ns)

    #resonant_flux.simulate()
    print("simulating res v...")
    flux.simulate()
    flux.decay_to_hnl()
    flux.propagate(Ualpha)

    return np.array(flux.hnl_energy), np.array(flux.hnl_angle), np.array(flux.hnl_flux)


def parameter_scan_bebc(mzp_to_mN_ratio=2.1, out_file="scans/hnl_BEBC_scan_mass-ratio-2.1_brem.txt", flavor=1):

    min_pbrem_HNL_mass = 250.0/mzp_to_mN_ratio

    mN_array = np.logspace(np.log10(min_pbrem_HNL_mass), np.log10(1500), 80)
    mZp_array = mzp_to_mN_ratio * mN_array
    ualpha_array = np.logspace(-7, 0, 50)

    file = open(out_file, "w")
    file.close()
    

    #flux = FluxHNLFromNeutralMeson(meson_flux=[[0.0, 0.0, 0.0, M_PI0]], flux_weight=0.08, coupling_BL=1e-4,
    #                               mixing_flavor=flavor, n_samples=100000, apply_angle_cut=False)
    
    #flux = FluxHNLFromElectronPositron(zprime_mass=mZp_array[0], mixing_angle=1.0, hnl_mass=mN_array[0],
    #                                   coupling_BL=1e-4, n_samples=50000, mixing_flavor=flavor, max_track_length=10.0)
    
    flux = FluxHNLFromProtonBrem(proton_energy=400000.0, target=Material("Fe"), det_dist=BEBC_DIST,
                                 det_length=BEBC_LENGTH, det_area=BEBC_AREA, zprime_mass=mZp_array[0],
                                 mixing_flavor=flavor, mixing_angle=1.0, hnl_mass=mN_array[0],
                                 coupling_BL=1e-4, n_samples=10000)
    
    for mN in mN_array:
        flux.set_new_params(zprime_mass=mzp_to_mN_ratio*mN, hnl_mass=mN)
        print("------------- Setting new mass mN = {} and simulating".format(mN))
        flux.simulate()
        flux.decay_to_hnl()

        for Ualpha in ualpha_array:
            flux.propagate(Ualpha, verbose=True)

            # flux is normalized to per POT, so scale to total POT
            # 100% efficiency
            wgts = BEBC_EFF * BEBC_POT * flux.decay_axion_weight / 3 # divide by 3 since only pi mu channel is measured

            # Energy cuts and mT cut
            pT_hnl = sqrt(np.array(flux.hnl_energy)**2 - flux.hnl_mass**2) * np.sin(flux.hnl_angle)
            mt = sqrt(pT_hnl**2 + flux.hnl_mass**2) + pT_hnl
            print(mt)
            wgt_mask = (np.array(flux.hnl_energy) >= BEBC_ET_CUT) & (mt < (1869.0 - M_MU))

            counts = np.sum(wgts[wgt_mask])

            print("Found {} counts at mN = {} with Ualpha={}".format(counts, mN, Ualpha))

            file = open(out_file, "a")
            file.write("{} {} {}\n".format(str(mN), str(Ualpha), str(counts)))
            file.close()




def main():


    #parameter_scan_bebc(mzp_to_mN_ratio=2.1,
    #                    out_file="scans/recast_scans/hnl_T2K_scan_mass-ratio-2_pbrem_umu.txt",
    #                    flavor=1)
    parameter_scan_bebc(mzp_to_mN_ratio=2.1,
                        out_file="scans/recast_scans/hnl_BEBC_scan_mass-ratio-2_pbrem_utau.txt",
                        flavor=2)
    #parameter_scan_bebc(mzp_to_mN_ratio=5,
    #                    out_file="scans/recast_scans/hnl_T2K_scan_mass-ratio-5_pbrem_umu.txt",
    #                    flavor=1)
    #parameter_scan_bebc(mzp_to_mN_ratio=5,
    #                    out_file="scans/recast_scans/hnl_T2K_scan_mass-ratio-5_pbrem_utau.txt",
    #                    flavor=2)


if __name__ == "__main__":
    main()