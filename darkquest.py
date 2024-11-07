import sys
sys.path.append("../")

from alplib.fluxes import *
from hnl_decays import *
from hnl_flux import *
from darkquest_constants import *



def parameter_scan_dq(mzp_to_mN_ratio=2.1, out_file="scans/hnl_DQ_scan_mass-ratio-2.1_pbrem.txt", flavor=1):

    min_pbrem_HNL_mass = 250.0/mzp_to_mN_ratio

    mN_array = np.logspace(np.log10(min_pbrem_HNL_mass), np.log10(5000), 80)
    mZp_array = mzp_to_mN_ratio * mN_array
    ualpha_array = np.logspace(-5, 0, 100)

    file = open(out_file, "w")
    file.close()
    

    #flux = FluxHNLFromNeutralMeson(meson_flux=[[0.0, 0.0, 0.0, M_PI0]], flux_weight=0.08, coupling_BL=1e-4,
    #                               mixing_flavor=flavor, n_samples=100000, apply_angle_cut=False)
    
    #flux = FluxHNLFromElectronPositron(zprime_mass=mZp_array[0], mixing_angle=1.0, hnl_mass=mN_array[0],
    #                                   coupling_BL=1e-4, n_samples=50000, mixing_flavor=flavor, max_track_length=10.0)
    
    flux = FluxHNLFromProtonBrem(proton_energy=120000.0, target=Material("Fe"), det_dist=DQ_DIST,
                                 det_length=DQ_LENGTH, det_area=DQ_AREA, zprime_mass=mZp_array[0],
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
            wgts = DQ_EFFICIENCY * DQ_POT_TOTAL * flux.decay_axion_weight

            counts = np.sum(wgts)

            print("Found {} counts at mN = {} with Ualpha={}".format(counts, mN, Ualpha))

            file = open(out_file, "a")
            file.write("{} {} {}\n".format(str(mN), str(Ualpha), str(counts)))
            file.close()





def main():

    parameter_scan_dq(mzp_to_mN_ratio=2.1,
                        out_file="scans/hnl_DQ_scan_mass-ratio-2_pbrem_umu.txt",
                        flavor=1)
    parameter_scan_dq(mzp_to_mN_ratio=2.1,
                        out_file="scans/hnl_DQ_scan_mass-ratio-2_pbrem_utau.txt",
                        flavor=2)
    parameter_scan_dq(mzp_to_mN_ratio=5,
                        out_file="scans/hnl_DQ_scan_mass-ratio-5_pbrem_umu.txt",
                        flavor=1)
    parameter_scan_dq(mzp_to_mN_ratio=5,
                        out_file="scans/hnl_DQ_scan_mass-ratio-5_pbrem_utau.txt",
                        flavor=2)


   

if __name__ == "__main__":
    main()