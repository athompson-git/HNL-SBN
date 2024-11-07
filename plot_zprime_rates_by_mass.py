from hnl_flux import *

import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

petite_res_dune_rates = np.genfromtxt("resonance_kk_repro/petite_resonance_dune_Nv_vs_mV.txt")
petite_brem_dune_rates = np.genfromtxt("resonance_kk_repro/petite_brem_dune_Nv_vs_mV.txt")
petite_comp_dune_rates = np.genfromtxt("resonance_kk_repro/petite_compton_dune_Nv_vs_mV.txt")

def plot_flux_by_mass():
    mass_list = np.logspace(1, 4, 200)

    pbrem_flux_list = []
    res_vec_flux_list = []
    brem_v_flux_list = []

    brem_flux_epem = FluxHNLFromElectronPositron(zprime_mass=mass_list[0], coupling_BL=1.0, n_samples=500000)
    res_flux_epem = FluxHNLFromElectronPositron(zprime_mass=mass_list[0], coupling_BL=1.0, n_samples=500000)
    brem_flux_proton = FluxHNLFromProtonBrem(zprime_mass=mass_list[0], coupling_BL=1.0, n_samples=200000)

    for ma in mass_list:
        print("simulating res for ma = ", ma)
        
        # vector fluxes
        brem_flux_epem.set_new_params(zprime_mass=ma)
        res_flux_epem.set_new_params(zprime_mass=ma)
        brem_flux_proton.set_new_params(zprime_mass=ma)

        #resonant_flux.simulate()
        print("simulating resonance + ebrem...")
        brem_flux_epem.simulate(simulate_res=False)
        res_flux_epem.simulate(simulate_brem=False)

        print("simulating brem...")
        brem_flux_proton.simulate()

        brem_v_flux_list.append(1.47e22 * np.sum(brem_flux_epem.axion_flux))
        res_vec_flux_list.append(1.47e22 * np.sum(res_flux_epem.axion_flux))
        pbrem_flux_list.append(1.47e22 * np.sum(brem_flux_proton.axion_flux))
    
    # PLOT OUR RATES
    plt.plot(mass_list, pbrem_flux_list, label=r"$p p \to p p Z^\prime$", color='r', linewidth=1.0)
    plt.plot(mass_list, brem_v_flux_list, label=r"$e^\pm A \to e^\pm A Z^\prime$", color="b", ls='dashed', linewidth=1.0)
    plt.plot(mass_list, res_vec_flux_list, label=r"$e^+ e^- \to Z^\prime$", color="b", linewidth=1.0)
    
    #plt.plot(petite_res_dune_rates[:,0], petite_res_dune_rates[:,1], label=r"PETITE $e^+ e^- \to V$", color='r', ls='dashed')
    #plt.plot(petite_brem_dune_rates[:,0], petite_brem_dune_rates[:,1], label="PETITE Vector brem", ls='dashed', color="b")
    #plt.plot(petite_comp_dune_rates[:,0], petite_comp_dune_rates[:,1], label="PETITE Compton brem", ls='dashed')

    plt.ylabel(r"$N_{Z^\prime}$ / ($1.47 \cdot 10^{22}$ POT) / $g^2$", fontsize=16)
    plt.xlabel(r"$m_{Z^\prime}$ [MeV]", fontsize=16)
    plt.yscale('log')
    plt.xscale('log')
    #plt.ylim((1e14, 1e25))
    plt.xlim((10,1e4))
    plt.title(r"DUNE-ND", loc="right", fontsize=12)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()



def main():
    plot_flux_by_mass()


if __name__ == "__main__":
    main()