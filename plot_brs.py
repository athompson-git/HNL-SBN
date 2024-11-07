from zprime_decay_widths import *

import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.pylab import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import sys
sys.path.append("../")
from alplib.constants import *

mZp_vals = np.logspace(1, 4, 3000)

widths_nunu = decay_width_zprime_nunu(mZp_vals)
widths_ee = decay_width_zprime_dilepton(mZp_vals, ml=M_E)
widths_mumu = decay_width_zprime_dilepton(mZp_vals, ml=M_MU)
widths_tautau = decay_width_zprime_dilepton(mZp_vals, ml=M_TAU)
widths_had = decay_width_zprime_hadrons(mZp_vals)

widths_hnl = decay_width_zprime_2HNL(mZp_vals, mN=mZp_vals/2.1)
widths_hnl_5 = decay_width_zprime_2HNL(mZp_vals, mN=mZp_vals/5)


total_width_noHNL = widths_nunu + widths_tautau + widths_ee + widths_mumu + widths_had
total_width_HNL_massRatio2 = widths_nunu + widths_tautau + widths_ee + widths_mumu + widths_had + widths_hnl
total_width_HNL_massRatio5 = widths_nunu + widths_tautau + widths_ee + widths_mumu + widths_had + widths_hnl_5


# plot BRs: mass ratio 2.1
plt.plot(mZp_vals, widths_nunu/total_width_HNL_massRatio2, label=r"$Z^\prime \to \nu \nu$")
plt.plot(mZp_vals, widths_ee/total_width_HNL_massRatio2, label=r"$Z^\prime \to e^+ e^-$")
plt.plot(mZp_vals, widths_mumu/total_width_HNL_massRatio2, label=r"$Z^\prime \to \mu^+ \mu^-$")
plt.plot(mZp_vals, widths_tautau/total_width_HNL_massRatio2, label=r"$Z^\prime \to \tau^+ \tau^-$")
plt.plot(mZp_vals, widths_had/total_width_HNL_massRatio2, label=r"$Z^\prime \to $ hadrons")
plt.plot(mZp_vals, widths_hnl/total_width_HNL_massRatio2, label=r"$Z^\prime \to NN$")

plt.legend(fontsize=11, framealpha=1.0)
plt.yscale('log')
plt.xscale('log')
plt.xlim((100.0, 1e4))
plt.ylim((2e-3, 1.0))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel(r"$m_{Z^\prime}$ [MeV]", fontsize=16)
plt.ylabel(r"Branching Ratio $\Gamma(Z^\prime \to X) /  \Gamma_{Z^\prime, {\rm total}}$", fontsize=16)
plt.title(r"$m_{Z^\prime}/m_N = 2.1$", loc="right", fontsize=11)
plt.tight_layout()
plt.show()
plt.close()


# plot BRs: mass ratio 5
plt.plot(mZp_vals, widths_nunu/total_width_HNL_massRatio5, label=r"$Z^\prime \to \nu \nu$")
plt.plot(mZp_vals, widths_ee/total_width_HNL_massRatio5, label=r"$Z^\prime \to e^+ e^-$")
plt.plot(mZp_vals, widths_mumu/total_width_HNL_massRatio5, label=r"$Z^\prime \to \mu^+ \mu^-$")
plt.plot(mZp_vals, widths_tautau/total_width_HNL_massRatio5, label=r"$Z^\prime \to \tau^+ \tau^-$")
plt.plot(mZp_vals, widths_had/total_width_HNL_massRatio5, label=r"$Z^\prime \to $ hadrons")
plt.plot(mZp_vals, widths_hnl_5/total_width_HNL_massRatio5, label=r"$Z^\prime \to NN$")

plt.legend(fontsize=11, framealpha=1.0)
plt.yscale('log')
plt.xscale('log')
plt.xlim((100.0, 1e4))
plt.ylim((6e-2, 1.0))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel(r"$m_{Z^\prime}$ [MeV]", fontsize=16)
plt.ylabel(r"Branching Ratio $\Gamma(Z^\prime \to X) / \Gamma_{Z^\prime, {\rm total}}$", fontsize=16)
plt.title(r"$m_{Z^\prime}/m_N = 5$", loc="right", fontsize=11)
plt.tight_layout()
plt.show()
plt.close()


# match to 1801.04847

# plot BRs
plt.plot(1e-3*mZp_vals, widths_nunu/total_width_noHNL, label=r"$Z^\prime \to \nu \nu$")
plt.plot(1e-3*mZp_vals, widths_ee/total_width_noHNL, label=r"$Z^\prime \to e^+ e^-$")
plt.plot(1e-3*mZp_vals, widths_mumu/total_width_noHNL, label=r"$Z^\prime \to \mu^+ \mu^-$")
plt.plot(1e-3*mZp_vals, widths_tautau/total_width_noHNL, label=r"$Z^\prime \to \tau^+ \tau^-$")
plt.plot(1e-3*mZp_vals, widths_had/total_width_noHNL, label=r"$Z^\prime \to $ hadrons")

plt.legend(fontsize=11)
plt.xlim((0.0, 2.0))
plt.ylim((0.0, 1.0))
plt.xlabel(r"$m_{Z^\prime}$ [GeV]", fontsize=16)
plt.ylabel(r"Branching Ratio $\Gamma_i / \Gamma$", fontsize=16)
plt.tight_layout()
plt.show()