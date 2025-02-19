import sys
sys.path.append("../")

from alplib.fluxes import *
from hnl_decays import *
from zprime_decay_widths import *

from dune_constants import *


# Grab the photon flux below 100 mrad
photon_flux_sub100mrad = np.genfromtxt("../DUNE/data/photon_flux/DUNE_target_photons_2d_sub100mrad_1e6POT.txt", delimiter=",")


pot_per_sample = 1e6
photon_flux_sub100mrad[:,2] *= 1/pot_per_sample  # converts to per POT
angle_cut = DUNE_SOLID_ANGLE
forward_photon_flux = np.array([photon_flux_sub100mrad[photon_flux_sub100mrad[:,1] <= angle_cut][:,0],
                                photon_flux_sub100mrad[photon_flux_sub100mrad[:,1] <= angle_cut][:,1],
                                photon_flux_sub100mrad[photon_flux_sub100mrad[:,1] <= angle_cut][:,2]]).transpose()

# electron/positron fluxes
pos_diff_flux_dune = np.genfromtxt("../DUNE/data/epem_flux/positron_DIFF_flux_dPhidE_20210621_TargetSim_QGSP_BIC_AllHP_POT1E6.txt")
pos_diff_flux_dune[:,1] *= 1/pot_per_sample  # per POT
el_diff_flux_dune = np.genfromtxt("../DUNE/data/epem_flux/electron_DIFF_flux_dPhidE_20210621_TargetSim_QGSP_BIC_AllHP_POT1E6.txt")
el_diff_flux_dune[:,1] *= 1/pot_per_sample  # per POT


class FluxHNLFromElectronPositron(AxionFlux):
    """
    Generator for Primakoff-produced axion flux
    Takes in the differential fluxes of electrons, positrons dN/dE in the target
    which are pointing within detector solid angle
    """
    def __init__(self, electron_flux=el_diff_flux_dune, positron_flux=pos_diff_flux_dune, target=Material("C"),
                 det_dist=DUNE_DIST, det_length=DUNE_LENGTH, det_area=DUNE_AREA, zprime_mass=0.1, coupling_BL=1e-3,
                 mixing_angle=1e-3, mixing_flavor=0, hnl_mass=0.05, n_samples=1000, max_track_length=5.0,
                 flux_interpolation="log"):
        super().__init__(zprime_mass, target, det_dist, det_length, det_area)
        self.electron_flux = electron_flux
        self.positron_flux = positron_flux
        self.hnl_mass = hnl_mass  # HNL mass
        self.ntarget_area_density = target.rad_length * AVOGADRO / (2*target.z[0])
        self.Ualpha = mixing_angle  # mixing angle
        self.gBL = coupling_BL  # B-L coupling
        self.alpha = mixing_flavor
        self.n_samples = n_samples
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []
        self.support = np.ones(n_samples)
        self.max_t = max_track_length
        self.flux_interp = flux_interpolation
    
    def set_new_params(self, zprime_mass=None, coupling_BL=None, mixing_angle=None,
                       mixing_flavor=None, hnl_mass=None):
        if hnl_mass is not None:
            self.hnl_mass = hnl_mass  # HNL mass
        if mixing_angle is not None:
            self.Ualpha = mixing_angle  # mixing angle
        if coupling_BL is not None:
            self.gBL = coupling_BL  # B-L coupling
        if mixing_flavor is not None:
            self.alpha = mixing_flavor
        if zprime_mass is not None:
            self.ma = zprime_mass
        
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []
    
    def electron_flux_dN_dE(self, energy):
        if self.flux_interp == "log":
            return np.exp(np.interp(np.log(energy), np.log(self.electron_flux[:,0]), np.log(self.electron_flux[:,1]))) \
                * np.heaviside(energy - min(self.electron_flux[:,0]), 1.0) * np.heaviside(max(self.electron_flux[:,0])-energy, 1.0)
        
        return np.interp(energy, self.electron_flux[:,0], self.electron_flux[:,1], left=0.0, right=0.0)

    def positron_flux_dN_dE(self, energy):
        if self.flux_interp == "log":
            return np.exp(np.interp(np.log(energy), np.log(self.positron_flux[:,0]), np.log(self.positron_flux[:,1]))) \
                * np.heaviside(energy - min(self.positron_flux[:,0]), 1.0) * np.heaviside(max(self.positron_flux[:,0])-energy, 1.0)
        
        return np.interp(energy, self.positron_flux[:,0], self.positron_flux[:,1], left=0.0, right=0.0)

    def positron_flux_attenuated(self, t, E0, E1):
        return self.positron_flux_dN_dE(E0) * track_length_prob(E0, E1, t)

    def electron_positron_flux_attenuated(self, t, E0, E1):
        return (self.electron_flux_dN_dE(E0) + self.positron_flux_dN_dE(E0)) * track_length_prob(E0, E1, t)
    
    def resonance_peak(self):
        return 2*pi*self.gBL**2 * ALPHA / M_E

    def simulate_res(self):
        pass

    def simulate_brem(self, electron, n_samples):
        el_energy = electron[0]
        el_wgt = electron[1]

        ea_max = el_energy * (1 - power(self.ma/el_energy, 2))
        if ea_max <= self.ma:
            return
        
        ea_rnd = power(10, np.random.uniform(np.log10(self.ma), np.log10(ea_max), n_samples))
        theta_rnd = np.zeros_like(ea_rnd)  # forward approx
        x_rnd = ea_rnd / el_energy
        mc_vol = np.log(10) * x_rnd *  (np.log10(ea_max/el_energy) - np.log10(self.ma/el_energy)) / n_samples
        diff_br = (self.ntarget_area_density * HBARC**2) * mc_vol * brem_dsigma_dx_vector(x_rnd, self.gBL, self.ma, self.target_z)

        self.axion_energy.extend(ea_rnd)
        self.axion_angle.extend(theta_rnd)
        self.axion_flux.extend(el_wgt * diff_br)
    
    def simulate(self, simulate_res=True, simulate_brem=True):
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []

        if simulate_brem:
            # simulate bremsstrahlung from the electron and positron fluxes!
            ep_min = max(self.ma, M_E)
            if ep_min > max(self.electron_flux[:,0]):
                return

            # setup electron flux grid
            epem_energy_grid = 10**np.random.uniform(log10(ep_min*1.01), log10(max(self.electron_flux[:,0])), int(sqrt(self.n_samples)))
            epem_flux_mc_vol = np.log(10) * (log10(max(self.electron_flux[:,0])) - log10(ep_min*1.01)) / int(sqrt(self.n_samples))

            for el in epem_energy_grid:
                t_depth = 10**np.random.uniform(-4, np.log10(self.max_t), 10)
                new_energy = np.random.uniform(ep_min, el, 10)
                for i in range(10):
                    flux_weight = self.electron_positron_flux_attenuated(t_depth[i], el, new_energy[i]) \
                        * np.log(10) * t_depth[i] * (np.log10(self.max_t*4)) * (el - ep_min) / 10
                    self.simulate_brem([new_energy[i], flux_weight*epem_flux_mc_vol], n_samples=int(sqrt(self.n_samples)))
        
        if simulate_res:
            # simulate resonance production and append to arrays
            resonant_energy = -M_E + self.ma**2 / (2 * M_E)
            if resonant_energy + M_E < self.ma:
                return

            if resonant_energy < M_E:
                return

            if resonant_energy > max(self.positron_flux[:,0]):
                return

            e_rnd = np.random.uniform(resonant_energy, max(self.positron_flux[:,0]), self.n_samples)
            t_rnd = np.random.uniform(0.0, self.max_t, self.n_samples)
            mc_vol = self.max_t*(max(self.positron_flux[:,0]) - resonant_energy)

            attenuated_flux = mc_vol*np.sum(self.positron_flux_attenuated(t_rnd, e_rnd, resonant_energy))/self.n_samples
            wgt = self.target_z * (self.ntarget_area_density * HBARC**2) * self.resonance_peak() * attenuated_flux

            self.axion_energy.append(self.ma**2 / (2 * M_E))
            self.axion_angle.append(self.det_sa()/2)
            self.axion_flux.append(wgt)
    
    def decay_to_hnl(self):
        if self.ma < 2*self.hnl_mass:
            return
        
        mc = Decay2Body(LorentzVector(self.ma, 0.0, 0.0, 0.0), self.hnl_mass, self.hnl_mass, n_samples=1)
        p_zprime = LorentzVector(self.ma, 0.0, 0.0, 0.0)

        phis = np.random.uniform(0.0, 2*pi, len(self.axion_flux))

        for i in range(len(self.axion_flux)):
            theta = self.axion_angle[i]
            pzp = sqrt(self.axion_energy[i]**2 - self.ma**2)
            p_zprime.set_p4(self.axion_energy[i], pzp*sin(theta)*cos(phis[i]), pzp*sin(theta)*sin(phis[i]), pzp*cos(theta))
            mc.set_new_decay(p_zprime, self.hnl_mass, self.hnl_mass)
            mc.decay()

            p4_hnl1 = mc.p1_lab_4vectors[0]
            p4_hnl2 = mc.p2_lab_4vectors[0]

            br_zprime = br_zprime_2HNL(self.ma, self.hnl_mass)

            # make sure it maps onto detector solid angle
            if p4_hnl1.theta() < self.det_sa():
                self.hnl_angle.append(p4_hnl1.theta())
                self.hnl_energy.append(p4_hnl1.energy())
                self.hnl_flux.append(self.axion_flux[i]*br_zprime)

            if p4_hnl2.theta() < self.det_sa():
                self.hnl_angle.append(p4_hnl2.theta())
                self.hnl_energy.append(p4_hnl2.energy())
                self.hnl_flux.append(self.axion_flux[i]*br_zprime)
    
    def propagate(self, Ualpha, timing_cut=None, is_isotropic=False, verbose=False):

        gamma_vis = decay_width_visible(self.hnl_mass, Ualpha, [self.alpha])
        gamma_total = total_decay_width_hnl(self.hnl_mass, Ualpha, self.alpha)

        br_visible = gamma_vis/gamma_total
        if verbose:
            print("BR visible = {}".format(br_visible))

        e_a = np.array(self.hnl_energy)
        wgt = np.array(self.hnl_flux)

        # Get axion Lorentz transformations and kinematics
        p_a = sqrt(e_a**2 - self.hnl_mass**2)
        v_a = p_a / e_a
        boost = e_a / self.hnl_mass
        tau = boost / gamma_total if gamma_total > 0.0 else np.inf * np.ones_like(boost)

        # Calculate time of flight
        if timing_cut is not None:
            tof = self.det_dist / (v_a * 1e-2 * C_LIGHT)
            delta_tof = abs((self.det_dist / (1e-2 * C_LIGHT)) - tof)
            in_timing_window_wgt = delta_tof < timing_cut
        else:
            in_timing_window_wgt = 1.0

        # Get decay and survival probabilities
        surv_prob = np.array([np.exp(-self.det_dist / METER_BY_MEV / v_a[i] / tau[i]) \
                     for i in range(len(v_a))])
        decay_prob = np.array([(1 - np.exp(-self.det_length / METER_BY_MEV / v_a[i] / tau[i])) \
                      for i in range(len(v_a))])

        self.decay_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob * decay_prob, dtype=np.float32)
        self.scatter_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob, dtype=np.float32)

        if is_isotropic:
            geom_accept = self.det_area / (4*pi*self.det_dist**2)
            self.decay_axion_weight *= geom_accept
            self.scatter_axion_weight *= geom_accept




"""
PROTON BREMSSTRAHLUNG
"""

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
                    * (z - z*(1-z)*(4*M_P**2 + mV**2)/Hfunc(z, pT, mV) \
                        + Hfunc(z, pT, mV)/(2*z*mV**2))

def sigma_NSD(s):
    # s in GeV
    # returns in mb
    return 1.76 + 19.8*power(s, 0.057)

S0_PBREM = 15.98
ETA1_PBREM = 0.45
ETA2_PBREM = 0.55
def sigma_total_proton(s):
    # s in GeV
    # returns in mb
    return 34.4 + 0.3*power(log(s/S0_PBREM), 2) + 13.1*power(s/S0_PBREM, -ETA1_PBREM) + 7.4*power(s/S0_PBREM, -ETA2_PBREM)

def dsigma(z, pT, Ep, mV, g_V):
    # returns dsigma / dz dpT in mb / MeV
    s = (2*M_P**2 + 2*Ep*M_P)*1e-6
    p_proton = sqrt(Ep**2 - M_P**2)
    sprime = 2*M_P*(p_proton*(1-z) + M_P)*1e-6
    phase_space = np.heaviside(0.1 - Hfunc(z, pT, mV) / (4*z*power((1-z)*Ep,2)), 0.0)

    # convert wV from MeV^-1 to mb
    return np.clip(2*pT*wV(z, pT, mV, g_V) * sigma_NSD(sprime) \
                                   * phase_space, a_min=0.0, a_max=np.inf)

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



class FluxHNLFromProtonBrem(AxionFlux):
    """
    Generator for Primakoff-produced axion flux
    Takes in the differential fluxes of electrons, positrons dN/dE in the target
    which are pointing within detector solid angle
    """
    def __init__(self, proton_energy=120000.0, target=Material("C"), det_dist=DUNE_DIST, det_length=DUNE_LENGTH,
                 det_area=DUNE_AREA, zprime_mass=0.1, coupling_BL=1e-3, mixing_angle=1e-3, mixing_flavor=0,
                 hnl_mass=0.05, n_samples=1000, max_track_length=5.0, flux_interpolation="log"):
        super().__init__(zprime_mass, target, det_dist, det_length, det_area)
        self.proton_energy = proton_energy
        self.p_proton = sqrt(proton_energy**2 - M_P**2)
        self.hnl_mass = hnl_mass  # HNL mass
        self.ntarget_area_density = target.rad_length * AVOGADRO / (2*target.z[0])
        self.Ualpha = mixing_angle  # mixing angle
        self.gBL = coupling_BL  # B-L coupling
        self.alpha = mixing_flavor
        self.n_samples = n_samples
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []
        self.support = np.ones(n_samples)
        self.max_t = max_track_length
        self.flux_interp = flux_interpolation
    
    def set_new_params(self, zprime_mass=None, coupling_BL=None, mixing_angle=None,
                       mixing_flavor=None, hnl_mass=None):
        if hnl_mass is not None:
            self.hnl_mass = hnl_mass  # HNL mass
        if mixing_angle is not None:
            self.Ualpha = mixing_angle  # mixing angle
        if coupling_BL is not None:
            self.gBL = coupling_BL  # B-L coupling
        if mixing_flavor is not None:
            self.alpha = mixing_flavor
        if zprime_mass is not None:
            self.ma = zprime_mass
        
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []

    def simulate(self):
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []

        z, pt, wgts = pbrem_mc_momentum_and_weights(self.p_proton, self.ma, self.gBL, n_samples=self.n_samples)
        sigma_total = sigma_total_proton((2*M_P**2 + 2*self.proton_energy*M_P)*1e-6)

        # z is energy fraction
        zp_energy = self.proton_energy * z
        zp_p = sqrt(zp_energy**2 - self.ma**2)

        sinTheta = pt / zp_p
        theta = arcsin(sinTheta)

        for i in range(len(z)):
            if wgts[i] <= 0.0:
                continue
            self.axion_angle.append(theta[i])
            self.axion_energy.append(zp_energy[i])
            self.axion_flux.append(wgts[i]/sigma_total)

    def decay_to_hnl(self):
        if self.ma < 2*self.hnl_mass:
            return
        
        mc = Decay2Body(LorentzVector(self.ma, 0.0, 0.0, 0.0), self.hnl_mass, self.hnl_mass, n_samples=1)
        p_zprime = LorentzVector(self.ma, 0.0, 0.0, 0.0)

        phis = np.random.uniform(0.0, 2*pi, len(self.axion_flux))

        for i in range(len(self.axion_flux)):
            theta = self.axion_angle[i]
            pzp = sqrt(self.axion_energy[i]**2 - self.ma**2)
            p_zprime.set_p4(self.axion_energy[i], pzp*sin(theta)*cos(phis[i]), pzp*sin(theta)*sin(phis[i]), pzp*cos(theta))
            mc.set_new_decay(p_zprime, self.hnl_mass, self.hnl_mass)
            mc.decay()

            p4_hnl1 = mc.p1_lab_4vectors[0]
            p4_hnl2 = mc.p2_lab_4vectors[0]

            br_zprime = br_zprime_2HNL(self.ma, self.hnl_mass)

            # make sure it maps onto detector solid angle
            if p4_hnl1.theta() < self.det_sa():
                self.hnl_angle.append(p4_hnl1.theta())
                self.hnl_energy.append(p4_hnl1.energy())
                self.hnl_flux.append(self.axion_flux[i]*br_zprime)

            if p4_hnl2.theta() < self.det_sa():
                self.hnl_angle.append(p4_hnl2.theta())
                self.hnl_energy.append(p4_hnl2.energy())
                self.hnl_flux.append(self.axion_flux[i]*br_zprime)
    
    def propagate(self, Ualpha, timing_cut=None, is_isotropic=False, verbose=False):

        gamma_vis = decay_width_visible(self.hnl_mass, Ualpha, [self.alpha])
        gamma_total = total_decay_width_hnl(self.hnl_mass, Ualpha, self.alpha)

        br_visible = gamma_vis/gamma_total
        if verbose:
            print("BR visible = {}".format(br_visible))

        e_a = np.array(self.hnl_energy)
        wgt = np.array(self.hnl_flux)

        # Get axion Lorentz transformations and kinematics
        p_a = sqrt(e_a**2 - self.hnl_mass**2)
        v_a = p_a / e_a
        boost = e_a / self.hnl_mass
        tau = boost / gamma_total if gamma_total > 0.0 else np.inf * np.ones_like(boost)

        # Calculate time of flight
        self.hnl_timing = self.det_dist / (v_a * 1e-2 * C_LIGHT)
        if timing_cut is not None:
            tof = self.det_dist / (v_a * 1e-2 * C_LIGHT)
            delta_tof = abs((self.det_dist / (1e-2 * C_LIGHT)) - tof)
            in_timing_window_wgt = delta_tof < timing_cut
        else:
            in_timing_window_wgt = 1.0

        # Get decay and survival probabilities
        surv_prob = np.array([np.exp(-self.det_dist / METER_BY_MEV / v_a[i] / tau[i]) \
                     for i in range(len(v_a))])
        decay_prob = np.array([(1 - np.exp(-self.det_length / METER_BY_MEV / v_a[i] / tau[i])) \
                      for i in range(len(v_a))])

        self.decay_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob * decay_prob, dtype=np.float32)
        self.scatter_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob, dtype=np.float32)

        if is_isotropic:
            geom_accept = self.det_area / (4*pi*self.det_dist**2)
            self.decay_axion_weight *= geom_accept
            self.scatter_axion_weight *= geom_accept



class FluxHNLFromNeutralMeson(AxionFlux):
    """
    Generator for Primakoff-produced axion flux
    Takes in the differential fluxes of electrons, positrons dN/dE in the target
    which are pointing within detector solid angle
    """
    def __init__(self, meson_flux=[[0.0, 0.0, 0.0, M_PI0]], flux_weight=1.0, meson_species="Pion", target=Material("C"),
                 det_dist=DUNE_DIST, det_length=DUNE_LENGTH, det_area=DUNE_AREA,
                 zprime_mass=0.1, coupling_BL=1e-3, mixing_angle=1e-3, mixing_flavor=0,
                 hnl_mass=0.05, n_samples=1000, off_axis_angle=0.0, apply_angle_cut=True):
        super().__init__(zprime_mass, target, det_dist, det_length, det_area)
        self.meson_flux = meson_flux
        self.flux_weight = flux_weight  # normalization

        meson_data = {
           "Pion": M_PI0,
            "Eta": M_ETA
        }
        self.meson_species = meson_species
        self.mass_meson = meson_data[meson_species]

        self.hnl_mass = hnl_mass  # HNL mass
        self.Ualpha = mixing_angle  # mixing angle
        self.gBL = coupling_BL  # B-L coupling
        self.alpha = mixing_flavor
        self.n_samples = n_samples
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []
        self.support = np.ones(n_samples)

        # Get detector solid angle
        self.off_axis_angle = off_axis_angle
        self.phi_range = 2*pi
        if off_axis_angle > 0.0:
            self.phi_range = self.det_sa()
        
        self.apply_angle_cut = apply_angle_cut
    
    def set_new_params(self, zprime_mass=None, coupling_BL=None, mixing_angle=None,
                       mixing_flavor=None, hnl_mass=None):
        if hnl_mass is not None:
            self.hnl_mass = hnl_mass  # HNL mass
        if mixing_angle is not None:
            self.Ualpha = mixing_angle  # mixing angle
        if coupling_BL is not None:
            self.gBL = coupling_BL  # B-L coupling
        if mixing_flavor is not None:
            self.alpha = mixing_flavor
        if zprime_mass is not None:
            self.ma = zprime_mass
        
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []
    
    def br(self):
        if self.ma > self.mass_meson:
            return 0.0
        if self.meson_species == "Eta":
            br_2gamma = 0.3941
        else:
            br_2gamma = 1.0
        return br_2gamma * 2 * (self.gBL)**2 * (1 - power(self.ma / self.mass_meson, 2))**3 / (4*pi*ALPHA)

    def simulate(self):
        self.axion_energy = []
        self.axion_angle = []
        self.axion_flux = []
        self.hnl_energy = []
        self.hnl_angle = []
        self.hnl_flux = []
        self.hnl_timing = []
        self.decay_axion_weight = []
        self.scatter_axion_weight = []

        m = self.meson_flux[0]
        pi0_p4 = LorentzVector(m[3], m[0], m[1], m[2])
        mc = Decay2Body(pi0_p4, m1=self.ma, m2=0.0, n_samples=self.n_samples)
        for m in self.meson_flux:
            pi0_p4.set_p4(m[3], m[0], m[1], m[2])
            mc.set_new_decay(p_parent=pi0_p4, m1=self.ma, m2=0.0)
            mc.decay()

            ap_energies = np.array([lv.energy() for lv in mc.p1_lab_4vectors])
            ap_thetas = np.array([lv.theta() for lv in mc.p1_lab_4vectors])
            weights = self.br() * mc.weights * self.flux_weight

            if self.apply_angle_cut:
                theta_mask = (ap_thetas < self.det_sa() + self.off_axis_angle) \
                    * (ap_thetas > self.off_axis_angle - self.det_sa())
                ap_energies = ap_energies[theta_mask]
                weights = weights[theta_mask] * self.phi_range/(2*pi)  # Assume azimuthal symmetry
                ap_thetas = ap_thetas[theta_mask]
            else:
                weights = self.det_area * weights / (4*pi*self.det_dist**2)

            self.axion_flux.extend(weights)
            self.axion_energy.extend(ap_energies)
            self.axion_angle.extend(ap_thetas)

    def decay_to_hnl(self):
        if self.ma < 2*self.hnl_mass:
            return
        
        mc = Decay2Body(LorentzVector(self.ma, 0.0, 0.0, 0.0), self.hnl_mass, self.hnl_mass, n_samples=1)
        p_zprime = LorentzVector(self.ma, 0.0, 0.0, 0.0)

        phis = np.random.uniform(0.0, 2*pi, len(self.axion_flux))

        for i in range(len(self.axion_flux)):
            theta = self.axion_angle[i]
            pzp = sqrt(self.axion_energy[i]**2 - self.ma**2)
            p_zprime.set_p4(self.axion_energy[i], pzp*sin(theta)*cos(phis[i]), pzp*sin(theta)*sin(phis[i]), pzp*cos(theta))
            mc.set_new_decay(p_zprime, self.hnl_mass, self.hnl_mass)
            mc.decay()

            p4_hnl1 = mc.p1_lab_4vectors[0]
            p4_hnl2 = mc.p2_lab_4vectors[0]

            zprime_br = br_zprime_2HNL(self.ma, self.hnl_mass)

            # make sure it maps onto detector solid angle
            if self.apply_angle_cut:
                if p4_hnl1.theta() < self.det_sa():
                    self.hnl_angle.append(p4_hnl1.theta())
                    self.hnl_energy.append(p4_hnl1.energy())
                    self.hnl_flux.append(self.axion_flux[i]*zprime_br)

                if p4_hnl2.theta() < self.det_sa():
                    self.hnl_angle.append(p4_hnl2.theta())
                    self.hnl_energy.append(p4_hnl2.energy())
                    self.hnl_flux.append(self.axion_flux[i]*zprime_br)
            else:
                self.hnl_angle.append(p4_hnl1.theta())
                self.hnl_energy.append(p4_hnl1.energy())
                self.hnl_flux.append(self.axion_flux[i]*zprime_br)
    
    def propagate(self, Ualpha, timing_cut=None, is_isotropic=False, verbose=False):


        gamma_vis = decay_width_visible(self.hnl_mass, Ualpha, [self.alpha])
        gamma_total = total_decay_width_hnl(self.hnl_mass, Ualpha, self.alpha)

        br_visible = gamma_vis/gamma_total
        if verbose:
            print("BR visible = {}".format(br_visible))

        e_a = np.array(self.hnl_energy)
        wgt = np.array(self.hnl_flux)

        # Get axion Lorentz transformations and kinematics
        p_a = sqrt(e_a**2 - self.hnl_mass**2)
        v_a = p_a / e_a
        boost = e_a / self.hnl_mass
        tau = boost / gamma_total if gamma_total > 0.0 else np.inf * np.ones_like(boost)

        # Calculate time of flight
        self.hnl_timing = self.det_dist / (v_a * 1e-2 * C_LIGHT)
        if timing_cut is not None:
            tof = self.det_dist / (v_a * 1e-2 * C_LIGHT)
            delta_tof = abs((self.det_dist / (1e-2 * C_LIGHT)) - tof)
            in_timing_window_wgt = delta_tof < timing_cut
        else:
            in_timing_window_wgt = 1.0

        # Get decay and survival probabilities
        surv_prob = np.array([np.exp(-self.det_dist / METER_BY_MEV / v_a[i] / tau[i]) \
                     for i in range(len(v_a))])
        decay_prob = np.array([(1 - np.exp(-self.det_length / METER_BY_MEV / v_a[i] / tau[i])) \
                      for i in range(len(v_a))])

        self.decay_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob * decay_prob, dtype=np.float32)
        self.scatter_axion_weight = np.asarray(br_visible * wgt * in_timing_window_wgt * surv_prob, dtype=np.float32)

        if is_isotropic:
            geom_accept = self.det_area / (4*pi*self.det_dist**2)
            self.decay_axion_weight *= geom_accept
            self.scatter_axion_weight *= geom_accept

