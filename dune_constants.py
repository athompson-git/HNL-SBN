import sys
sys.path.append("../")
from alplib.constants import *
from alplib.fmath import *


##################### UNIVERSAL CONSTANTS #####################
DUNE_MASS = 50000  # kg
DUNE_ATOMIC_MASS = 37.211e3  # mass of target atom in MeV
DUNE_NTARGETS = DUNE_MASS * MEV_PER_KG / DUNE_ATOMIC_MASS
DUNE_ATOMS_PER_VOLUME = AVOGADRO * 1.3982 / 40
LAR_Z = 18  # atomic number
EXPOSURE_YEARS = 7.0
EXPOSURE_DAYS = EXPOSURE_YEARS*365  # days of exposure
DUNE_AREA = 7.0*3.0  # cross-sectional det area
DUNE_THRESH = 1.0  # energy threshold [MeV]
DUNE_LENGTH=5.0
DUNE_DIST=574
DUNE_SOLID_ANGLE = np.arctan(sqrt(DUNE_AREA / pi) / DUNE_DIST)
DUNE_POT_PER_YEAR = 1.1e21
DUNE_POT_TOTAL = 1.47e22
DUNE_EFFICIENCY = 0.2


# Colors for plottting
COLOR_NUE = 'silver'
COLOR_NUMU = 'tan'
COLOR_NUEBAR = 'teal'
COLOR_NUMUBAR = 'rosybrown'


FLUX_SCALING = EXPOSURE_DAYS*S_PER_DAY * (DUNE_POT_PER_YEAR / 365.0 / S_PER_DAY)  # gets (exposure in s) * (# / s)
E2GAMMA_MISID = 0.18

NU_MU_REWGT = (EXPOSURE_YEARS/3.5) * (2.59e8/1e7)  # based on 3.5 year reweighting
NU_MUBAR_REWGT = (EXPOSURE_YEARS/3.5) * (1.2e7/1e7)
NU_E_REWGT = (EXPOSURE_YEARS/3.5) * (3.59e6/1e7)
NU_EBAR_REWGT = (EXPOSURE_YEARS/3.5) * (4.65e5/1e7)
