import sys
sys.path.append("../")

from alplib.constants import *
from alplib.materials import *

# all distances with respect to BNB

# SBND
sbnd_det = Material("Ar")
SBND_POT = 6.6e20  # not certain
SBND_DIST = 110
SBND_AREA = 16.0  # 4x4
SBND_LENGTH = 5.0
SBND_NE = SBND_LENGTH * 100 * 1.39 * 6.022e23 / 40
SBND_THRESH = 100.0
# equivalent solid angle: ~ 18 mrad

# MicroBooNE
mub_det = Material("Ar")
MUB_POT = 1.36e21  # 2x previous POT = 2 * 6.80e20 = 1.36e21
MUB_DIST = 470
MUB_AREA = 2.3*2.6
MUB_LENGTH = 10.4
MUB_NE = MUB_LENGTH * 100 * 1.39 * 6.022e23 / 40
MUB_THRESH = 100.0  # MeV
# equivalent solid angle: ~ 3.5 mrad

# ICARUS
# https://cds.cern.ch/record/1640260/files/arXiv:1312.7252.pdf
icarus_det = Material("Ar")
ICARUS_POT = 6.6e20  # not certain
ICARUS_DIST = 600
ICARUS_AREA = 3.0*3.16
ICARUS_LENGTH = 17.95
ICARUS_NE = ICARUS_LENGTH * 100 * 1.39 * 6.022e23 / 40
ICARUS_THRESH = 100.0  # MeV
# equivalent solid angle: ~ 2.5 mrad

# ICARUS-NuMI [see e.g. 1909.11670]
ICARUS_NUMI_DIST = 803.0
ICARUS_NUMI_ANGLE_RAD = 0.097
ICARUS_NUMI_DUMP_ANGLE_RAD = 0.7455  # calculated approx
ICARUS_NUMI_DUMP_DIST = 114.64  # calculated approx


# fluxes for mesons
bnb_target_pi0 = np.genfromtxt("../M3DM/data/mb_target_mode/bnb_pi_zero.txt")
bnb_target_pi0[:,:] *= 1e3
bnb_target_pi0[:,3] += M_PI0