import sys
sys.path.append("../")

from alplib.constants import *
from alplib.materials import *

# all distances with respect to BNB

# BEBC
# see e.g. https://arxiv.org/pdf/2208.00416
# volume = 3.57 × 2.52 × 1.85
BEBC_POT = 2.72e18  # not certain
BEBC_DIST = 404.0
BEBC_AREA = 3.57*2.52
BEBC_LENGTH = 1.85
BEBC_THRESH = 1000.0
BEBC_EFF = 0.96
BEBC_ET_CUT = 1000.0
# equivalent solid angle: ~ 18 mrad

# NA62

# CHARM

# T2K to ND280
# see https://arxiv.org/pdf/1902.07598
T2K_POT = 12.34e20  # not certain
T2K_DIST = 280.0
T2K_AREA = 3.4
T2K_LENGTH = 1.84
T2K_THRESH = 1000.0
T2K_EFF = 0.20
T2K_ANGLE = 0.0436

