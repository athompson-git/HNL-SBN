import sys
sys.path.append("../")
from alplib.constants import *
from alplib.fmath import *


##################### UNIVERSAL CONSTANTS #####################
DQ_AREA = 4.0  # rough area of first tracking plane, +/- 1 meter tall, accounts for spread to longer z
DQ_THRESH = 1.0  # energy threshold [MeV]
DQ_LENGTH=14.0
DQ_DIST=4.0
DQ_SOLID_ANGLE = np.arctan(sqrt(DQ_AREA / pi) / DQ_DIST)
DQ_POT_PER_YEAR = 1.0e18
DQ_POT_TOTAL = 1.0e18
DQ_EFFICIENCY = 0.2
