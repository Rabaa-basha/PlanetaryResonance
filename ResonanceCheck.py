__author__ = "Samantha Lawler"
__copyright__ = "Copyright 2020"
__version__ = "1.0.0"
__maintainer__ = "Rabaa"
__email__ = "beborabaa@gmail.com"

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from decimal import Decimal

tpnum = sys.argv[1]

timetp, atp, etp, inctp, Omegatp, omegatp, Mtp, LN = np.genfromtxt(str(tpnum) + "bary.out", unpack=True)

f = open("resinfoGF1234.aei", "a")
# print>>f,"#tpnum, propa, prope, propinc, rescen, resamp"

# get proper elements, look for n:1 resonance and n:2 resonance
Ltp = (Omegatp + omegatp + Mtp) % 360.  # The Lambda for test Particles in degrees
pomegatp = (Omegatp + omegatp) % 360.  # The longitude of pericenter in degrees

# Flags
isitres = 0  # Is it resonance?
isitresamp = 999  # Its amplitude.
isitrescen = 999  # The resonance center
resname = 999  # The name 'tag'

# list of resonances to check: pp and qq for pp:qq resonance
pp = [2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9, 10]
qq = [1, 1, 2, 1, 3, 1, 2, 3, 4, 1, 1, 2, 3, 4, 1, 3, 1, 2, 4, 1]

for jj in np.arange(0, len(pp)):
    # only search near a of the MMR
    ares = 30.1 * (float(pp[jj]) / float(qq[jj])) ** (2. / 3.)
    print(tpnum, pp[jj], qq[jj], ares) # Erase before release
    # search within 2 AU of the res center
    if isitres == 0 and np.average(atp) < (ares + 2.) and np.average(atp) > (ares - 2.):
        phi = (float(pp[jj]) * Ltp - float(qq[jj]) * LN - (float(pp[jj]) - float(qq[jj])) * pomegatp) % 360

        # build a window to look for resonance, loop through integration
        window = 1000  # 10000 timesteps are output, so each window is 10% of the run length
        resyes = np.zeros(int(len(phi)/window))  # Array of 10 binary elements to check for resonance each step '10%' set to zero
        rescen = np.zeros(int(len(phi)/window))  # Array of 10 binary elements to check the res angle each step '10%' set to zero
        q = 0
        c = 0
        angles = np.arange(0, 360, 5)  # Array of angles 5 degrees increment each step

        while q + window < len(phi):  # loop through timesteps in a few big windows
            windowa = np.average(atp[q:q + window])  # Average of the semi-major axis from q -> q + 1000
            if windowa > (ares - 2.) and windowa < (ares + 2.):  # Window within 2 AUs
                windowphi = phi[q:q + window]  # Next window
                resamp = np.zeros(len(angles)) + 1
                for m in np.arange(0, len(angles) - 1):  # find out where the res angle doesn't go, proxy for resamp
                    if len(windowphi[(windowphi > angles[m]) * (windowphi < (angles[m + 1]))]) == 0:  # Why multiplication
                        resamp[m] = 0
                resyes[c] = np.average(resamp) * 180.
                rescen[c] = np.average(resamp[resamp != 0] * angles[resamp != 0])  # cannot understand angles[resamp != 0]
            else:
                resyes[c] = 180.
            q = q + window
            c = c + 1
        if len(resyes[resyes < 180.]) > 8:  # changed it, 0 instead of 8
            isitres = 1
            isitresamp = np.average(resyes)
            isitrescen = np.average(rescen)
            resname = str(pp[jj]) + ':' + str(qq[jj])
        else:
            isitresamp = 999.

    print(f, tpnum, np.average(atp), np.average(etp), np.average(inctp), resname, isitresamp, isitrescen, isitres)

# This section is for development purposes, visualising the data to assure functionality.
plt.scatter(timetp, phi, s=7)

plt.xlabel('Time / Myr')
plt.ylabel('Φ / °')

TestPAverageDist = np.average(atp)
TestPAverageDistExp = TestPAverageDist ** 1.5
NeptuneAverageDist = 30.1
NeptuneAverageDistExp = NeptuneAverageDist ** 1.5
ResonanceDecimal = TestPAverageDistExp / NeptuneAverageDistExp
ResonanceRatio = Decimal(ResonanceDecimal).as_integer_ratio()
plt.title(str(tpnum) + ' ' + str(ResonanceRatio))

plt.show()
