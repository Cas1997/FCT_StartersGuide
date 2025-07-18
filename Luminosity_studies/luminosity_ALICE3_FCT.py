import numpy as np
import matplotlib.pyplot as plt

c = 299792458 # m/s
lhcCircumference = 26659 # m
nLHCOrbitsPerS = c / lhcCircumference # Number of LHC orbits per second (Hz)
interactionRate = 24 * 1e6 # Desired interaction rate of ALICE 3 (Hz)
nDesiredCollisions = interactionRate / nLHCOrbitsPerS # number of desired collisions per orbit

nCollidingBunchCrossing = 1842 # Number of colliding bunch crossing as per Run 3
inBunchPileup = nDesiredCollisions/nCollidingBunchCrossing
print("In-bunch pileup:", inBunchPileup)

# Flat distribution - Wrong because bunches are in trains and inbetween trains there are gaps
nBunchCrossingsLHCPerOrbit = 3564 # Number of bunch crossing per LHC orbit
probColPerBunchCrossing = nCollidingBunchCrossing / nBunchCrossingsLHCPerOrbit # Probability of a collision happening per bunch crossing

# Filling scheme for 500 ns ROF
nCollidingsBunchesPerROF = np.array([23, 4, 2, 4, 2, 3, 5, 5, 15, 12, 7, 7, 23, 14, 11, 10, 6, 1, 1, 3, 22])
probNCollidingBunches = nCollidingsBunchesPerROF / np.sum(nCollidingsBunchesPerROF)

rof = 500 # readout frame duration (ns)
bunchCrossingTimeSeparation = 25 # (ns)
nBunchCrossingsPerROF = int(rof / bunchCrossingTimeSeparation)
nColPerROF = [0 for _ in range(0, 2*nBunchCrossingsPerROF)] # Factor 2 is there to store the statistically unlikely ROFs
nColPerROF_BT2T_nBunchCrossingsPerROF = 0 # Number of collisions per ROF bigger than 2 times the number of bunch crossings per ROF

nSimulatedROFs = int(1e6)
for _ in range(0, nSimulatedROFs):
	nCol = 0
	
	randFlat = np.random.rand()
	cumProb = 0
	counter = -1
	while(cumProb < randFlat):
		counter += 1
		cumProb += probNCollidingBunches[counter]
	nCollidingBunches = counter
	
	for __ in range(0, nCollidingBunches):
		nCol += np.random.poisson(inBunchPileup)

	if(nCol > 2*nBunchCrossingsPerROF - 1):
		nColPerROF_BT2T_nBunchCrossingsPerROF+=1
		continue

	nColPerROF[nCol] += 1

nColPerROFProb = np.array(nColPerROF)/nSimulatedROFs
nCollisions = []
for i in range(0, 2*nBunchCrossingsPerROF):
	nCollisions.append(i)
nCollisions = np.array(nCollisions)

# plt.figure(figsize=(10,10))
# plt.plot(nCollisions, nColPerROFProb, 'bo')
# plt.grid()
# plt.xlabel("n collisions")
# plt.ylabel("prob.")
# plt.show()

# Calculate the luminosity collected per n collisions per rof
lumIntPerYear = 3 # fb^-1
expectedNCollisionsPerROF = 0
for nCol, prob in enumerate(nColPerROFProb):
	expectedNCollisionsPerROF += prob * nCol

lumIntPerNColPerROF = lumIntPerYear * nCollisions * nColPerROFProb/expectedNCollisionsPerROF
print("Integrated luminosity:", np.sum(lumIntPerNColPerROF), "fb^-1")
n = 1
print("Integrated luminosity for n =", n, "per ROF:",lumIntPerNColPerROF[n] , "fb^-1")
print("Integrated luminosity for n =", n, "per ROF:",lumIntPerNColPerROF[n] * 1000 , "pb^-1")

# plt.figure(figsize=(11,10))
# plt.plot(nCollisions, lumIntPerNColPerROF, 'bo')
# plt.xlabel("N Collisions in ROF", fontsize=26)
# plt.ylabel("$L_{\\text{Int}} \\, (\\text{fb}^{-1})$", fontsize=26)
# plt.tick_params(axis='both', which='major', labelsize=22)
# plt.yscale("log")
# plt.grid()
# plt.show()

plt.figure(figsize=(11,10))
plt.plot(nCollisions, lumIntPerNColPerROF/lumIntPerYear, 'bo')
plt.xlabel("N Collisions in ROF", fontsize=26)
plt.ylabel("Ratio of total integrated luminosity", fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=22)
# plt.yscale("log")
plt.grid()
plt.show()
