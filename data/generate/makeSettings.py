import numpy as np

NCELL = 86

weight = np.zeros([NCELL, NCELL])
delay = 20 * np.random.rand(NCELL,NCELL)
rev_e = 10 * np.ones([NCELL, NCELL])

print weight
print delay
print rev_e

np.savetxt('weight.dat',weight)
np.savetxt('delay.dat', delay)
np.savetxt('rev_potential.dat',rev_e)
