import numpy as np

a = np.loadtxt('delay.dat')
b = np.loadtxt('weight.dat')
c = np.loadtxt('rev_potential.dat')
print a
print a.shape
print b
print b.shape
print c
print c.shape

print np.min(a)
