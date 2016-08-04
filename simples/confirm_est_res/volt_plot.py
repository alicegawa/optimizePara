import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("./voltage_result.dat")
#data = np.loadtxt("./voltage_result.dat.crescendo")

x = data[:, 0]
y = data[:, 1]

plt.plot(x, y)
plt.savefig('graph_v.png')
plt.show()
