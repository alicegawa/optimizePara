import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("./voltage_result.dat")

x = data[:, 0]
y = data[:, 1]

plt.plot(x, y)
plt.show()
