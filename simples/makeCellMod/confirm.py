import math as m
import sys
import matplotlib.pyplot as plt
import numpy as np

def q_inf(x, coef1, coef2, coef3):
    return 1 / (1 + np.exp(-coef1 - coef2 * np.log(x / coef3)))

if len(sys.argv) < 4:
    coef1 = 1.120
    coef2 = 2.508
    coef3 = 0.05
    print 'enter coef you like'
else:
    coef1 = sys.argv[1]
    coef2 = sys.argv[2]
    coef3 = sys.argv[3]
coef1 = float(coef1)
coef2 = float(coef2)
coef3 = float(coef3)

filename = 'graph_for_%f_%f_%f.png' % (coef1, coef2, coef3)

print coef1, coef2, coef3
x = np.linspace(1, 0.001)
y = q_inf(x, coef1, coef2, coef3)
y = y*y
plt.plot(x, y, color="k", marker="*")
plt.savefig(filename)
plt.show()
