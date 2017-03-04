import numpy as np

ans = np.loadtxt('answer.txt', delimiter='\n')
result = np.loadtxt('result.txt', delimiter='\n')
diff = 0.0
eps = 0.0001
for i in range(0, len(ans)):
    diff = (ans[i] - result[i]) * (ans[i] - result[i])
if diff < eps:
    print("OK")
else:
    print("NG")
