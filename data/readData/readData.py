import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('./data/991109_3_p_Bol_kj.dat', skiprows=1)

time = data[:,0]
response =  data[:,1]
olfactory = data[:,2]

response = 1000 * response
max_v = np.max(response)
min_v = np.min(response)
ave_v = np.average(response)

threshold = 50

def spike_detect(res, index):
    spike_counter = 0
    for i in range(0, index):
        if res[i] > threshold:
            if res[i-1] < res[i] and res[i] > res[i+1]:
                spike_counter = spike_counter + 1
    return spike_counter

t_start_olf = 0
t_start_olf_index = -1

for i in range(0, len(time)):
    if time[i] >= 0:
        t_start_olf = time[i]
        t_start_olf_index = i
        break

print 't_start_olf_index = %d' % t_start_olf_index
        
spike_counter_spontaneous = spike_detect(response, t_start_olf_index)
spike_counter_full = spike_detect(response, len(response))
print 'num of spike till 0 points = %d' % spike_counter_spontaneous
print 'num of spike of full time = %d'% spike_counter_full

fig = plt.figure(figsize=(15,6))
ax1 = fig.add_subplot(211)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax2 = fig.add_subplot(212)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax1.plot(time, response)
ax2.plot(time, olfactory)
ax1.axhline(y=max_v, color='red')
ax1.axhline(y=50, color='yellow')
ax1.axhline(y=min_v, color='green')
ax1.axhline(y=ave_v, color='black')
ax1.set_title('time-response(upper) and time-olfactory(downer)', fontsize=20)
ax1.set_ylabel('response [mV]', fontsize=18)
ax2.set_xlabel('time [sec]', fontsize=18)
ax2.set_ylabel('olfactory [V]', fontsize=18)
#plt.savefig("./data/data.png")
print 'max_v = %f' % max_v
print 'min_v = %f' % min_v
print 'ave_v = %f' % ave_v
#plt.show()
