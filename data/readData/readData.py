import numpy as np
import matplotlib.pyplot as plt

filename_tmp = './data/991109_3_p_Bal_kj'
filename = filename_tmp + '.dat'
data = np.loadtxt(filename, skiprows=1)

time = data[:,0]
response =  data[:,1]
olfactory = data[:,2]

response = 1000 * response
max_v = np.max(response)
min_v = np.min(response)
ave_v = np.average(response)

threshold = 50

olf_dur = 0.5
bin_width = 0.25

def show_olf_strength(olf):
    ave_olf = np.average(olf)
    for i in range(0, len(olf)):
        if olf[i] > ave_olf and olf[i-1]==olf[i] and olf[i]==olf[i+1]:
                print '%f' %olf[i]
                return olf[i]
    print 'failure'
    return -1

olf_judge_strength = (show_olf_strength(olfactory) + np.average(olfactory)) * 0.5
#olf_strength = -0.0083

        
def spike_detect(res, start, end):
    spike_counter = 0
    for i in range(start, end):
        if res[i] > threshold:
            if res[i-1] < res[i] and res[i] > res[i+1]:
                spike_counter = spike_counter + 1
    return spike_counter

def get_index_of_time(time, target_time):
    for i in range(0, len(time)):
        if time[i] >= target_time:
            return i

def detect_olf_timing(olf, start):
    for i in range(start, len(olf)):
        if olf_judge_strength == -1:
            return -1
        else:
            if olf[i] > olf_judge_strength:
                return i
    return -1

def calc_spike_and_frequence(time, res, olf, start_index, spike_counter_spon_tmp):
    t_start_olf_index = detect_olf_timing(olf, start_index)
    if t_start_olf_index != -1:
        t_start_olf = time[t_start_olf_index]
        t_end_olf_index = get_index_of_time(time, t_start_olf + olf_dur)
        t_after_olf_part1 = get_index_of_time(time, t_start_olf + olf_dur + bin_width)
        t_after_olf_part2 = get_index_of_time(time, t_start_olf + olf_dur + bin_width * 2)
        t_after_olf_part3 = get_index_of_time(time, t_start_olf + olf_dur + bin_width * 3)
        t_after_olf_part4 = get_index_of_time(time, t_start_olf + olf_dur + bin_width * 4)

        if start_index == 0:
            spike_counter_spontaneous = spike_detect(res, 0, t_start_olf_index)
        else:
            spike_counter_spontaneous = spike_counter_spon_tmp
        spike_counter_part1 = spike_detect(res, t_end_olf_index, t_after_olf_part1)
        spike_counter_part2 = spike_detect(res, t_after_olf_part1, t_after_olf_part2)
        spike_counter_part3 = spike_detect(res, t_after_olf_part2, t_after_olf_part3)
        spike_counter_part4 = spike_detect(res, t_after_olf_part3, t_after_olf_part4)

        ave_f_spontaneous = spike_counter_spontaneous / (0 - time[0])
        ave_f_part1 = spike_counter_part1 / bin_width
        ave_f_part2 = spike_counter_part2 / bin_width
        ave_f_part3 = spike_counter_part3 / bin_width
        ave_f_part4 = spike_counter_part4 / bin_width
        
        #print 'firing rate: spontaneous\tsection1\tsection2\tsection3\tsection4\n'
        print '%f\t%f\t%f\t%f\t%f\n' %(ave_f_spontaneous, ave_f_part1, ave_f_part2, ave_f_part3, ave_f_part4)
        #print 'spike num: spontaneous\tsection1\tsection2\tsection3\tsection4\n'
        print '%d\t%d\t%d\t%d\t%d\n' %(spike_counter_spontaneous, spike_counter_part1, spike_counter_part2, spike_counter_part3, spike_counter_part4)
        ave_f_full = (ave_f_part1 + ave_f_part2 + ave_f_part3 + ave_f_part4) * 0.25
        print 'firing rate of full section after 1(sec) of olfactory stimulus = %f\n' % ave_f_full

        return (t_after_olf_part4, spike_counter_spontaneous)
    else:
        return (-1, spike_counter_spon_tmp)

start_index = 0
spike_counter_spon_tmp = 0
while True:
    start_index, spike_counter_spon_tmp = calc_spike_and_frequence(time, response, olfactory, start_index, spike_counter_spon_tmp)
    if start_index == -1:
        break
        
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
ax2.axhline(y=olf_judge_strength, color='black')
ax1.set_title('time-response(upper) and time-olfactory(downer)', fontsize=20)
ax1.set_ylabel('response [mV]', fontsize=18)
ax2.set_xlabel('time [sec]', fontsize=18)
ax2.set_ylabel('olfactory [V]', fontsize=18)
filename = filename_tmp + '.png'
plt.savefig(filename)
print 'max_v = %f' % max_v
print 'min_v = %f' % min_v
print 'ave_v = %f' % ave_v
#plt.show()
