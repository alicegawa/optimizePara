import numpy as np
import matplotlib.pyplot as plt
import os
import sys

bin_width = 0.25

def make_filelist(filepath):
    def find_all_files(directory):
        for root, dirs, files in os.walk(directory):
            yield root
            for file in files:
                yield os.path.join(root, file)
    filelist = []
    for file in find_all_files(filepath):
        if '.dat' in file:
            filelist.append(file)
    return filelist

def show_olf_strength(olf):
    ave_olf = np.average(olf)
    max_olf = np.max(olf)
    for i in range(0, len(olf)):
        #if olf[i] > ave_olf and olf[i-1]==olf[i] and olf[i]==olf[i+1]:
        if olf[i] == max_olf:
            #print '%f' %olf[i+600]
            print '%f' %olf[i]
            #return olf[i+600]
            return olf[i]
    print 'failure'
    return -1

def get_spike_threshold(res, start, end):
    #print 'start = %d, end = %d, ref(lenght of res) = %d\n' %(start, end, len(res))
    max_v_in_range = np.max(res[start:(end + 1)])
    min_v_in_range = np.min(res[start:(end + 1)])
    threshold_tmp_in_range = (max_v_in_range + min_v_in_range) * 0.5
    threshold_tmp_whole = (np.max(res) + np.min(res)) * 0.5
    # voltage_range = np.max(res) - np.min(res)
    # if voltage_range < 15:
    #     return np.max(res) + 1
    if threshold_tmp_in_range > threshold_tmp_whole:
        return threshold_tmp_in_range
    else:
        return threshold_tmp_whole

def spike_detect(res, start, end):
    if start == -1 or end == -1:
        return -1
    spike_counter = 0
    threshold = get_spike_threshold(res, start, end)
    for i in range(start, end):
        if res[i] > threshold:
            if res[i-1] < res[i] and res[i] > res[i+1]:
                spike_counter = spike_counter + 1
    return spike_counter

def get_index_of_time(time, target_time):
    for i in range(0, len(time)):
        if time[i] >= target_time:
            return i
    print 'too short experiment time'
    return -1

def detect_olf_timing(olf, start, olf_judge_strength):
    def detect_olf_start_timing(olf, start, olf_judge_strength):
        for i in range(start, len(olf)):
            if olf_judge_strength == -1:
                return -1
            else:
                if olf[i] > olf_judge_strength:
                    return i
        return -1

    def detect_olf_end_timing(olf, start, olf_judge_strength):
        for i in range(start, len(olf)):
            if olf_judge_strength == -1:
                return -1
            else:
                if olf[i] < olf_judge_strength:
                    return i
        return -1
    
    start_index = detect_olf_start_timing(olf, start, olf_judge_strength)
    end_index = detect_olf_end_timing(olf, start_index, olf_judge_strength)
    return (start_index, end_index)

def calc_spike_and_frequence(time, res, olf, olf_judge_strength, start_index, num_of_call, spike_counter_spon_tmp, spon_time_tmp, filename_tmp):
    t_start_olf_index, t_end_olf_index = detect_olf_timing(olf, start_index, olf_judge_strength)
    #print 't_end_olf_index = %d, olf_judge_strength = %f\n' % (t_end_olf_index, olf_judge_strength)
    print 'start_olf_index = %d, end_olf_inex = %d\n' %(t_start_olf_index, t_end_olf_index)
    print 'start_olf_time = %f, end_olf_time = %f\n' %(time[t_start_olf_index], time[t_end_olf_index])
    if t_start_olf_index != -1:
        f1 = open("freqs.csv", "a")
        f2 = open("spikenum.csv", "a")
        t_start_olf = time[t_start_olf_index]
        if t_end_olf_index != -1:
            continue_loop_flag = 1
            t_end_olf = time[t_end_olf_index]
            # print 't_end_olf_index = %d, time is %f' % (t_end_olf_index, t_start_olf)
            # t_end_olf_index_ref = get_index_of_time(time, t_start_olf + 0.5)
            # print 'reference = %d, time is %f' % (t_end_olf_index_ref, t_end_olf)
            
        else:
            continue_loop_flag = 0
            t_end_olf_index = get_index_of_time(time, t_start_olf + bin_width)
            t_end_olf = time[t_end_olf_index]

        t_between_olf = get_index_of_time(time, t_end_olf - bin_width)
        t_after_olf_part1 = get_index_of_time(time, t_end_olf + bin_width)
        t_after_olf_part2 = get_index_of_time(time, t_end_olf + bin_width * 2)
        t_after_olf_part3 = get_index_of_time(time, t_end_olf + bin_width * 3)
        t_after_olf_part4 = get_index_of_time(time, t_end_olf + bin_width * 4)

        if start_index == 0:
            spike_counter_spontaneous = spike_detect(res, 0, t_start_olf_index)
            spon_time = t_start_olf - time[0]
        else:
            spike_counter_spontaneous = spike_counter_spon_tmp
            spon_time = spon_time_tmp

        spike_counter_between_stim = spike_detect(res, t_between_olf, t_end_olf_index)
        spike_counter_part1 = spike_detect(res, t_end_olf_index, t_after_olf_part1)
        spike_counter_part2 = spike_detect(res, t_after_olf_part1, t_after_olf_part2)
        spike_counter_part3 = spike_detect(res, t_after_olf_part2, t_after_olf_part3)
        spike_counter_part4 = spike_detect(res, t_after_olf_part3, t_after_olf_part4)

        ave_f_between_stim = spike_counter_between_stim / bin_width
        ave_f_part1 = spike_counter_part1 / bin_width
        ave_f_part2 = spike_counter_part2 / bin_width
        ave_f_part3 = spike_counter_part3 / bin_width
        ave_f_part4 = spike_counter_part4 / bin_width
        ave_f_spontaneous = spike_counter_spontaneous / spon_time
        
        #print 'firing rate: spontaneous\tsection1\tsection2\tsection3\tsection4\n'
        print '%f\t%f\t%f\t%f\t%f\t%f\n' %(ave_f_spontaneous, ave_f_between_stim, ave_f_part1, ave_f_part2, ave_f_part3, ave_f_part4)
        #print 'spike num: spontaneous\tsection1\tsection2\tsection3\tsection4\n'
        print '%d\t%d\t%d\t%d\t%d\t%d\n' %(spike_counter_spontaneous, spike_counter_between_stim, spike_counter_part1, spike_counter_part2, spike_counter_part3, spike_counter_part4)
        ave_f_full = (ave_f_between_stim + ave_f_part1 + ave_f_part2 + ave_f_part3 + ave_f_part4) * 0.2
        print 'firing rate of full section after 1(sec) of olfactory stimulus = %f\n' % ave_f_full

        f1.write('%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n' %(filename_tmp, num_of_call, ave_f_spontaneous, ave_f_between_stim, ave_f_part1, ave_f_part2, ave_f_part3, ave_f_part4))
        f2.write('%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n' %(filename_tmp, num_of_call, spon_time, spike_counter_spontaneous, spike_counter_between_stim, spike_counter_part1, spike_counter_part2, spike_counter_part3, spike_counter_part4))
        f1.close()
        f2.close()
        if continue_loop_flag == 0:
            t_after_olf_part4 = -1
        return (t_after_olf_part4, spike_counter_spontaneous, spon_time)
    else:
        return (-1, spike_counter_spon_tmp, spon_time_tmp)

def main():
    argvs = sys.argv
    argc = len(argvs)
    if argc != 2:
        print 'Usage: #python %s <directory name>' % argvs[0]
        quit()
    filename = argvs[1]
    filelist = make_filelist(filename)
    print filelist
    filelist_tmp = []
    for i in range(0, len(filelist)):
        filelist_tmp.append(filelist[i].replace('.dat',''))
    print filelist_tmp
    for i in range(0, len(filelist)):
        print 'start %s' % filelist[i]
        data = np.loadtxt(filelist[i], skiprows=1)
        if data.shape[1] != 3:
            print 'incorrect data size'
            continue
        time = data[:,0]
        response =  data[:,1]
        olfactory = data[:,2]
        if response[0] < 1:
            response = 1000 * response
        #for graph drawing
        max_v = np.max(response)
        min_v = np.min(response)
        ave_v = np.average(response)

        olf_judge_strength = (show_olf_strength(olfactory) + np.average(olfactory)) * 0.5
        #olf_judge_strength = np.average(olfactory)
 
        start_index = 0
        spike_counter_spon_tmp = 0
        spon_time_tmp = 1
        counter = 0
        while True:
            counter += 1
            start_index, spike_counter_spon_tmp, spon_time_tmp = calc_spike_and_frequence(time, response, olfactory, olf_judge_strength, start_index, counter, spike_counter_spon_tmp, spon_time_tmp, filelist_tmp[i])
            if start_index == -1:
                break
        fig = plt.figure(figsize=(22,9))
        ax1 = fig.add_subplot(211)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        ax2 = fig.add_subplot(212)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        ax1.plot(time, response)
        ax2.plot(time, olfactory)
        ax1.axhline(y=max_v, color='red')
        ax1.axhline(y=min_v, color='green')
        ax1.axhline(y=ave_v, color='black')
        ax2.axhline(y=olf_judge_strength, color='black')
        ax1.set_title('time-response(upper) and time-olfactory(downer)', fontsize=20)
        ax1.set_ylabel('response [mV]', fontsize=18)
        ax2.set_xlabel('time [sec]', fontsize=18)
        ax2.set_ylabel('olfactory [V]', fontsize=18)
        graphname = filelist_tmp[i] + '.png'
        plt.savefig(graphname)
        plt.close()
        print 'max_v = %f' % max_v
        print 'min_v = %f' % min_v
        print 'ave_v = %f' % ave_v

    
if __name__ == '__main__':
    main()
