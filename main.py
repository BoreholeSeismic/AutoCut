from matplotlib import pyplot as plt
import pandas as pd
import scipy.io
import numpy as np


class Info(object):

    def __init__(self):
        self.file_path
        self.noRec
        self.noTrace
        self.noWell
        self.noSample
        self.

        self.datetime

class AutoCut(object):
    def __init__(self):
        month, day, year = date_str.split('/')
        hour, minute, second = time_str.split(':')
        miliSecond = int(secondFraction_str)
        microSecond = int(nanoSecond / 1e3) % 1000
        dt = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second),
                               miliSecond * 1000 + microSecond)

    def readMatFile(name = 'STG Time/Merged/Cont_Ludwig_20160818_020540_617.mat'):
        stg_file_name = 'STG Time/Merged/Cont_Ludwig_20160818_020540_617.mat'
        traces_5s = scipy.io.loadmat(stg_file_name)
        noRec = traces_5s['headerData']['noRec'][0][0][0][0]
        noTrace = traces_5s['headerData']['noTrace'][0][0][0][0]
        noWell = traces_5s['headerData']['childHeaders'][0][0].shape[1]
        noSample = traces_5s['headerData']['noSample'][0][0][0][0]
        date_str = traces_5s['headerData']['acquisition'][0][0][0][0][0][0]
        time_str = traces_5s['headerData']['acquisition'][0][0][0][0][1][0]
        secondFraction_str = traces_5s['headerData']['acquisition'][0][0][0][0][2][0]
        nanoSecond = traces_5s['headerData']['acquisition'][0][0][0][0][3][0][0]

import datetime



# Core logic! Any std larger than std of silent time is a signal.
trigger_windows = []
for i_trace in range(noTrace):
    rolstd = pd.rolling_std(traces_5s['traceData'][:,i_trace], window=300)
    rolstd_nonan = rolstd[~np.isnan(rolstd)]
#     print np.where(rolstd > 2.1*np.nanmean(rolstd_nonan[:100]))[0]
    trigger_windows.append( np.where(rolstd > 4*np.nanmean(rolstd_nonan[:400]))[0] )


# Scan each well. If most traces detect signal, cut that period.
for i_well in range(noWell):
    count = 0
    for j_trace in range(noTrace/noWell):
        no_trace = i_well * noTrace/noWell + j_trace
        if len(trigger_windows[no_trace]) > 1:
            count += 1
    if count < 0.7*noTrace/noWell:
        # remove the noise
        for j_trace in range(noTrace/noWell):
            trigger_windows[i_well * noTrace/noWell + j_trace] = np.array([])

signal_window = set()
for window in trigger_windows:
    for pos in window:
        signal_window.add(pos)
signal_window = list(signal_window)
signal_window.sort()


from copy import deepcopy
l = 0
r = l
step = 30
while r+step < len(signal_window) and signal_window[r+step] - signal_window[r] < step*1.1:
    r += step

l_cut = max(0, signal_window[l]-800)
r_cut = min(noSample-1, signal_window[r]+800)
traces_cut = traces_5s['traceData'][l_cut: r_cut, :]



tdelta = datetime.timedelta(microseconds=l_cut*250)
new_dt = dt + tdelta
new_obj = deepcopy(traces_5s)
# modify metadata
new_noSample = r_cut-l_cut
new_obj['headerData']['noSample'][0][0][0][0] = long(new_noSample) # try to save in long type, but saved as uint16 (hacked matlab script)
new_obj['headerData']['acquisition'][0][0][0][0][0][0] = '%02d/%02d/%04d'%((new_dt.month, new_dt.day, new_dt.year))
new_obj['headerData']['acquisition'][0][0][0][0][1][0] = '%02d:%02d:%02d'%((new_dt.hour, new_dt.minute, new_dt.second))
new_miliSecond = new_dt.microsecond / 1000
new_obj['headerData']['acquisition'][0][0][0][0][2][0] = new_miliSecond
new_nanoSecond = (new_dt.microsecond % 100000) * 1000 + nanoSecond%1000
new_obj['headerData']['acquisition'][0][0][0][0][3][0][0] = new_nanoSecond
new_obj['headerData']['fileName'][0][0][0] = getFileName(stg_file_name, dt)
new_obj['traceData'] = traces_cut
new_stg_file_name = getPathName(stg_file_name, dt)
# scipy.io.savemat('test/test_old.mat', {'old_obj':traces_5s}, do_compression=True)
scipy.io.savemat('test/test_new.mat', {'headerData':new_obj['headerData'], 'traceData':new_obj['traceData'], 'auxData':new_obj['auxData']}, do_compression=True, format='5')
#