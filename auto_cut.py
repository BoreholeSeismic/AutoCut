import scipy.io
import numpy as np
import datetime
from copy import deepcopy

from lib import *

class TraceInfo(object):
    """
    contains traces info in a single .mat file
    """
    def __init__(self, noRec, noTrace, noWell, noSample, dt, nanoSecond):
        self.noRec = noRec
        self.noTrace = noTrace
        self.noWell = noWell
        self.noSample = noSample
        self.datetime = dt
        self.nanoSecond = nanoSecond


class AutoCut(object):
    """
    work horse
    """
    def __init__(self):
        self.traces_info = None
        self.traces = None
        self.input_dir_path = None
        self.out_dir_path = None
        self.filename = None

    def cut(self, file_name=None):
        """
        if file_name provided, only cut this file
        else cut all file in the directory
        """
        if file_name is not None:
            self.filename = file_name
            self.traces, self.traces_info = self.read_matfile(file_name)
            window = self.get_signal_window(self.traces_info, self.traces['traceData'])
            start, end = self.slice_window(window)
            self.save_matfile(start, end, self.traces_info.datetime, self.traces)
        else:  # TODO support multiple files
            raise ValueError('need a filename')
        print('Cut Done!')

    def read_matfile(self, file_name):
        traces = scipy.io.loadmat(self.input_dir_path+'/'+file_name)

        noRec = traces['headerData']['noRec'][0][0][0][0]
        noTrace = traces['headerData']['noTrace'][0][0][0][0]
        noWell = traces['headerData']['childHeaders'][0][0].shape[1]
        noSample = traces['headerData']['noSample'][0][0][0][0]
        if noSample < 500:
            raise ValueError('trace noSample too short. ')
        date_str = traces['headerData']['acquisition'][0][0][0][0][0][0]
        time_str = traces['headerData']['acquisition'][0][0][0][0][1][0]
        secondFraction_str = traces['headerData']['acquisition'][0][0][0][0][2][0]
        nanoSecond = traces['headerData']['acquisition'][0][0][0][0][3][0][0]
        dt = self.get_datetime(date_str, time_str, secondFraction_str, nanoSecond)
        traces_info = TraceInfo(noRec, noTrace, noWell, noSample, dt, nanoSecond)
        return traces, traces_info

    def get_signal_window(self, traces_info, data):
        # Core logic! Any std larger than std of silent time is a signal.
        trigger_windows = []
        # single trace trigger detection. trigger_windows[noTrace][list_pos] stores triggered positions
        for i_trace in range(traces_info.noTrace):
            rolstd = pd.rolling_std(data[:, i_trace], window=300)
        #     TODO better std_silent sampling, try median of median
            rolstd_nonan = rolstd[~np.isnan(rolstd)]
            std_silent = np.nanmean(rolstd_nonan[:400])
            trigger_windows.append(np.where(rolstd > 4 * std_silent)[0])

        # Multi-trace collaborative filtering.
        # Scan each well. If most receivers detect signal, true-positive; else false-positive.
        for i_well in range(traces_info.noWell):
            count = 0
            for j_trace in range(traces_info.noTrace / traces_info.noWell):
                no_trace = i_well * traces_info.noTrace / traces_info.noWell + j_trace
                if len(trigger_windows[no_trace]) > 1:
                    count += 1
            if count < 0.7 * traces_info.noTrace / traces_info.noWell:
                # remove the noise
                for j_trace in range(traces_info.noTrace / traces_info.noWell):
                    no_trace = i_well * traces_info.noTrace / traces_info.noWell + j_trace
                    trigger_windows[no_trace] = np.array([])
        #
        signal_window = set()
        for window in trigger_windows:
            for pos in window:
                signal_window.add(pos)
        signal_window = list(signal_window)
        signal_window.sort()
        return signal_window

    def slice_window(self, window):
        # TODO support multiple cuts in a trace batch
        info = self.traces_info
        l = 0
        r = l
        step = 30
        while r+step < len(window) and window[r+step] - window[r] < step*1.1:
            r += step

        l_cut = max(0, window[l] - 800)
        r_cut = min(info.noSample - 1, window[r] + 800)
        return l_cut, r_cut

    def save_matfile(self, start_pos, end_pos, old_date_time, old_traces):
        dt_delta = datetime.timedelta(microseconds=start_pos * 250)
        new_dt = old_date_time + dt_delta
        short_traces = deepcopy(old_traces)
        # edit metadata
        new_noSample = end_pos - start_pos
        short_traces['headerData']['noSample'][0][0][0][0] = long(new_noSample) # try to save in long type, but saved as uint16 (hacked matlab script)
        short_traces['headerData']['acquisition'][0][0][0][0][0][0] = '%02d/%02d/%04d' % (new_dt.month, new_dt.day, new_dt.year)
        short_traces['headerData']['acquisition'][0][0][0][0][1][0] = '%02d:%02d:%02d' % (new_dt.hour, new_dt.minute, new_dt.second)
        new_miliSecond = new_dt.microsecond / 1000
        short_traces['headerData']['acquisition'][0][0][0][0][2][0] = new_miliSecond
        new_nanoSecond = (new_dt.microsecond % 100000) * 1000 + self.traces_info.nanoSecond % 1000
        short_traces['headerData']['acquisition'][0][0][0][0][3][0][0] = new_nanoSecond
        # edit filename
        new_filename = getFileName(self.filename, new_dt)
        short_traces['headerData']['fileName'][0][0][0] = new_filename
        # slice data
        new_data = old_traces['traceData'][start_pos: end_pos, :]
        short_traces['traceData'] = new_data
        new_stg_file_name = getPathName(self.filename, new_dt)
        scipy.io.savemat(self.out_dir_path+'/'+new_stg_file_name, {'headerData': short_traces['headerData'], 'traceData': short_traces['traceData'], 'auxData': short_traces['auxData']}, do_compression=True, format='5')

    def set_input_dir(self, dir_path):
        self.input_dir_path = dir_path

    def set_output_dir(self, out_dir_path=None):
        """
        default output dir is the same as input dir
        """
        if out_dir_path is None:
            self.out_dir_path = self.input_dir_path
        else:
            self.out_dir_path = out_dir_path

    def get_datetime(self, date_str, time_str, secondFraction_str, nanoSecond):
        month, day, year = date_str.split('/')
        hour, minute, second = time_str.split(':')
        miliSecond = int(secondFraction_str)
        microSecond = int(nanoSecond / 1e3) % 1000
        return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), miliSecond * 1000 + microSecond)
