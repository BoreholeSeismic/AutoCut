"""autocut.auto_cut: provides entry point main()."""

__version__ = "0.0.4"

import sys
import os
import scipy.io
import numpy as np
import datetime
from copy import deepcopy
import pandas as pd
import warnings
import argparse
import colorlog
from tqdm import tqdm
import logging
logger = colorlog.getLogger("AutoCut")
from lib import getFileName




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

class TqdmHandler(logging.StreamHandler):
    """
    progressor bar
    """
    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg)


class AutoCut(object):
    """
    work horse
    """
    def __init__(self):
        self.traces_info = None
        self.traces = None
        self.input_path = None
        self.output_path = None

    def cut(self, file_name=None):
        """
        if file_name provided, only cut this file
        else recursively cut all files in the directory
        """
        if file_name is not None:
            self.cut_single_matfile(self.input_path, file_name)
        else:
            for path, subdirs, files in os.walk(self.input_path):
                with tqdm(total=100) as pbar:
                    path_diff = path[len(self.input_path):]
                    out_path = self.output_path + path_diff
                    for i, file_name in enumerate(files):
                        pbar.update(100.0/len(files))  # update progress bar
                        self.cut_single_matfile(path, file_name, out_path)


    def cut_single_matfile(self, input_path, file_name, output_path):
        res = self.read_matfile(input_path + '/' + file_name)
        if res is None:
            return
        self.traces, self.traces_info = res
        self.traces['traceData'] = self.traces['traceData'].astype('f8')
        window = self.get_signal_window(self.traces_info, self.traces['traceData'])
        segments = self.slice_multi_window(window)
        self.save_multi_matfiles(segments, self.traces_info.datetime, self.traces, file_name, output_path)

    def read_matfile(self, file_path):
        """
        :return: if .mat file, return data; otherwise return None
        """
        filename, file_ext = os.path.splitext(file_path)
        if file_ext != '.mat':
            logger.warning('Not a mat file: %s. File skipped. ' % file_path)
            return

        logger.info('Reading mat file %s' % file_path)
        traces = scipy.io.loadmat(file_path)

        noRec = int(traces['headerData']['noRec'][0][0][0][0])
        noTrace = int(traces['headerData']['noTrace'][0][0][0][0])
        noWell = int(traces['headerData']['childHeaders'][0][0].shape[1])
        noSample = int(traces['headerData']['noSample'][0][0][0][0])
        if noSample < 500:
            raise ValueError('trace noSample too short. ')
        date_str = traces['headerData']['acquisition'][0][0][0][0][0][0]
        time_str = traces['headerData']['acquisition'][0][0][0][0][1][0]
        secondFraction_str = traces['headerData']['acquisition'][0][0][0][0][2][0]
        nanoSecond = int(traces['headerData']['acquisition'][0][0][0][0][3][0][0])
        dt = self.get_datetime(date_str, time_str, secondFraction_str, nanoSecond)
        traces_info = TraceInfo(noRec, noTrace, noWell, noSample, dt, nanoSecond)
        return traces, traces_info

    def get_signal_window(self, traces_info, data):
        logger.info('Calculating trigger points... ')
        # Core logic! Any std larger than std of silent time is a signal.
        trigger_windows = []
        # single trace trigger detection. trigger_windows[noTrace][list_pos] stores triggered positions
        for i_trace in range(traces_info.noTrace):
            rolstd = pd.rolling_std(data[:, i_trace], window=300)
            std_silent = self.get_std_silent(rolstd)
            trigger_windows.append(np.where(rolstd > 4 * std_silent)[0])

        # Multi-trace collaborative filtering.
        # Scan each well. If most receivers detect signal, true-positive; else false-positive.
        for i_well in range(traces_info.noWell):
            count = 0
            for j_trace in range(traces_info.noTrace / traces_info.noWell):
                no_trace = i_well * traces_info.noTrace / traces_info.noWell + j_trace
                if len(trigger_windows[no_trace]) > 1:
                    count += 1
            if count < 0.5 * traces_info.noTrace / traces_info.noWell:
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

    def get_std_silent(self, rolstd):
        rolstd_nonan = rolstd[~np.isnan(rolstd)]
        # (1) naive method
        # std_silent = np.mean(rolstd_nonan[:400])
        # (2) discard the top 1%. But not accurate.
        # df = pd.Series(rolstd_nonan)
        # std_silent = np.mean(df[df < df.quantile(.99)])
        # (3) Sampling. median of means
        means = []
        chunk_size = 400
        n = len(rolstd_nonan)/chunk_size
        for i in range(n):
            means.append(np.mean(rolstd_nonan[i:i+chunk_size]))
        std_silent = np.median(means)
        return std_silent

    def slice_window(self, window, l=0):
        """
        :return: right boundary of one cut
        """
        r = l
        step = 30
        while r+step < len(window) and window[r+step] - window[r] < step*1.1:
            r += step

        return r

    def slice_multi_window(self, window, l=0):
        """
        :return: list of start and end index
        """
        logger.info('Slicing into multiple event windows')
        info = self.traces_info
        step = 30
        traces_cuts = []
        while l < len(window):
            r = self.slice_window(window, l)

            if l < r:
                l_cut = max(0, window[l] - 800)
                r_cut = min(info.noSample - 1, window[r] + 800)
                traces_cuts.append((l_cut, r_cut))

            l = r + step

        logger.debug('Event Windows:'+str(traces_cuts))
        return traces_cuts

    def save_multi_matfiles(self, segments, old_date_time, old_traces, old_file_name, output_path):
        logger.info('Saving into mat files...')
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        for segment in segments:
            self.save_matfile(segment, old_date_time, old_traces, old_file_name, output_path)

    def save_matfile(self, segment, old_date_time, old_traces, old_file_name, output_path):
        start_pos, end_pos = segment
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
        new_filename = getFileName(old_file_name, new_dt)
        short_traces['headerData']['fileName'][0][0][0] = new_filename
        # slice data
        new_data = old_traces['traceData'][start_pos: end_pos, :]
        short_traces['traceData'] = new_data
        scipy.io.savemat(output_path + '/' + new_filename, {'headerData': short_traces['headerData'], 'traceData': short_traces['traceData'], 'auxData': short_traces['auxData']}, do_compression=True, format='5')

    def set_input_dir(self, dir_name=None):
        """
        default input path is './input'
        """
        if dir_name is None:
            self.set_input_path(os.getcwd() + '/input')
        else:
            self.set_input_path(os.getcwd()+'/'+dir_name)

    def set_output_dir(self, dir_name=None):
        """
        default output path is './output'
        """
        if dir_name is None:
            self.set_output_path(os.getcwd() + '/output')
        else:
            self.set_output_path(os.getcwd() + '/' + dir_name)

    def set_input_path(self, input_path):
        if input_path[0] != '/':
            raise ValueError('need absolute input path!')
        if os.path.exists(input_path) is False:
            raise ValueError('Input path does not exists!')
        else:
            self.input_path = input_path

    def set_output_path(self, output_path):
        if output_path[0] != '/':
            raise ValueError('need absolute output path!')
        if os.path.exists(output_path) is False:
            raise ValueError('Output path does not exists!')
        self.output_path = output_path

    def get_datetime(self, date_str, time_str, secondFraction_str, nanoSecond):
        month, day, year = date_str.split('/')
        hour, minute, second = time_str.split(':')
        miliSecond = int(secondFraction_str)
        microSecond = int(nanoSecond / 1e3) % 1000
        return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), miliSecond * 1000 + microSecond)




def main():
    parser = argparse.ArgumentParser(description='AutoCut utility for slicing long seismic traces into short events. Executing AutoCut version %s' % __version__)
    parser.add_argument('-i', '--input', required=True, help='input absolute path')
    parser.add_argument('-o', '--output', required=True, help='output absolute path')
    parser.add_argument('-d', '--debug', dest='debug', nargs='?', const=True, default=False, help='debug mode (default: not debug)')
    args = vars(parser.parse_args())

    if len(args) < 2:
        print("Example command: autocut input output")
        print("First argument is input path, second argument is output path. (All paths needs to be absolute path. )")
        return

    print("List of argument strings: %s" % args)
    input = args['input']
    output = args['output']
    debug = args['debug']

    assert input[0] == '/', "input path should be absolute. "
    assert output[0] == '/', "output path should be absolute. "

    logger.setLevel(logging.WARNING)
    if debug:
        logger.setLevel(logging.DEBUG)

    handler = TqdmHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(name)s | %(asctime)s | %(levelname)s | %(message)s',
        datefmt='%Y-%d-%d %H:%M:%S',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'white',
            'SUCCESS:': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white'}, ))

    logger.addHandler(handler)

    warnings.filterwarnings("ignore", category=FutureWarning)

    autocut = AutoCut()

    autocut.set_input_path(input)
    autocut.set_output_path(output)

    print("Input path %s" % autocut.input_path)
    print("Output path %s" % autocut.output_path)

    autocut.cut()
    print("AutoCut finished. ")
