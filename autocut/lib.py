# from matplotlib import pyplot as plt
# import pandas as pd
import os

def getFileName(old_file_name, dt):
    filename, file_ext = os.path.splitext(old_file_name)
    if file_ext != '.mat':
        raise ValueError('Not a mat file')

    filename_parts = filename.split('_')
    if filename_parts[-1] == 'trace3C':
        filename_parts[-2] = '%03d' % (dt.microsecond / 1000)
        filename_parts[-3] = '%02d%02d%02d' % (dt.hour, dt.minute, dt.second)
        filename_parts[-4] = '%04d%02d%02d' % (dt.year, dt.month, dt.day)

    else:
        # TODO support other file name format
        # raise ValueError('file name error: expecting name ends with trace3C')
        filename_parts[-1] = '%03d' % (dt.microsecond / 1000)
        filename_parts[-2] = '%02d%02d%02d' % (dt.hour, dt.minute, dt.second)
        filename_parts[-3] = '%04d%02d%02d' % (dt.year, dt.month, dt.day)

    filename = '_'.join(filename_parts)
    new_file_name = filename+file_ext
    return new_file_name

#
# # for plotting
# def rolling_std(timeseries, window=1000):
#     # Determing rolling statistics
#     rolmean = pd.rolling_mean(timeseries, window=window)
#     rolstd = pd.rolling_std(timeseries, window=window)
#
#     # Plot rolling statistics:
#     orig = plt.plot(timeseries, color='blue', label='Original')
#     mean = plt.plot(rolmean, color='black', label='Rolling Mean')
#     std = plt.plot(rolstd, color='red', label='Rolling Std')
#     plt.legend(loc='best')
#     plt.title('Rolling Standard Deviation')
#     plt.show(block=False)
#
# # for plotting
# def plot3C(timeseries, i_rec):
#     plt.plot(timeseries['traceData'][:, i_rec * 3], color='red')
#     plt.plot(timeseries['traceData'][:, i_rec * 3 + 1], color='gray')
#     plt.plot(timeseries['traceData'][:, i_rec * 3 + 2], color='blue')
