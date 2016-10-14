from matplotlib import pyplot as plt
import pandas as pd

def getPathName(stg_file_name, dt):
    paths = stg_file_name.split('/')
    filename = paths[-1]
    filename_parts = filename.split('_')
    filename_parts[-1] = str(dt.microsecond / 1000)+'.mat'
    filename_parts[-2] = '%02d%02d%02d' % (dt.hour, dt.minute, dt.second)
    filename_parts[-3] = '%04d%02d%02d' % (dt.year, dt.month, dt.day)
    paths[-1] = '_'.join(filename_parts)
    new_stg_file_name = '/'.join(paths)
    return new_stg_file_name

def getFileName(stg_file_name, dt):
    paths = stg_file_name.split('/')
    filename = paths[-1]
    filename_parts = filename.split('_')
    filename_parts[-1] = str(dt.microsecond / 1000)+'.mat'
    filename_parts[-2] = '%02d%02d%02d' % (dt.hour, dt.minute, dt.second)
    filename_parts[-3] = '%04d%02d%02d' % (dt.year, dt.month, dt.day)
    paths[-1] = '_'.join(filename_parts)
    new_stg_file_name = '/'.join(paths)
    return paths[-1]


# for plotting only
def rolling_std(timeseries, window=1000):
    # Determing rolling statistics
    rolmean = pd.rolling_mean(timeseries, window=window)
    rolstd = pd.rolling_std(timeseries, window=window)

    # Plot rolling statistics:
    orig = plt.plot(timeseries, color='blue', label='Original')
    mean = plt.plot(rolmean, color='black', label='Rolling Mean')
    std = plt.plot(rolstd, color='red', label='Rolling Std')
    plt.legend(loc='best')
    plt.title('Rolling Standard Deviation')
    plt.show(block=False)


def plot3C(timeseries, i_rec):
    plt.plot(timeseries['traceData'][:, i_rec * 3], color='red')
    plt.plot(timeseries['traceData'][:, i_rec * 3 + 1], color='gray')
    plt.plot(timeseries['traceData'][:, i_rec * 3 + 2], color='blue')
