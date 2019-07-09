"""Module for processing sync pulses."""

import os
import glob
import numpy as np
import neo
import scipy.interpolate as interp
import scipy.linalg as linalg

def convert_nlx_ttl(nlx_dir, bids_dir, sub, ses, task):
    """Convert NLX TTLs to standard format."""

    event_files = glob(os.path.join(nlx_dir, '*.nev'))
    
    

def read_nlx_ttl(nlx_dir):
    """Read TTL signals from a Neuralynx data directory."""

    reader = neo.io.NeuralynxIO(nlx_dir)
    ttl_signals = {}
    for i in range(reader.event_channels_count()):
        times, durations, labels = reader.get_event_timestamps(
            event_channel_index=i)

        if 'ttl' in labels[0].lower():
            signal = labels[0].split('(')[1][:-2]
            ttl_signals[signal] = times

    if '0x0000' not in ttl_signals:
        raise IOError('No zero signals found in {}'.format(reader.dirname))
    if '0x0001' not in ttl_signals:
        raise IOError('No one signals found in {}'.format(reader.dirname))

    return ttl_signals


def binary2analog(signal_on, signal_off, interval):
    """Generate a signal from binary event times."""

    # combine on and off times
    signal_times = np.hstack((signal_off.flatten(),
                              signal_on.flatten()))
    signal_value = np.hstack((np.zeros(len(signal_off)),
                              np.ones(len(signal_on))))
    ind = np.argsort(signal_times)

    # sort by total time
    sort_times = signal_times[ind]
    sort_value = signal_value[ind]

    # initialize the signal
    min_time = np.min(signal_times)
    max_time = np.max(signal_times)
    n_samp = int(np.ceil((max_time - min_time) / interval))
    signal = np.zeros(n_samp)
    times = np.linspace(min_time, max_time, n_samp)

    # set each sample based on the most recent event
    for i in range(n_samp):
        last_ind = np.nonzero(times[i] >= sort_times)[0][-1]
        signal[i] = sort_value[last_ind]

    return times, signal


def align(send, recv):

    ns = len(send)
    nr = len(recv)
    nshift = nr - ns + 1
    m = np.zeros(nshift)
    for i in range(nshift):
        m[i] = np.corrcoef(send, recv[i:(i+ns)])[0,1]
    return m


def align_series(par, send, send_time, recv, recv_time):
    """Calculate error in aligning a timeseries with a shorter timeseries."""
    
    # calculate the shifted time for the short series
    time_hat = par[0] + par[1] * send_time

    # interpolate to get the predicted signal
    f = interp.interp1d(recv_time, recv)
    recv_hat = f(time_hat)
    import pdb; pdb.set_trace()
    err = linalg.norm(recv_hat - recv, 2)
    
    return err


def align_ttl(ttl_on, ttl_off, sync_on, sync_off):

    ttl_min = np.min(np.hstack((ttl_on, ttl_off)))
    ttl_max = np.max(np.hstack((ttl_on, ttl_off)))
    #ttl_signal = np.
