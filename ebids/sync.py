"""Module for processing sync pulses."""

import os
import glob
import numpy as np
import pandas as pd
import neo
import scipy.interpolate as interp
import scipy.linalg as linalg
from bids import BIDSLayout

def prep_nlx_ttl(nlx_dir, bids_dir, sub, ses, task):
    """Prepare NLX TTLs for task alignment."""

    # read all TTL signals from the NLX directory
    recv_times, recv_signals = read_nlx_ttl(nlx_dir)
    data = pd.DataFrame({'onset':recv_times, 'signal':recv_signals})

    # write to a file for each run in the task (will be same for each)
    layout = BIDSLayout(bids_dir)
    runs = layout.get(subject=sub, session=ses, task=task)
    for events_file in runs:
        recv_file = events_file.path.replace('_events.tsv', '_recv.tsv')
        data.to_csv(recv_file, sep='\t', index=False)
        

def read_nlx_ttl(nlx_dir):
    """Read TTL signals from a Neuralynx data directory."""

    # neo has problems when there are multiple events files, so need
    # to hack into their data model a little
    reader = neo.io.NeuralynxIO(nlx_dir)
    event_times = []
    event_signals = []
    for chan_id in reader._nev_memmap.keys():
        data = reader._nev_memmap[chan_id]
        ids = np.unique(data['event_id'])

        # find a channel with up and down pulses
        ttl = np.zeros(len(ids), dtype=bool)
        for i, event_id in enumerate(ids):
            sigs = data['ttl_input'][data['event_id']==event_id]
            if 0 in sigs and 1 in sigs:
                ttl[i] = 1

        # check to see if there are potential matches
        if not np.any(ttl):
            continue

        # get just the ttl pulses
        isttl = np.isin(data['event_id'], ids[np.nonzero(ttl)[0]])
        event_times.append(data['timestamp'][isttl])
        event_signals.append(data['ttl_input'][isttl])

    times = np.hstack(event_times)
    signals = np.hstack(event_signals)

    return times, signals


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
