"""Module for processing sync pulses."""

import os
import glob
import numpy as np
import pandas as pd
import neo
import scipy.interpolate as interp
import scipy.linalg as linalg
import scipy.optimize as optim
from bids import BIDSLayout

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
        

def binary2analog(event_times, event_signal, interval):
    """Generate a signal from binary event times."""

    # initialize the signal
    min_time = np.min(event_times)
    max_time = np.max(event_times)
    n_samp = int(np.ceil((max_time - min_time) / interval))

    # set each sample based on the most recent event
    times = np.linspace(min_time, max_time, n_samp)
    f = interp.interp1d(event_times, event_signal, 'previous')
    signal = f(times)
    return times, signal


def find_signal_blocks(times, values, signal, blink, minwait, maxwait, tol):
    """Find blocks where a specified signal was sent."""

    # find blinks (periods of a given signal of a given duration)
    d_time = np.hstack((np.diff(times), 0))
    blink_start = ((d_time > blink - tol) &
                   (d_time < blink + tol) & (values == signal))
    blink_times = times[blink_start]

    # find blinks that are followed by another blink within the
    # correct time range
    d_blink = np.hstack((np.diff(blink_times), 0))
    inseq = (d_blink > minwait-.1) & (d_blink < maxwait+.1)

    # start and finish of each block of blinks
    breaks = np.nonzero(inseq==False)[0]
    start = blink_times[np.hstack((0, breaks[:-1] + 1))]
    finish = blink_times[breaks]

    # account for blink duration at the end of each block
    adjust = np.array([times[np.nonzero(times==x)[0][0]+1] for x in finish])
    return start, adjust


def align1(offset, send, recv):
    return align([offset, 1], send, recv)


def align(par, send, recv):
    """Calculate error in aligning a timeseries with a shorter timeseries."""
    
    # calculate the shifted time for the short series
    time_hat = par[0] + par[1] * send['times']
    if (np.any(time_hat > recv['times'][-1]) or
        np.any(time_hat < recv['times'][0])):
        err = 1e9
        return err

    # get bounds for the window
    start = np.argmin(np.abs(recv['times'] - time_hat[0]))
    finish = start + len(time_hat)
    if finish > len(recv['signal']):
        err = 1e9
        return err

    # interpolate to get the predicted signal
    time_recv = recv['times'][start:finish]
    time_interp = np.hstack((time_recv[0], time_hat, time_recv[-1]))
    signal_interp = np.hstack((send['signal'][0], send['signal'],
                               send['signal'][-1]))
    f = interp.interp1d(time_interp, signal_interp)
    recv_hat = f(time_recv)

    # calculate error
    err = linalg.norm(recv_hat - recv['signal'][start:finish], 2)
    return err


def load_sync_signal(sync_file, interval=0.01, scale=1):
    """Get events and samples for send and receive signals."""

    # load file
    if not os.path.exists(sync_file):
        warnings.warn('Sync file not found: {}'.format(sync_file),
                      RuntimeWarning)
    sync = pd.read_csv(sync_file, delimiter='\t')

    # translate events into continous signals
    sig_sync_times, sig_sync = binary2analog(sync.onset.values * scale,
                                             sync.signal.values, interval)
    d_sync = {'times':sig_sync_times, 'signal':sig_sync,
              'event_times':sync.onset.values,
              'event_signal':sync.signal.values}
    return d_sync
