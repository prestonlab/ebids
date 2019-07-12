"""Module for processing sync pulses."""

import os
import glob
import numpy as np
import pandas as pd
import json
import neo
import scipy.interpolate as interp
import scipy.linalg as linalg
import scipy.optimize as optim
import matplotlib.pyplot as plt
from bids import BIDSLayout

def binary2analog(event_times, event_signal, interval):
    """Generate a signal from binary event times."""

    # initialize the signal
    min_time = np.min(event_times)
    max_time = np.max(event_times)
    n_samp = int(np.ceil((max_time - min_time) / interval))
    if n_samp > 10e8:
        raise ValueError('Too many samples to hold in memory.')

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
    """Get events and samples for a set of sync signals."""

    # load file
    if not os.path.exists(sync_file):
        warnings.warn('Sync file not found: {}'.format(sync_file),
                      RuntimeWarning)
    sync = pd.read_csv(sync_file, delimiter='\t')

    # translate events into continous signals
    sig_sync_times, sig_sync = binary2analog(sync.onset.values * scale,
                                             sync.signal.values, interval)
    d_sync = {'times':sig_sync_times, 'signal':sig_sync,
              'event_times':sync.onset.values * scale,
              'event_signal':sync.signal.values}
    return d_sync


def align_run_reg(send_file, recv_file, sync_file,
                  send_scale=1, recv_scale=1, spacing=0.05):
    """Align sync pulses for a run using regression."""

    # load send and receive pulses
    send = load_sync_signal(send_file, scale=send_scale)
    recv = load_sync_signal(recv_file, scale=recv_scale)

    # initial brute-force search for the correct offset
    print(f'Aligning sync pulses: {send_file}')
    s_range = send['times'][-1] - send['times'][0]
    ranges = (slice(recv['times'][0], recv['times'][-1]-s_range, spacing),)
    x0 = optim.brute(align1, ranges, (send, recv))

    # refining search that allows slope to change
    x = optim.fmin(align, (x0, 1), (send, recv), disp=False)

    # display a summary of the alignment
    sse = align(x, send, recv)
    mse = sse / len(send['times'])
    sync_start = recv['event_times'][0]
    start = x[0] + x[1] * send['event_times'][0]
    finish = x[0] + x[1] * send['event_times'][-1]
    print(f'Offset: {x[0]:.0f}')
    print(f'Scale:  {x[1]:.8f}')
    print(f'MSE:    {mse:.8f}')
    print(f'Start:  {start-sync_start:.0f} s')
    print(f'Stop:   {finish-sync_start:.0f} s')

    # write to a sync file
    d = {'offset':x[0], 'slope':x[1], 'scale':1/recv_scale}
    with open(sync_file, 'w') as f:
        json.dump(d, f)


def align_session(bids_dir, sub, ses, rec, send_scale=1, recv_scale=1,
                  spacing=0.05):
    """Align sync pulses for all runs in a session."""

    layout = BIDSLayout(bids_dir)
    runs = layout.get(subject=sub, session=ses)
    for events in runs:
        # get send and receive files
        events_file = events.path
        send_file = events_file.replace('_events.tsv', '_send.tsv')
        recv_file = events_file.replace('_events.tsv', '_recv.tsv')

        # save sync information in a separate file
        sync_file = events_file.replace('_events.tsv', f'_{rec}.json')
        align_run_reg(send_file, recv_file, sync_file,
                      send_scale, recv_scale, spacing)


def sync_events(events_file, sync_file, rec):
    """Add recording times to events."""

    # load events
    events = pd.read_csv(events_file, delimiter='\t')

    with open(sync_file, 'r') as f:
        sync = json.load(f)

    # get equivalent times in the recording
    rec_times = (sync['offset'] + sync['slope'] * events.onset) * sync['scale']
    events[rec] = rec_times.astype(int)

    # write events back out with the new field
    events.to_csv(events_file, sep='\t', float_format='%.3f',
                  index=False, na_rep='n/a')


def sync_session(bids_dir, sub, ses, rec):
    """Add recording times to all events in a session."""

    layout = BIDSLayout(bids_dir)
    runs = layout.get(subject=sub, session=ses)
    for events in runs:
        events_file = events.path
        sync_file = events_file.replace('_events.tsv', f'_{rec}.json')
        sync_events(events_file, sync_file, rec)


def plot_sync_session(bids_dir, sub, ses, out_file=None, scale=1):
    """Plot all sync pulses and events for a session."""

    layout = BIDSLayout(bids_dir)
    run_sync = []
    run_signal = []
    run_event = []
    runs = layout.get(subject=sub, session=ses)
    for events in runs:
        events_file = events.path
        recv_file = events_file.replace('_events.tsv', '_recv.tsv')
        recv = load_sync_signal(recv_file, scale=scale)
        events = pd.read_csv(events_file, delimiter='\t')
        
        run_sync.append(recv['event_times'])
        run_signal.append(recv['event_signal'])
        run_event.append(events.ieeg)

    sync = np.hstack(run_sync)
    signal = np.hstack(run_signal)
    event = np.hstack(run_event)
    sig_sync_times, sig_sync = binary2analog(sync, signal, 0.01)

    fig, ax = plt.subplots(figsize=(20,4), dpi=300)
    ax.plot(sig_sync_times, sig_sync, linewidth=0.1)
    ax.vlines(event * scale, 0, 1, colors='r', linewidth=0.5)
    plt.tight_layout()

    # print to file
    if out_file is None:
        out_file = os.path.join(runs[0].dirname,
                                f'sub-{sub}_ses-{ses}_sync.pdf')
    fig.savefig(out_file)
