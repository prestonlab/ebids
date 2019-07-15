"""Module for processing Neuralynx data."""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import neo
import neo.rawio.neuralynxrawio as nlxio
from bids import BIDSLayout
import sync


def timediff(t1, t2):
    """Difference in seconds between datetimes."""

    if t1 > t2:
        delta = t1 - t2
    else:
        delta = t2 - t1
    return delta.seconds


def read_rec_info(nlx_dir):
    """Read timing information for recordings in a directory."""

    HEADER_SIZE = 2 ** 14  # file have a txt header of 16kB
    extensions = ['nse', 'ncs', 'nev', 'ntt']

    # read all file headers to get opening and closing times
    rec = {}
    d_header = {}
    for filename in sorted(os.listdir(nlx_dir)):
        filename = os.path.join(nlx_dir, filename)

        _, ext = os.path.splitext(filename)
        ext = ext[1:]  # remove dot
        if ext not in extensions:
            continue

        if (os.path.getsize(filename) <= HEADER_SIZE) and (ext in ['ncs']):
            continue

        # All file have more or less the same header structure
        try:
            info = nlxio.read_txt_header(filename)
        except ValueError:
            print('Cannot open file: {}'.format(filename))
            continue

        recid = (info['recording_opened'], info['recording_closed'])

        if recid not in rec:
            rec[recid] = []
        rec[recid].append(filename)
        d_header[filename] = info

    # map together files that were close to each other
    tol = 10
    dts = []
    dt_dict = {}
    for dtpair in rec.keys():
        for dt in dtpair:
            if (dt not in dts) and dts:
                # find the closest date
                delta = [timediff(dt, x) for x in dts]
                ind = np.argmin(delta)

                if delta[ind] < tol:
                    # close enough
                    dt_dict[dt] = dts[ind]
                else:
                    # make a new cluster
                    dts.append(dt)
                    dt_dict[dt] = dt
            else:
                dts.append(dt)
                dt_dict[dt] = dt

    # re-sort files based on the binned times
    filt = {}
    for dtpair in rec.keys():
        newpair = (dt_dict[dtpair[0]], dt_dict[dtpair[1]])
        if newpair not in filt:
            filt[newpair] = []

        for filename in rec[dtpair]:
            filt[newpair].append(filename)

    # output files with start and stop times
    sortpairs = sorted(set(filt.keys()))
    start = [x[0] for x in sortpairs]
    finish = [x[1] for x in sortpairs]
    files = [filt[x] for x in sortpairs]
    headers = [[d_header[f] for f in flist] for flist in files]

    rec = []
    for i in range(len(files)):
        names = [h['channel_names'][0] for h in headers[i]]
        d = {'names': names, 'start': start[i], 'finish': finish[i],
             'files': files[i], 'headers': headers[i]}
        rec.append(d)

    return rec


def split_rec_dir(nlx_dir, out_dir):
    """Split a directory with multiple recordings."""

    nlx_rec = read_rec_info(nlx_dir)
    for rec in nlx_rec:
        # make outer directory
        dirname = rec['start'].strftime('%Y-%m-%d_%H-%M-%S')
        dirpath = os.path.join(out_dir, dirname)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        # use heuristics to sort the different channel types
        rectype = []
        for header in rec['headers']:
            chan_name = header['channel_names'][0]
            if 'Spike' in chan_name or 'Micro' in chan_name:
                chan_type = 'micro'
            elif 'DC' in chan_name:
                chan_type = 'dc'
            else:
                chan_type = 'ieeg'
            rectype.append(chan_type)

        for i, src in enumerate(rec['files']):
            destdir = os.path.join(dirpath, rectype[i])
            if not os.path.exists(destdir):
                os.makedirs(destdir)

            filename = rec['names'][i] + os.path.splitext(src)[1]
            dest = os.path.join(destdir, filename)
            shutil.move(src, dest)


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
            sigs = data['ttl_input'][data['event_id'] == event_id]
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


def prep_nlx_ttl(bids_dir, sub, ses):
    """Prepare NLX TTLs for task alignment."""

    # get the session directory
    layout = BIDSLayout(bids_dir)
    runs = layout.get(subject=sub, session=ses)

    # read TTL signals from the NLX directory
    nlx_dir = os.path.join(runs[0].dirname, f'sub-{sub}_ses-{ses}_ieeg')
    recv_times, recv_signals = read_nlx_ttl(nlx_dir)
    data = pd.DataFrame({'onset': recv_times, 'signal': recv_signals})

    # write to a file for each run in the task (will be same for each;
    # there will be a sync solution for each run separately, correcting
    # for any clock drift across runs)
    for events_file in runs:
        recv_file = events_file.path.replace('_events.tsv', '_recv.tsv')
        data.to_csv(recv_file, sep='\t', index=False)


def plot_nlx_ttl(nlx_dir, out_file=None, interval=0.01):
    """Plot received sync signals."""

    # load file
    times, signals = read_nlx_ttl(nlx_dir)

    # translate events into continous signals
    sig_sync_times, sig_sync = sync.binary2analog(times / 10e5, signals, interval)

    # plot sync signal
    fig, ax = plt.subplots(figsize=(20, 4), dpi=300)
    ax.plot(sig_sync_times, sig_sync, linewidth=0.1)
    plt.tight_layout()

    # print to file
    if out_file is None:
        out_file = os.path.join(nlx_dir, 'Events.pdf')
    fig.savefig(out_file)
