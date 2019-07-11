"""Module for processing Neuralynx data."""

import os
import shutil
import numpy as np
import neo.rawio.neuralynxrawio as nlxio

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
        except ValueError as err:
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
        d = {'names':names, 'start':start[i], 'finish':finish[i],
             'files':files[i], 'headers':headers[i]}
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
