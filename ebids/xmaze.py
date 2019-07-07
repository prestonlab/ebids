"""Module for exporting xmaze data to BIDS."""

import sys
import os
import glob
import warnings
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import json

def file_len(fname):
    """Number of lines in a file."""
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def log2events(log_file, events_file):
    """Convert an X-Maze log file into BIDS events."""

    # initialize different data outputs
    n_max = file_len(log_file)
    data = pd.DataFrame({'onset':np.zeros(n_max),
                         'duration':np.zeros(n_max),
                         'trial_type':np.empty(n_max, dtype=object),
                         'horiz_dist':np.zeros(n_max),
                         'diag_dist':np.zeros(n_max),
                         'pose':np.zeros(n_max),
                         'x_pos':np.zeros(n_max),
                         'z_pos':np.zeros(n_max),
                         'trial':np.zeros(n_max, dtype=int),
                         'onset_trial':np.zeros(n_max),
                         'chamber':np.zeros(n_max, dtype=int),
                         'reward':np.zeros(n_max, dtype=int),
                         'score':np.zeros(n_max, dtype=int)})
    
    # read all data lines
    n = 0
    trial = 0
    reward = 0
    score = 0
    with open(log_file, 'r') as log:
        for line in log:
            l = line.split()
            ttype = l[0]
            if ttype == 'Frame' and len(l) == 8 and l[1] != 'frameNum:':
                # just a frame logging position in the environment
                d = l[2:]
                data.at[n,'onset'] = d[3]
                data.at[n,'trial_type'] = 'frame'
                data.at[n,'horiz_dist'] = d[0]
                data.at[n,'diag_dist'] = d[1]
                data.at[n,'pose'] = d[2]
                data.at[n,'x_pos'] = d[4]
                data.at[n,'z_pos'] = d[5]
                data.at[n,'trial'] = trial
                data.at[n,'reward'] = reward
                data.at[n,'score'] = score

                n += 1
            elif ttype == 'Segment:' and len(l) == 10 and l[1] != 'distHori':
                # a specific event in the experiment
                d = l[1:]
                trial = d[7]
                
                data.at[n,'onset'] = d[3]
                data.at[n,'trial_type'] = d[8].lower()
                data.at[n,'horiz_dist'] = d[0]
                data.at[n,'diag_dist'] = d[1]
                data.at[n,'pose'] = d[2]
                data.at[n,'x_pos'] = d[4]
                data.at[n,'z_pos'] = d[5]
                data.at[n,'trial'] = trial
                data.at[n,'onset_trial'] = d[6]
                data.at[n,'reward'] = reward
                data.at[n,'score'] = score

                n += 1
            elif ttype == 'Selection' and len(l) == 5 and l[1] != 'trialNum:':
                # additional information about a selection; add to the
                # previous reward event and do not advance the trial counter
                d = l[2:]
                reward = d[1]
                score = d[2]
                data.at[n-1,'chamber'] = d[0]
                data.at[n-1,'reward'] = reward
                data.at[n-1,'score'] = score

    # remove unused trial space
    data = data[:n]
    data.to_csv(events_file, sep='\t', index=False)

def convert_sync(sync_file, out_file):
    """Convert a sync XML file to a BIDS-compatible table."""

    # read all sync start events
    tree = ET.parse(sync_file) 
    events = tree.findall('Event')
    sync = np.zeros(len(events))
    for i, event in enumerate(events):
        sync[i] = event.text

    # write out events to a text file
    data = pd.DataFrame({'sync':sync})
    data.to_csv(out_file, sep='\t', index=False, header=False)

    # write a sidecar file to identify the column
    base = os.path.basename(out_file).split('.')[0]
    sidecar = os.path.join(os.path.dirname(out_file), base + '.json')
    #with open(sidecar, 'w') as f:
    #    json.dump({'Columns':['onset']}, f)

def convert_session(raw_dir, bids_dir, sub, ses):
    """Convert raw data for a session to bids format."""

    if not os.path.exists(raw_dir):
        raise IOError('Raw directory does not exist: {}'.format(raw_dir))
    if not os.path.exists(bids_dir):
        raise IOError('BIDS directory does not exist: {}'.format(bids_dir))
    
    # create the standard session directory
    ses_dir = os.path.join(bids_dir, 'sub-{}'.format(sub),
                           'ses-{:02d}'.format(ses), 'func')
    if not os.path.exists(ses_dir):
        os.makedirs(ses_dir)

    n_run = 6
    for i in range(1, n_run+1):
        # find the log file for this run
        files = glob.glob(os.path.join(raw_dir, '*_task_run{}.xml'.format(i)))
        if len(files) > 1:
            warnings.warn('Found multiple log files for run {}. Skipping.'.format(i),
                          RuntimeWarning)
            continue
        elif not files:
            warnings.warn('Found no log files for run {}. Skipping.'.format(i),
                          RuntimeWarning)
            continue

        # convert the log file to a BIDS-style table
        log_file = files[0]
        run_name = 'sub-{}_ses-{:02d}_task-xmaze_run-{:02d}'.format(
            sub, ses, i)
        events_file = os.path.join(ses_dir, run_name + '_events.tsv')
        log2events(log_file, events_file)

        # convert the sync file to a regular text file, with sidecar
        raw_file = os.path.join(raw_dir, 'soundoutput_{}.xml'.format(i))
        if not os.path.exists(raw_file):
            warnings.warn('No sync file found for run {}.'.format(i),
                          RuntimeWarning)
            continue
        sync_file = os.path.join(ses_dir, run_name + '_sync.tsv')
        convert_sync(raw_file, sync_file)
